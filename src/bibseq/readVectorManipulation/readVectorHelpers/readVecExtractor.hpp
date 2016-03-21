#pragma once
//
//  extraction.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/11/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "bibseq/alignment.h"
#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp"
#include "bibseq/seqToolsUtils/seqToolsUtils.hpp"
#include "bibseq/objects/helperObjects/motif.hpp"
namespace bibseq {

class readVecExtractor {

 public:
  template <typename T>
  static std::vector<T> extractReadsByPrimers(
       std::vector<T>& inputReads,
       std::map<std::string, std::map<std::string, std::vector<readObject>>>& readsByPrimer, std::map<std::string, std::pair<std::string, std::string>>& primersBySeq,
      uint32_t primerPos, bool lowerCasePrimer) {
    std::vector<T> outputReads;
    for (auto& rIter : inputReads) {
    	bool foundAPrimer = false;
    	for(const auto & primer : primersBySeq){
    		if(rIter.seqBase_.seq_.size() > (primerPos + primer.second.first.size())){
    			//check forward direction
    			if(primer.second.first == rIter.seqBase_.seq_.substr(primerPos, primer.second.first.size())){
    				if(lowerCasePrimer){
    					changeSubStrToLower(rIter.seqBase_.seq_, primerPos, primer.second.first.size());
    				}
    				readsByPrimer[primer.first]["forward"].emplace_back(rIter);
    				foundAPrimer = true;
    				break;
    			}
    			//check complement direction
    			if(primer.second.second == rIter.seqBase_.seq_.substr(primerPos, primer.second.second.size())){
    				if(lowerCasePrimer){
    					changeSubStrToLower(rIter.seqBase_.seq_, primerPos, primer.second.second.size());
    				}
    				readsByPrimer[primer.first]["reverse"].emplace_back(rIter);
    				foundAPrimer = true;
    				break;
    			}
    		}
    	}
    	if(!foundAPrimer){
    		outputReads.emplace_back(rIter);
    	}
    }
    return outputReads;
  }

	template <typename T>
  static std::vector<T> extractReadsByMID(
      const std::vector<T>& inputReads,
      std::map<std::string, std::vector<T>>& readsByMID, MapStrStr& midsByTags,
      int midPos, int midSize) {
    std::vector<T> outputReads;
    for (const auto& rIter : inputReads) {
      std::string currentMidSearch =
          rIter.seqBase_.seq_.substr(midPos, midSize);
      if (midsByTags.find(currentMidSearch) != midsByTags.end()) {
        if (readsByMID.find(midsByTags[currentMidSearch]) == readsByMID.end()) {
          readsByMID.insert(std::make_pair(midsByTags[currentMidSearch],
                                           std::vector<T>{rIter}));
        } else {
          readsByMID[midsByTags[currentMidSearch]].push_back(rIter);
        }
      } else {
        outputReads.push_back(rIter);
      }
    }
    return outputReads;
  }
	template <typename T>
  static std::vector<T> extractReadsByMIDVariableStart(
      const std::vector<T>& inputReads,
      std::map<std::string, std::vector<T>>& readsByMID, MapStrStr& midsByTags,
      uint32_t variableStart) {
		VecStr barcodesStr = getVectorOfMapKeys(midsByTags);
    std::vector<T> outputReads;
  	for(auto read: inputReads){
  		for(const auto & str : barcodesStr){
  			size_t pos = read.seqBase_.seq_.find(str);
  			if(pos < variableStart){
  				if(pos > 0){
  					read.trimFront(pos);
  				}
  				readsByMID[midsByTags[str]].emplace_back(read);
  				read.remove = true;
  				break;
  			}
  		}
  		if(!read.remove){
  			outputReads.emplace_back(read);
  		}
  	}
    return outputReads;
  }
	template <typename T>
  static std::vector<T> extractReadsByUnknownMID(
      const std::vector<T>& inputReads,
      std::map<std::string, std::vector<T>>& readsByMID,
      uint32_t midPos, uint32_t midSize, bool lowerBarcode) {
    std::vector<T> outputReads;
    for (const auto& rIter : inputReads) {
    	if(rIter.seqBase_.seq_.size() > midPos + midSize){
    		readsByMID[rIter.seqBase_.seq_.substr(midPos, midSize)].emplace_back(rIter);
    		changeSubStrToLower(readsByMID[rIter.seqBase_.seq_.substr(midPos, midSize)].back().seqBase_.seq_, midPos, midSize);
    	} else {
        outputReads.emplace_back(rIter);
      }
    }
    return outputReads;
  }

  template <typename T>
  static std::vector<T> extractReadsByMIDFull(
      const std::vector<T>& inputReads,
      std::map<std::string, std::vector<T>>& readsByMID,
      const std::string& midFilename, const std::string& midFileDelim,
      int midPos, bool checkComplement, bool keepBarcode, bool barcodeLower,
      bool variableStart, uint32_t variableStop, uint32_t errors, bool barcodeBothEnds) {

    std::vector<T> unrecognizedReads;
    int midSize = 10;

    // read in mids and barcodes
    MapStrStr midsByTags;
    table mids = seqUtil::readBarcodes(midFilename, midFileDelim, midSize);
    mids.outPutContentOrganized(std::cout);
    for (const auto& fIter : mids.content_) {
      midsByTags.insert(std::make_pair(fIter[1], fIter[0]));
    }
    if(errors > 0){
    	auto primers = seqUtil::readPrimers(midFilename, midFileDelim, false);
    	readObject forPrimer = readObject(seqInfo(primers.content_[0][0],primers.content_[0][1]));

    	unrecognizedReads = extractByMidsVariableErrors<T>(inputReads, readsByMID,
      		midsByTags, checkComplement, errors, variableStop, keepBarcode, midSize,
      		forPrimer);

    }else if( variableStart){
    	unrecognizedReads =
    	    	    	        extractReadsByMIDVariableStart(inputReads, readsByMID, midsByTags, variableStop);
    } else{
    	unrecognizedReads =
    	        extractReadsByMID(inputReads, readsByMID, midsByTags, midPos, midSize);
    }
    if (checkComplement && errors == 0) {
      // reverse complement and check again
      for (auto& read : unrecognizedReads) {
      	read.seqBase_.reverseComplementRead();
        read.seqBase_.name_.append("_Comp");
      }
      if(variableStart){
      	unrecognizedReads =
      	    	        extractReadsByMIDVariableStart(unrecognizedReads, readsByMID, midsByTags, variableStop);
      }else{
      	unrecognizedReads = extractReadsByMID(unrecognizedReads, readsByMID,
      	                                            midsByTags, midPos, midSize);
      }

      // reverse complement again to leave the unrecognized reads as they were
      // found
      for (auto& read : unrecognizedReads) {
        read.seqBase_.reverseComplementRead();
        read.seqBase_.name_ = replaceString(read.seqBase_.name_, "_Comp", "");
      }
    }

    if(barcodeBothEnds){
      // read in mids and barcodes
      MapStrStr midsByTagsRev;
      table midsRev = seqUtil::readBarcodes(midFilename, midFileDelim, midSize);
      //mids.outPutContentOrganized(std::cout);
      for (const auto& fIter : mids.content_) {
      	midsByTagsRev.insert(std::make_pair(seqUtil::reverseComplement(fIter[1], "DNA"), fIter[0]));
      }
      std::vector<T> revReadsUnRec;
      std::map<std::string, std::vector<T>> readsByMIDRev;
      if(errors > 0){
      	auto primers = seqUtil::readPrimers(midFilename, midFileDelim, false);
      	readObject revPrimer = readObject(seqInfo(primers.content_[0][0],primers.content_[0][2]));

      	revReadsUnRec = extractByMidsVariableErrors<T>(unrecognizedReads, readsByMIDRev,
        		midsByTags, checkComplement, errors, variableStop, keepBarcode, midSize,
						revPrimer);

      }else if( variableStart){
      	revReadsUnRec =
      	    	    	        extractReadsByMIDVariableStart(unrecognizedReads, readsByMIDRev, midsByTagsRev, variableStop);
      } else {
      	revReadsUnRec =
      	        extractReadsByMID(unrecognizedReads, readsByMIDRev, midsByTagsRev, midPos, midSize);
      }

      if (checkComplement && errors == 0) {
        // reverse complement and check again
        for (auto& read : revReadsUnRec) {
        	read.seqBase_.reverseComplementRead();
          read.seqBase_.name_.append("_Comp");
        }
        if(variableStart){
        	revReadsUnRec =
        	    	        extractReadsByMIDVariableStart(revReadsUnRec, readsByMIDRev, midsByTagsRev, variableStop);
        }else{
        	revReadsUnRec = extractReadsByMID(revReadsUnRec, readsByMIDRev,
        			midsByTagsRev, midPos, midSize);
        }
        // reverse complement again to leave the unrecognized reads as they were
        // found
        for (auto& read : revReadsUnRec) {
          read.seqBase_.reverseComplementRead();
          read.seqBase_.name_ = replaceString(read.seqBase_.name_, "_Comp", "");
        }
      }
      if(errors < 1){
        if (keepBarcode) {
          if (barcodeLower) {
            for (auto& rIter : readsByMIDRev) {
              readVec::changeFrontEndToLowerCase(rIter.second, midSize -1);
            }
          }
        } else {
          for (auto& rIter : readsByMIDRev) {
            readVecTrimmer::trimOffForwardBases(rIter.second, midSize);
          }
        }
      }
      if (keepBarcode) {
        if (barcodeLower) {
          for (auto& rIter : readsByMID) {
            readVec::changeBackEndToLowerCase(rIter.second, midSize);
          }
        }
      } else {
        for (auto& rIter : readsByMID) {
          readVecTrimmer::trimOffEndBases(rIter.second, midSize);
        }
      }
      //now insert the reads in the correct orientation becuase they were found with the reverse primer
      //and add to the out reads
      for (auto& rIter : readsByMIDRev) {
      	readVec::allReverseComplement(rIter.second, true);
      	if(errors > 0){
          if (keepBarcode) {
            if (barcodeLower) {
              for (auto& rIter : readsByMID) {
                readVec::changeFrontEndToLowerCase(rIter.second, midSize -1);
              }
            }
          } else {
            for (auto& rIter : readsByMID) {
              readVecTrimmer::trimOffForwardBases(rIter.second, midSize);
            }
          }
      	}
      	addOtherVec(readsByMID[rIter.first],rIter.second);
      }

      unrecognizedReads = revReadsUnRec;
    }
    if(errors < 1){
      if (keepBarcode) {
        if (barcodeLower) {
          for (auto& rIter : readsByMID) {
            readVec::changeFrontEndToLowerCase(rIter.second, midSize -1);
          }
        }
      } else {
        for (auto& rIter : readsByMID) {
          readVecTrimmer::trimOffForwardBases(rIter.second, midSize);
        }
      }
    }
    return unrecognizedReads;
  }



  template <typename T>
  static std::vector<T> extractSameReadsByName(
      const std::vector<T> inputReads, const std::vector<T>& compareReads) {
    // untested nick, 8/7/13
    std::vector<T> ans;
    MapStrStr names;
    for (const auto& compIter : compareReads) {
      names.insert({compIter.seqBase_.name_, compIter.seqBase_.name_});
    }
    for (const auto inputIter : inputReads) {
      if (names.find(inputIter.seqBase_.name_) != names.end()) {
        ans.push_back(inputIter);
      }
    }
    return ans;
  }
  template <typename T>
  static std::vector<T> extractReadsWithNames(const std::vector<T>& reads,
                                              const VecStr& names) {
    std::vector<T> ans;
    for (const auto& read : reads) {
      if (contains(names, read.seqBase_.name_)) {
        ans.emplace_back(read);
      }
    }
    return ans;
  }
  template <typename T>
  static std::vector<T> extractReadsExcludingNames(const std::vector<T>& reads,
                                              const VecStr& names) {
    std::vector<T> ans;
    for (const auto& read : reads) {
      if (!contains(names, read.seqBase_.name_)) {
        ans.emplace_back(read);
      }
    }
    return ans;
  }

  template<typename T>
  static std::vector<T> extractByMidsVariableErrors(std::vector<T> reads,
  		std::map<std::string, std::vector<T>>& readsByMID, MapStrStr midsByTags,
  		bool checkComplement,uint32_t maxErrors, uint64_t stop, bool keepBarcode,
  		uint32_t barcodeSize,
  		readObject forwardPrimer){
  	//std::cout <<"here?" << "\n";
  	std::unordered_map<std::string, motif> barcodeMotifs;
  	for(const auto & bar : midsByTags){
  		barcodeMotifs.emplace(bar.second, motif(bar.first));
  	}
  	bool debug = false;
  	std::vector<T> unmatched;
  	for(auto & read : reads){
  		if(read.seqBase_.seq_.size() < stop + barcodeSize){
  			unmatched.emplace_back(read);
  			continue;
  		}
  		bool found = false;
  		for(const auto & bar : barcodeMotifs){
  			auto positions = bar.second.findPositionsFull(read.seqBase_.seq_, 0, 0, stop);
  			if(!positions.empty()){
  				if(keepBarcode){
  					read.trimFront(positions.back());
  				}else{
  					read.trimFront(positions.back() + bar.second.motifOriginal_.size());
  				}
  				readsByMID[bar.first].emplace_back(read);
  				found = true;
  				break;
  			}
  		}
  		if(checkComplement && !found){
  			T readComp = read;//(read.seqBase_);
  			readComp.seqBase_.reverseComplementRead(true);
  			for(const auto & bar : barcodeMotifs){
  				auto positions = bar.second.findPositionsFull(readComp.seqBase_.seq_, 0, 0, stop);
  				if(!positions.empty()){
  					//printVector(positions, ", ", std::cout);
  					//std::cout << "\t" << readComp.seqBase_.seq_.size() << std::endl;;
  					if(keepBarcode){
  						readComp.trimFront(positions.back());
  					}else{
  						readComp.trimFront(positions.back() + bar.second.motifOriginal_.size());
  					}
  					readsByMID[bar.first].emplace_back(readComp);
  					found = true;
  					break;
  				}
  			}
  		}
  		if(!found){
  			unmatched.emplace_back(read);
  		}
  	}
  	std::unordered_map<uint32_t, uint32_t> multimatched;
  	std::unordered_map<uint32_t, uint32_t> errorMatched;
  	if(maxErrors > 1){
  		uint64_t maxReadLen = 0;
  		readVec::getMaxLength(reads, maxReadLen);
  		substituteMatrix scoring(2,-2);
  		scoring.setWithDegenScoringCaseSense(2,-2);
  		aligner alignerObj(maxReadLen, gapScoringParameters(7,1,7,1,7,1), scoring);
  		std::vector<uint64_t> eraseThese;
  		//std::cout << "in the errors" << std::endl;
  		for(auto readPos : iter::range(unmatched.size())){
  			if(unmatched[readPos].seqBase_.seq_.size() < stop + barcodeSize){
  				continue;
  			}
  			//std::cout << readPos << ":" << unmatched.size()<< "\n";
  			if(readPos % 100 == 0){
  				std::cout << readPos << ":" << unmatched.size()<< "\n";
  			}

  			auto & read = unmatched[readPos];
  			//read.printDescription(std::cout, true);
  			auto readComp = unmatched[readPos];
  			readComp.seqBase_.reverseComplementRead(true);
  			//readComp.printDescription(std::cout, true);
  			auto forwardPosition = alignerObj.findReversePrimer(read.seqBase_.seq_,forwardPrimer.seqBase_.seq_);
  			uint32_t forStop = 0;
  			if(forwardPosition.first > 20){
  				forStop = forwardPosition.first - 20;
  			}
  			if(forwardPosition.first > stop && !checkComplement){
  				if(debug){
  					std::cout << "forwardPosition.first > stop" << std::endl;
  					std::cout << forwardPosition.first << std::endl;
  				}
  				//if can find the forward primer within the search distance break out but wait for the reverse check
  				continue;
  			}
  			std::pair<uint32_t,uint32_t> forwardPosComp;
  			uint32_t forStopComp = 0;
  			if(checkComplement){
  				forwardPosComp = alignerObj.findReversePrimer(readComp.seqBase_.seq_, forwardPrimer.seqBase_.seq_);
  				if(forwardPosComp.first > 20){
  					forStopComp = forwardPosComp.first - 20;
  				}
  				if(forwardPosition.first > stop && forwardPosComp.first > stop){
  					if(debug){
  						std::cout << "forwardPosition.first > stop && forwardPosComp.first > stop" << std::endl;
  						std::cout << forwardPosition.first << std::endl;
  						std::cout << forwardPosComp.first << std::endl;
  					}
  					//if can find the forward primer within the search distance break out
  					continue;
  				}
  			}

  			bool found = false;
  			bool foundComp = false;
  			std::string appendStr = "";
  			std::string matchBar = "";
  			if(debug){
  				std::cout << forStop << std::endl;
  				std::cout << forwardPosition.first << std::endl;
  				std::cout << forStopComp << std::endl;
  				std::cout << forwardPosComp.first << std::endl;
  			}
  			for(auto err : iter::range<uint32_t>(1,maxErrors)){
  				for(const auto & bar : barcodeMotifs){
  					auto positions = bar.second.findPositionsFull(read.seqBase_.seq_, err, forStop, forwardPosition.first);
  					if(!positions.empty()){
  						appendStr = "_erNum:" + to_string(err);
  						matchBar = bar.first;
  						if(found){
  							//found multiple possible matches
  							found = false;
  							read.seqBase_.name_.append("_multiMatched");
  							++multimatched[err];
  							break;
  						}
  						found = true;
  					}else if (checkComplement){
  						auto positionsComp = bar.second.findPositionsFull(readComp.seqBase_.seq_, err, forStopComp, forwardPosComp.first);
  						if(!positionsComp.empty()){
  							appendStr = "_erNum:" + to_string(err) + "_Comp";
  							matchBar = bar.first;
  							if(found){
  								//found multiple possible matches
  								read.seqBase_.name_.append("_multiMatched");
  								found = false;
  								++multimatched[err];
  								break;
  							}
  							found = true;
  							foundComp = true;
  						}
  					}
  				}
  				if(found){
  					++errorMatched[err];
  					if(foundComp){
  						auto positions = barcodeMotifs.at(matchBar).findPositionsFull(readComp.seqBase_.seq_, err, forStopComp, forwardPosComp.first);
  						if(keepBarcode){
  							//printVector(positions, ", ", std::cout);
  							readComp.trimFront(positions.back());
  						}else{
  							readComp.trimFront(positions.back() + barcodeMotifs.at(matchBar).motifOriginal_.size());
  						}
  						readComp.seqBase_.name_ = replaceString(readComp.seqBase_.name_, "_Comp", appendStr);
  						readsByMID[matchBar].emplace_back(readComp);
  					}else{
  						auto positions = barcodeMotifs.at(matchBar).findPositionsFull(read.seqBase_.seq_, err, forStop, forwardPosition.first);
  						if(keepBarcode){
  							read.trimFront(positions.back());
  						}else{
  							read.trimFront(positions.back() + barcodeMotifs.at(matchBar).motifOriginal_.size());
  						}
  						read.seqBase_.name_.append(appendStr);
  						readsByMID[matchBar].emplace_back(read);
  					}
  					eraseThese.emplace_back(readPos);
  					break;
  				}
  			}
  		}
  		for(const auto & rem : iter::reversed(eraseThese)){
  			unmatched.erase(unmatched.begin() + rem);
  		}
  	}
  	return unmatched;
  }
};
}  // namespace bib
