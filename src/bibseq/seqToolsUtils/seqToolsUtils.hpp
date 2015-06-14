#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//
//  seqToolsUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 4/27/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/seqObjects/sffObject.hpp"

#include "bibseq/alignment.h"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp"
#include "bibseq/simulation/simulationCommon.hpp"

#include "bibseq/simulation/errorProfile.hpp"
#include "bibseq/simulation/profile.hpp"
#include "bibseq/objects/helperObjects/probabilityProfile.hpp"

#include <bibcpp/graphics/colorUtils.hpp>


namespace bibseq {




void getReadCnt(const std::string & filename, std::string & format,
		bool processed, uint64_t & currentCount) ;

void setUpSampleDirs(
    const std::string& barcodeFilename,
		const std::string& mainDirectoryName) ;
std::pair<uint32_t,uint32_t> getMaximumRounds(double cnt);
std::vector<char> getAAs(bool noStopCodon);
class cluster;

template<typename T>
simulation::profile simProfileWithReads(const std::vector<T> & reads, aligner & alignerObj,
		randomGenerator & gen, uint32_t sampNum, bool local, bool ignoreGaps){
  std::vector<uint32_t> randomFirstPos = gen.unifRandVector<uint32_t>(0, len(reads),sampNum);
  std::vector<uint32_t> randomSecondPos = gen.unifRandVector<uint32_t>(0, len(reads),sampNum);
  simulation::profile allProfile;
  for(const auto & pos : iter::range(sampNum)){
  	alignerObj.alignVec(reads[randomFirstPos[pos]],
  			reads[randomSecondPos[pos]], local);
  	allProfile.increaseCountAmount(alignerObj.alignObjectA_.seqBase_.seq_,
  	  				alignerObj.alignObjectB_.seqBase_.seq_, 1);
  }
  allProfile.setProfile(true);
  return allProfile;
}
simulation::errorProfile getErrors(const simulation::profile & allProfile,
		const std::vector<char> & alphabet);
substituteMatrix getMatrixFromProfile(simulation::profile allProfile,
		const std::vector<char> & alphabet);

std::string genHtmlStrForPsuedoMintree(std::string jsonFileName);




template <typename T>
void findForwardPrimerAndSort(
    std::vector<T>& inputReads,
    const std::map<std::string, std::pair<std::string, std::string>>& primers,
    std::map<std::string, std::vector<T>>& reads,
    std::vector<T>& unrecognizedPrimer, bool condensed) {
  size_t maxSize = 0;
  for (const auto& currentPrimer : primers) {
    if (currentPrimer.second.first.size() > maxSize) {
      maxSize = currentPrimer.second.first.size();
    }
  }
  for (auto& read : inputReads) {
    if (read.seqBase_.seq_.size() < maxSize) {
      unrecognizedPrimer.push_back(read);
    } else {

      bool foundMatch = false;
      for (const auto& currentPrimer : primers) {
        if (currentPrimer.first == "barcodeSize") {
          continue;
        }
        // find forward primer or if it isn't found to the unrecognized primer
        // vector you go
        std::string primerComparer;
        std::string readComparer;
        if (condensed) {
          primerComparer = condenseSeqSimple(currentPrimer.second.first);
          read.createCondensedSeq();
          readComparer = read.condensedSeq.substr(0, primerComparer.size());
        } else {
          primerComparer = currentPrimer.second.first;
          readComparer =
              read.seqBase_.seq_.substr(0, currentPrimer.second.first.size());
        }
        if (primerComparer == readComparer) {
          changeSubStrToLowerFromBegining(read.seqBase_.seq_,
                                          currentPrimer.second.first.length());
          if (reads.find(currentPrimer.first) == reads.end()) {
            reads.insert({currentPrimer.first, {read}});
          } else {
            reads[currentPrimer.first].push_back(read);
          }
          foundMatch = true;
          break;
        }
      }
      if (!foundMatch) {
        unrecognizedPrimer.push_back(read);
      }
    }
  }
}

template <typename T>
void findForwardPrimerAndSortDegenerative(
    std::vector<T>& inputReads,
    const std::map<std::string, std::pair<std::string, std::string>>& primers,
    std::map<std::string, std::vector<T>>& reads,
    std::vector<T>& unrecognizedPrimer,uint32_t variableStart,
    double queryCoverage,
    uint32_t numOfMismatches, double percentageGaps) {
  uint64_t maxSize = 0;
  for (const auto& currentPrimer : primers) {
    if (currentPrimer.second.first.size() > maxSize) {
      maxSize = currentPrimer.second.first.size();
    }
  }
  uint64_t maxReadLength = 0;
  readVec::getMaxLength(inputReads, maxReadLength);
  auto scoringMatrixMap = substituteMatrix::createDegenScoreMatrix(1,-1);

  // create alinger class object
  auto gapPars = gapScoringParameters (7, 1, 7, 1, 7, 1);
  aligner alignerObj =
      aligner(maxReadLength, gapPars, scoringMatrixMap);

  int count = 1;
  std::cout << "Finding forward primer with degenerative bases in mind"
            << std::endl;
  for (auto& read : inputReads) {
    if (count % 5000 == 0) {
      std::cout << "Curently on " << count << " of " << inputReads.size()
                << std::endl;
    }
    ++count;
    if (read.seqBase_.seq_.size() < maxSize) {
      unrecognizedPrimer.push_back(read);
    } else {
      bool foundMatch = false;
      for (const auto& currentPrimer : primers) {
        if (currentPrimer.first == "barcodeSize") {
          continue;
        }
        // find forward primer or if it isn't found to the unrecognized primer
        // vector you go
        auto primer =
            baseReadObject(seqInfo("primer", currentPrimer.second.first));
        auto readBegin =
                    baseReadObject(seqInfo("readBegin", read.seqBase_.seq_.substr(0, variableStart + maxSize)));
        /*auto forwardPosition = alignerObj.findReversePrimer(read.seqBase_.seq_, primer.seqBase_.seq_);
        alignerObj.rearrangeObjs(read.seqBase_, primer.seqBase_, true);
        alignerObj.profilePrimerAlignment(read, primer, false);*/
        auto forwardPosition = alignerObj.findReversePrimer(readBegin.seqBase_.seq_, primer.seqBase_.seq_);
				alignerObj.rearrangeObjs(readBegin, primer.seqBase_, true);
				alignerObj.profilePrimerAlignment(readBegin, primer, false);
        if (forwardPosition.first <= variableStart
            && alignerObj.comp_.distances_.queryCoverage_ >= queryCoverage
            && alignerObj.comp_.distances_.percentageGaps_ <= percentageGaps
            && alignerObj.comp_.hqMismatches_ <= numOfMismatches) {
          changeSubStrToLowerFromBegining(read.seqBase_.seq_,
                                          currentPrimer.second.first.length());
          reads[currentPrimer.first].emplace_back(read);
          foundMatch = true;
          break;
        }
      }
      if (!foundMatch) {
        unrecognizedPrimer.push_back(read);
      }
    }
  }
}

template <typename T>
void checkForReversePrimer(T& read, readObject& reversePrimer,
                           runningParameters runParmas, aligner& alignObj,
                           bool reverserPrimerToLowerCase,
                           bool weighHomopolyers, double queryCoverageCutOff,
                           double percentIdentityCutoff) {
  auto rPos = alignObj.findReversePrimer(read, reversePrimer);
  alignObj.rearrangeObjs(read.seqBase_, reversePrimer.seqBase_, true);
  alignObj.profilePrimerAlignment(read, reversePrimer, weighHomopolyers);
  bool primerGood = true;
  if (alignObj.comp_.largeBaseIndel_ > runParmas.errors_.largeBaseIndel_ ||
      alignObj.comp_.oneBaseIndel_ > runParmas.errors_.oneBaseIndel_ ||
      alignObj.comp_.twoBaseIndel_ > runParmas.errors_.twoBaseIndel_) {
    primerGood = false;
  }
  if (alignObj.comp_.distances_.percentIdentity_ < percentIdentityCutoff) {
    primerGood = false;
  }
  if (alignObj.comp_.distances_.queryCoverage_ < queryCoverageCutOff) {
    primerGood = false;
  }

  read.remove = !primerGood;
  read.seqBase_.on_ = primerGood;
  if (primerGood) {
    if (reverserPrimerToLowerCase) {
      while (read.seqBase_.seq_[rPos.first] ==
             read.seqBase_.seq_[rPos.first - 1]) {
        --rPos.first;
      }
      changeSubStrToLowerToEnd(read.seqBase_.seq_, rPos.first);
    }
    read.setClip(0, rPos.second);
  }
}

template<typename T>
void checkForReversePrimer(std::vector<T>& reads, readObject& reversePrimer,
		runningParameters runParmas, aligner& alignObj,
		bool reverserPrimerToLowerCase, bool weighHomopolyers,
		double queryCoverageCutOff, double percentIdentityCutoff) {
	for (auto& read : reads) {
		checkForReversePrimer(read, reversePrimer, runParmas, alignObj,
				reverserPrimerToLowerCase, weighHomopolyers, queryCoverageCutOff,
				percentIdentityCutoff);
	}
}





template <typename T>
void sortReadsBySeqFront(const std::vector<T>& inputReads,
                         std::map<std::string, std::vector<readObject>>& reads,
                         size_t compareSize, std::vector<T>& unrecognizedPrimer,
                         size_t cutOffMinSize) {
  std::map<std::string, std::vector<readObject>> tempReads;

  for (const auto& read : inputReads) {
    if (read.seqBase_.seq_.size() < compareSize) {
      unrecognizedPrimer.push_back(read);
    } else {
      if (tempReads.find(read.seqBase_.seq_.substr(0, compareSize)) ==
          tempReads.end()) {
        tempReads.insert({read.seqBase_.seq_.substr(0, compareSize), {read}});
      } else {
        tempReads[read.seqBase_.seq_.substr(0, compareSize)].push_back(read);
      }
    }
  }
  for (const auto& read : tempReads) {
    if (read.second.size() > cutOffMinSize) {
      reads.insert(read);
    } else {
      addOtherVec(unrecognizedPrimer, read.second);
    }
  }
}

template <typename T>
void renameReadNames(std::vector<T>& reads, const std::string& stub,
                     bool processed, bool keepChimeraFlag = true,
                     bool keepCompFlag = true,
                     const std::string& sortBy = "none") {
	uint64_t count = 0;
  if (sortBy != "none") {
    readVecSorter::sortReadVector(reads, sortBy);
  }
  uint64_t maxSize = reads.size();
  for (auto& read : reads) {
    bool chimera = read.seqBase_.name_.find("CHI") != std::string::npos;
    bool comp = read.seqBase_.name_.find("_Comp") != std::string::npos;
    read.seqBase_.name_ = stub;
    if (chimera && keepChimeraFlag) {
      read.seqBase_.markAsChimeric();
    }
    if (comp && keepCompFlag) {
    	read.seqBase_.name_ += "_Comp";
    }
    read.seqBase_.name_ += "." + leftPadNumStr(count, maxSize);
    if (processed) {
      read.seqBase_.name_ += "_t" + to_string(read.seqBase_.cnt_);
    }
    ++count;
  }
}

template <typename T>
void renameReadNamesNewClusters(std::vector<T>& reads, const std::string& stub,
                                bool processed, bool keepChimeraFlag = true,
                                bool keepCompFlag = true,
                                const std::string& sortBy = "none") {
  renameReadNames(reads, stub, processed, keepChimeraFlag, keepCompFlag, sortBy);
  for (auto& read : reads) {
    read.reads_[0].seqBase_.name_ = read.seqBase_.name_;
  }
}



void processRunCutoff(uint32_t& runCutOff, const std::string& runCutOffString, int counter);

uint32_t processRunCutoff(const std::string& runCutOffString, uint64_t counter);





/////////tools for finding additional location output
std::string makeIDNameComparable(const std::string& idName);

std::string findIdNameFromFileName(const std::string& filename);

std::string processFileNameForID(const std::string& fileName);

std::string findAdditonalOutLocation(const std::string& locationFile,
                                     const std::string& fileName);




void processAlnInfoInput(aligner& alignerObj,
                         const std::string& alnInfoDirName);
template <typename T>
uint32_t seqQualSizeAgreementCheck(const std::vector<T>& reads) {
  uint32_t count = 0;
  for (const auto& read : reads) {
    if (read.seqBase_.seq_.size() != read.seqBase_.qual_.size()) {
      std::cout << "name:" << read.seqBase_.name_ << std::endl;
      std::cout << "seq:" << read.seqBase_.seq_ << std::endl;
      std::cout << "qual:" << read.seqBase_.qual_ << std::endl;
      std::cout << "sSize:" << read.seqBase_.seq_.size() << std::endl;
      std::cout << "qSize:" << read.seqBase_.qual_.size() << std::endl;
      ++count;
    }
  }
  std::cout << getPercentageString(count, reads.size()) << std::endl;
  return count;
}
template <typename T>
std::unordered_map<std::string, bib::color> getColorsForNames(
    const std::vector<T>& reads, double sat, double lum) {
  std::unordered_map<std::string, bib::color> colorsForName;
  uint32_t count = 0;
  std::vector<bib::color> colors = bib::evenHuesAll(sat, lum, len(reads));
  for (const auto& read : reads) {
    colorsForName[read.getReadId()] = colors[count];
    ++count;
  }
  return colorsForName;
}




std::vector<uint32_t> getWindowQuals(const std::vector<uint32_t>& qual,
                                     uint32_t qWindowSize, uint32_t pos);
std::vector<double> likelihoodForBaseQ(
    const std::vector<uint32_t>& qual,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMeanQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMedianQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMinQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);


cluster createClusterFromSimiliarReads(std::vector<readObject> similiarReads,
                                       aligner& alignerObj);

template<typename READ>
std::vector<READ> vecStrToReadObjs(const VecStr & strs, const std::string & stubName){
	std::vector<READ> ans;
	for(const auto & strPos : iter::range(strs.size())){
		ans.emplace_back(READ(seqInfo(stubName + "." + leftPadNumStr(strPos, strs.size()), strs[strPos])));
	}
	return ans;
}

template<typename T>
VecStr readObjsToVecStr(const std::vector<T> & vec){
	VecStr ans;
	for(const auto & read : vec){
		ans.emplace_back(read.seqBase_.seq_);
	}
	return ans;
}



template<typename READ, typename REF>
bool checkPossibleChiByRefs(const READ & read, const std::vector<REF> & refSeqs, table& outInfo, aligner & alignerObj,
		const comparison & chiOverlap, bool breakAtFirst,   bool weightHomopolymers ){
	bool foundAnExactMatch = false;
	bool foundAChimera = false;
	std::string firstChiName = "";
	std::string secondChiName = "";
	uint32_t inflectionPoint = UINT32_MAX;
	uint32_t inflectionPointPar1 = UINT32_MAX;
	uint32_t inflectionPointPar2 = UINT32_MAX;
	for(const auto & refPos : iter::range(len(refSeqs))){
		auto & ref = refSeqs[refPos];
		alignerObj.alignVec(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_, weightHomopolymers);
		if(alignerObj.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
			return foundAChimera;
		}
	}
	for(const auto & refPos : iter::range(len(refSeqs))){
		auto & ref = refSeqs[refPos];
		alignerObj.alignVec(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_, weightHomopolymers);
		if(alignerObj.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
		} else if(alignerObj.mismatches_.size() > 0 && !foundAnExactMatch){
			auto savedMismatches = alignerObj.mismatches_;
			//check front
			alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, 11, true, false, true, false,
					weightHomopolymers, 0,
					alignerObj.getAlignPosForSeqBPos(savedMismatches.begin()->second.seqBasePos));
			bool passFront = chiOverlap.passErrorProfile(alignerObj.comp_);

			//check back
			alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, 11, true, false, true, false,
					weightHomopolymers,
								alignerObj.getAlignPosForSeqBPos(savedMismatches.rbegin()->second.seqBasePos + 1));
			bool passBack = chiOverlap.passErrorProfile(alignerObj.comp_);

			auto firstRefAlignA = alignerObj.alignObjectA_;
			auto firstRefAlignB = alignerObj.alignObjectB_;
			//std::cout << "pass front: " << passFront << std::endl;
			//std::cout << "pass back: " << passBack << std::endl;
			if(passFront){
				for(const auto & secondRefPos : iter::range(refPos + 1, len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}

					alignerObj.alignVec(secondRef.seqBase_, read.seqBase_, false);
					//check to see if from the mismatch on from the other one ref matches
					//to another ref
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers);
					if(alignerObj.mismatches_.size() < 1){
						continue;
					}
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers,
												alignerObj.getAlignPosForSeqBPos(savedMismatches.begin()->second.seqBasePos));

					//alignerObj.errors_.printDescription(std::cout, true);
					if(chiOverlap.passErrorProfile(alignerObj.comp_)){
						firstChiName = ref.seqBase_.name_;
						secondChiName = secondRef.seqBase_.name_;
						inflectionPoint = savedMismatches.begin()->second.seqBasePos;
						inflectionPointPar1 = savedMismatches.begin()->second.refBasePos;
						std::string refTwoPortion = alignerObj.alignObjectA_.seqBase_.seq_.substr(0, alignerObj.getAlignPosForSeqAPos(savedMismatches.begin()->second.refBasePos) + 1);
						inflectionPointPar2 = removeCharReturn(refTwoPortion, '-').size() - 1;

						foundAChimera = true;
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, to_string(read.seqBase_.frac_),
							to_string(inflectionPoint),
							firstChiName, "1", to_string(inflectionPointPar1),
							secondChiName, "1", to_string(inflectionPointPar2)});
						if(foundAChimera && breakAtFirst){
							//VecStr{"readName", "fraction", "par1", "par1Frac", "par2", "par2Frac", "inflectionPoint"}
							break;
						}
					}
				}
			}
			if(foundAChimera && breakAtFirst){
				break;
			}
			if(passBack){
				for(const auto & secondRefPos : iter::range(refPos + 1,len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}
					alignerObj.alignVec(secondRef.seqBase_, read.seqBase_, false);
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers);
					if(alignerObj.mismatches_.size() < 1){
						continue;
					}
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers,
												0, alignerObj.getAlignPosForSeqBPos(savedMismatches.rbegin()->second.seqBasePos + 1));
					if(chiOverlap.passErrorProfile(alignerObj.comp_)){
						foundAChimera = true;
						firstChiName = ref.seqBase_.name_;
						secondChiName = secondRef.seqBase_.name_;
						inflectionPoint = savedMismatches.rbegin()->second.seqBasePos;
						inflectionPointPar1 = savedMismatches.rbegin()->second.refBasePos;
						std::string refTwoPortion = alignerObj.alignObjectA_.seqBase_.seq_.substr(0, alignerObj.getAlignPosForSeqAPos(savedMismatches.rbegin()->second.refBasePos) + 1);
						inflectionPointPar2 = removeCharReturn(refTwoPortion, '-').size() - 1;
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, to_string(read.seqBase_.frac_),
							to_string(inflectionPoint),
							firstChiName, "1", to_string(inflectionPointPar1),
							secondChiName, "1", to_string(inflectionPointPar2)});
						if(foundAChimera && breakAtFirst){
							break;
						}
					}
				}
			}
		}
	}
	return foundAChimera;
}



uint32_t processCutOffStr(const std::string& runCutOffString,
  uint64_t readCount);

uint32_t getAlnPos(const std::string & seq, uint32_t realSeqPos);

uint32_t getRealPos(const std::string & seq, uint32_t seqAlnPos);


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "seqToolsUtils.cpp"
#endif
