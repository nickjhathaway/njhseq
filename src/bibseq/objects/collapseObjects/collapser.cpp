//
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
//  collapser.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/1/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "collapser.hpp"
#include "bibseq/objects/seqObjects/seqWithKmerInfo.hpp"
#include "bibseq/objects/helperObjects/nucCompCluster.hpp"
#include "bibseq/helpers/profiler.hpp"

namespace bibseq {

void collapser::runFullClustering(std::vector<cluster> & clusters,
		bool onPerId, std::map<int, std::vector<double>> iteratorMap,
		bool useNucComp,bool useMinLenNucComp, bool findBestNuc, const std::vector<double> & diffCutOffVec,
		bool useKmerBinning, uint32_t kCompareLen, double kmerCutOff,
		aligner & alignerObj, seqSetUp & setUp,
		bool snapShots, const std::string & snapShotsDirName){
  clusterVec::allSetFractionClusters(clusters);
  uint32_t minLen = readVec::getMinLength(clusters);
  if (useNucComp) {
  	//get the alphabet of the current cluster set and use that in the nuc comp setting
  	charCounterArray allCounter;
  	for(const auto & read : clusters){
  		allCounter.increaseCountByString(read.seqBase_.seq_, read.seqBase_.cnt_);
  	}
  	allCounter.resetAlphabet(false);
  	auto alphabet = allCounter.alphabet_;
  	//set the verbosity to debug level
  	bool oldVerbose = opts_.verbose_;
  	opts_.verbose_ = setUp.debug_;
  	//now cluster on nucleotide composition clusters
  	for(const auto & diffCutOff : diffCutOffVec){
    	std::vector<nucCompCluster> comps;
    	if(useMinLenNucComp){
    		comps = clusterOnNucComp(clusters, minLen, alphabet, diffCutOff, findBestNuc, true, setUp.debug_);
    	}else{
    		comps = clusterOnNucComp(clusters, alphabet, diffCutOff, findBestNuc, true, setUp.debug_);
    	}
			if (oldVerbose) {
				std::cout << "On Nucleotide Composition Bin difference of "
						<< diffCutOff << " in " << vectorToString(diffCutOffVec, ", ")
						<< std::endl;
			}
			//now run clustering on the nucleotide clusters
    	uint32_t compCount = 0;
    	for(auto & comp : comps){
    		++compCount;
    		if(oldVerbose){
      		std::cout << "On " << compCount << " of " << comps.size() << "\r";
      		std::cout.flush();
    		}
    		if(comp.readPositions_.size() < 2){
    			continue;
    		}
    		double currentReadCnt = 0;
    		for(const auto & pos : comp.readPositions_){
    			currentReadCnt += clusters[pos].seqBase_.cnt_;
    		}
    		//set the run cut off of low kmer frequency to the current count
    		auto oldRunCutoff = alignerObj.kMaps_.runCutOff_ ;
    		uint32_t currentRunCutoff = processCutOffStr(setUp.runCutOffString_,currentReadCnt);
    		alignerObj.kMaps_.runCutOff_ = currentRunCutoff;
    		/**@todo need to sort here*/
    		if(onPerId){
    			runClusteringOnId(clusters, comp.readPositions_, iteratorMap, alignerObj);
    		}else{
    			runClustering(clusters, comp.readPositions_, iteratorMap, alignerObj);
    		}
    		alignerObj.kMaps_.runCutOff_ = oldRunCutoff;
    	}
    	if(oldVerbose){
    		std::cout << std::endl;
    	}
  	}
  	opts_.verbose_ = oldVerbose;
  	//remove all the clusters marked removed
  	clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
  } else if (useKmerBinning){
  	//set verbosity at debug level
  	bool oldVerbose = opts_.verbose_;
  	opts_.verbose_ = setUp.debug_;
  	//create vector of seqWtihKmerInfo reads to use to clusters by kmer distance
  	std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  	reads.reserve(clusters.size());
  	for(const auto & read : clusters){
  		reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  	}
  	allSetKmers(reads, kCompareLen, false);
  	std::vector<kmerClusterPos> kClusters;
  	//now cluster reads based on kmers
  	for(const auto & readPos : iter::range(reads.size())){
  		if(setUp.debug_ && readPos % 50 == 0){
  			std::cout << "currently on " << readPos << " of "
  					<< reads.size() << std::endl;
  		}
  		bool foundMatch = false;
  		for(auto & kClus : kClusters){
  			foundMatch = kClus.compareRead(reads[readPos], readPos, kmerCutOff,
  					false);
  			if(foundMatch){
  				break;
  			}
  		}
  		if(!foundMatch){
  			kClusters.emplace_back(kmerClusterPos(reads[readPos], readPos));
  		}
  	}
  	//run clustering on the different kmer clustering
  	uint32_t clusterCount = 0;
  	for(auto & kClus : kClusters){
  		++clusterCount;
  		if(clusterCount % 50 == 0){
  			std::cout << "On " << clusterCount << " of " << kClusters.size() << std::endl;
  		}
  		if(kClus.readPositions_.size() < 2){
  			continue;
  		}
  		double currentReadCnt = 0;
  		for(const auto & pos : kClus.readPositions_){
  			currentReadCnt += clusters[pos].seqBase_.cnt_;
  		}
  		auto oldRunCutoff = alignerObj.kMaps_.runCutOff_ ;
  		uint32_t currentRunCutoff = processCutOffStr(setUp.runCutOffString_,currentReadCnt );
  		alignerObj.kMaps_.runCutOff_ = currentRunCutoff;
  		/**@todo need to sort here*/
  		if(onPerId){
  			runClusteringOnId(clusters, kClus.readPositions_, iteratorMap, alignerObj);
  		}else{
  			runClustering(clusters, kClus.readPositions_, iteratorMap, alignerObj);
  		}
  		alignerObj.kMaps_.runCutOff_ = oldRunCutoff;
  	}
  	opts_.verbose_ = oldVerbose;
  	//remove the clusters marked removed
  	clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
  }

  readVecSorter::sortReadVector(clusters, "totalCount");
  //create snapshots directory if taking snap shots
  std::string snapShotsDirectoryName = "";
  if (snapShots) {
    snapShotsDirectoryName = bib::files::makeDir(setUp.directoryName_, snapShotsDirName);
  }
  for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {
  	if(opts_.verbose_){
  		std::cout << std::endl;
  	}
  	uint32_t sizeOfReadVector = readVec::getReadVectorSize(clusters);
    runningParameters runPars(iteratorMap[i],i,  sizeOfReadVector, onPerId);
    if(opts_.verbose_){
    	runPars.printIterInfo(std::cout, true,onPerId);
    }
    // if collapse is only on indels, create cluerters based on condensed seq and run clustering on those
    if (!onPerId && runPars.errors_.hqMismatches_ == 0 &&
        runPars.errors_.lqMismatches_ == 0 &&
        runPars.errors_.largeBaseIndel_ == 0) {
      auto byCondensed = readVec::organizeByCondensedSeq(clusters);
      bool oldVerbose = opts_.verbose_;
    	uint32_t startingSizeOfReadVector = 0;
    	for(const auto & read : clusters){
    		if(!read.remove){
    			++startingSizeOfReadVector;
    		}
    	}
      if (opts_.verbose_){
        std::cout << "Starting with " << startingSizeOfReadVector << " clusters"
                  << std::endl;
      }
      opts_.verbose_ = setUp.debug_;
  		for (auto & condensedReads : byCondensed){
  			if (onPerId) {
  				collapseWithPerId(condensedReads.second, runPars, alignerObj, false);
  			} else {
  				collapseWithParameters(condensedReads.second, runPars, alignerObj, false);
  			}
  		}
  		opts_.verbose_ = oldVerbose;
      clusters.clear();
      for (const auto& condensedReads : byCondensed) {
        addOtherVec(clusters, condensedReads.second);
      }
  		bib::stopWatch watch;
  		watch.setLapName("sorting vector");
  		readVecSorter::sortReadVector(clusters, "totalCount");
  		if(opts_.debug_){
  			watch.logLapTimes(std::cout, true, 6,true);
  		}
      int stopSizeOfReadVector = readVec::getReadVectorSize(clusters);
      if(opts_.verbose_){
      	std::cout << "Collapsed down to " << stopSizeOfReadVector << std::endl;
      }
    } else {
			if (onPerId) {
				collapseWithPerId(clusters, runPars, alignerObj, true);
			} else {
				collapseWithParameters(clusters, runPars, alignerObj, true);
			}
    }
    //write out current iteration if taking snapshots
    if (snapShots) {
      std::string iterDir =
          bib::files::makeDir(snapShotsDirectoryName, std::to_string(i));
      std::vector<cluster> currentClusters =
          readVecSplitter::splitVectorOnRemove(clusters).first;
      std::string seqName = bib::files::getFileName(setUp.ioOptions_.firstName_);
      renameReadNames(currentClusters, seqName, true, false);
      readObjectIO::write(currentClusters, snapShotsDirectoryName + std::to_string(i),
      		setUp.ioOptions_);
      clusterVec::allWriteOutClusters(currentClusters, iterDir, setUp.ioOptions_);
      if (setUp.refFilename_ == "") {
        profiler::getFractionInfoCluster(currentClusters, snapShotsDirectoryName,
                                  std::to_string(i) + ".tab.txt");
      } else {
        profiler::getFractionInfoCluster(currentClusters, snapShotsDirectoryName,
                                  std::to_string(i) + ".tab.txt",
																	setUp.refFilename_, alignerObj, false,
                                  true);
      }
    }
    //log info
    if(opts_.verbose_){
    	std::cout << "Current duration: ";
    	setUp.logRunTime(std::cout);
    }
    setUp.rLog_ << "Current iteration: " << i << "\n";
    setUp.rLog_ << "\tCurrent duration: "; setUp.rLog_.logTotalTime(6);
    setUp.rLog_ << "\tCurrent clusters size " << sizeOfReadVector << "\n";
  }
  clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
}

void collapser::markChimerasAdvanced(std::vector<cluster> &processedReads,
		aligner &alignerObj, double parentFreqs, int runCutOff,
		const comparison &chiOverlap, uint32_t overLapSizeCutoff,
		uint32_t &chimeraCount, uint32_t allowableError) {
	//assumes clusters are coming in sorted by total count
	if(processedReads.size() <=2){
		return;
	}
	if(opts_.verbose_){
		std::cout << "Marking Chimeras" << std::endl;
		std::cout << "Initial Pass" << std::endl;
	}
	for(const auto & subPos : iter::range<size_t>(2, processedReads.size())){
  	if(opts_.verbose_){
  		std::cout << subPos << ":" << processedReads.size()<< "\r";
  		std::cout.flush();
  	}
  	for(const auto &pos : iter::range<size_t>(0,subPos)){
			if(opts_.verbose_ && opts_.debug_){
				std::cout << std::endl;
				std::cout << processedReads[pos].seqBase_.name_ << std::endl;
				std::cout << processedReads[subPos].seqBase_.name_ << std::endl;
				std::cout << "processedReads[pos].seqBase_.cnt_" << processedReads[pos].seqBase_.cnt_ << std::endl;
				std::cout << "processedReads[pos].seqBase_.frac_" << processedReads[pos].seqBase_.frac_ << std::endl;
				std::cout << "processedReads[subPos].seqBase_.cnt_" << processedReads[subPos].seqBase_.cnt_ << std::endl;
				std::cout << "processedReads[subPos].seqBase_.frac_" << processedReads[subPos].seqBase_.frac_ << std::endl;
			}
      if ((processedReads[pos].seqBase_.cnt_ /
           processedReads[subPos].seqBase_.cnt_) < parentFreqs ||
          processedReads[subPos].seqBase_.cnt_ <= runCutOff) {
        continue;
      }
      alignerObj.alignVec(processedReads[pos], processedReads[subPos], false);
			alignerObj.profileAlignment(processedReads[pos], processedReads[subPos],
					11, true, false, true, false, opts_.weighHomopolyer_);
			if(opts_.verbose_ && opts_.debug_){
				std::cout << "alignerObj.mismatches_.size()" << alignerObj.mismatches_.size() << std::endl;
				std::cout << "alignerObj.alignmentGaps_.size()" << alignerObj.alignmentGaps_.size() << std::endl;
			}
      if (alignerObj.mismatches_.size() > allowableError
      		|| alignerObj.comp_.twoBaseIndel_ > allowableError
					|| alignerObj.comp_.largeBaseIndel_ > allowableError) {
      	std::map<uint32_t, mismatch> savedMismatches = alignerObj.mismatches_;
      	auto savedGapInfo = alignerObj.alignmentGaps_;
      	if(savedMismatches.size() > allowableError){
        	bool passEnd = true;
        	bool passFront = true;

        	//check front until first mismatch
          if (savedMismatches.begin()->first + 1 <= overLapSizeCutoff){
          	passFront = false;
          }else{
            alignerObj.profileAlignment(
                processedReads[pos], processedReads[subPos], 11, true, false, true,
                false, opts_.weighHomopolyer_, 0, savedMismatches.begin()->first);
            passFront = chiOverlap.passErrorProfile(alignerObj.comp_);
            if(passFront){
            	processedReads[subPos].frontChiPos.insert(
            	                        {savedMismatches.begin()->second.seqBasePos, pos});
            }
          }
          //check end from last mismatch
          if(overLapSizeCutoff >=
              (alignerObj.alignObjectA_.seqBase_.seq_.size() -
              		savedMismatches.rbegin()->first)){
          	passEnd = false;
          }else{
            alignerObj.profileAlignment(
                processedReads[pos], processedReads[subPos], 11, true, false, true,
                false, opts_.weighHomopolyer_, savedMismatches.rbegin()->first + 1);
            passEnd = chiOverlap.passErrorProfile(alignerObj.comp_);
            if(passEnd){
            	processedReads[subPos].endChiPos.insert(
            	                        {savedMismatches.rbegin()->second.seqBasePos, pos});
            }
          }
      	}

      	if(savedGapInfo.size() > allowableError){
          //check front until first large gap
          for(const auto & g : savedGapInfo){
          	if(g.second.size_ <=1){
          		continue;
          	}
  					alignerObj.profileAlignment(processedReads[pos],
  							processedReads[subPos], 11, true, false, true, false,
  							opts_.weighHomopolyer_, 0, g.second.startPos_);
  					if (chiOverlap.passErrorProfile(alignerObj.comp_)) {
  						processedReads[subPos].frontChiPos.insert(
  								{ getRealPos(alignerObj.alignObjectB_.seqBase_.seq_, g.second.startPos_), pos });
  					}
          	break;
          }

  				//check end from last large gap
  				for (const auto & g : iter::reverse(savedGapInfo)) {
  					if (g.second.size_ <= 1) {
  						continue;
  					}

  					alignerObj.profileAlignment(processedReads[pos],
  							processedReads[subPos], 11, true, false, true, false,
  							opts_.weighHomopolyer_, g.second.startPos_ + g.second.size_);

  					if (chiOverlap.passErrorProfile(alignerObj.comp_)) {
  						if (g.second.ref_) {
  							processedReads[subPos].endChiPos.insert(
  									{ getRealPos(alignerObj.alignObjectB_.seqBase_.seq_,
  											g.second.startPos_) + g.second.size_ - 1, pos });
  						} else {
  							/**@todo something else should probably happen if this is 0 */
  							if (0
  									== getRealPos(alignerObj.alignObjectB_.seqBase_.seq_,
  											g.second.startPos_)) {
  								processedReads[subPos].endChiPos.insert(
  										{ getRealPos(alignerObj.alignObjectB_.seqBase_.seq_,
  												g.second.startPos_), pos });
  							} else {
  								processedReads[subPos].endChiPos.insert(
  										{ getRealPos(alignerObj.alignObjectB_.seqBase_.seq_,
  												g.second.startPos_) - 1, pos });
  							}
  						}
  					}
  					break;
  				}
      	}
      }
    }
  }
  if(opts_.verbose_){
  	std::cout << std::endl;
  }

  if(opts_.debug_ && opts_.debug_){
    for (const auto &clusPos : iter::range(processedReads.size())) {
    	auto & clus = processedReads[clusPos];
    	if(clus.endChiPos.empty() || clus.frontChiPos.empty()){
    		continue;
    	}
    	auto endPos = getVectorOfMapKeys(clus.endChiPos);
    	auto frontPos = getVectorOfMapKeys(clus.frontChiPos);
    	auto endMin = std::min_element(endPos.begin(), endPos.end());
    	auto frontMax = std::max_element(frontPos.begin(), frontPos.end());
    	std::cout << clusPos << ":" << processedReads.size()<< std::endl;;
    	std::cout << clus.seqBase_.name_ << std::endl;
    	std::cout << "\tends" << std::endl;
    	std::cout << "\t";
    	for(const auto & endPosition : clus.endChiPos){
    		std::cout << endPosition.first << ":" << endPosition.second << ", ";
    	}
    	std::cout << std::endl;
    	std::cout << "\t" << *endMin << std::endl;
    	std::cout << "\tfronts" << std::endl;
    	std::cout << "\t";
    	for(const auto & frontPosition : clus.frontChiPos){
    		std::cout << frontPosition.first << ":" << frontPosition.second << ", ";
    	}
    	std::cout << std::endl;
    	std::cout << "\t" << *frontMax << std::endl;

    }
  }

	if(opts_.verbose_){
		std::cout << "Second Pass" << std::endl;
	}

  for (const auto &clusPos : iter::range(processedReads.size())) {
  	auto & clus = processedReads[clusPos];
  	if(opts_.verbose_){
  		std::cout << clusPos << ":" << processedReads.size()<< "\r";
  		std::cout.flush();
  	}
  	if(clus.endChiPos.empty() || clus.frontChiPos.empty()){
  		continue;
  	}
  	auto endPos = getVectorOfMapKeys(clus.endChiPos);
  	auto frontPos = getVectorOfMapKeys(clus.frontChiPos);
  	auto endMin = std::min_element(endPos.begin(), endPos.end());
  	auto frontMax = std::max_element(frontPos.begin(), frontPos.end());

  	if((*endMin) < (*frontMax) ){
      clus.seqBase_.markAsChimeric();
      ++chimeraCount;
  	}
  }
}


}  // namespace bibseq



