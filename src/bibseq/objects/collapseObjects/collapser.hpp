#pragma once
//
//  collapser.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/1/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/objects/collapseObjects/collapserOpts.hpp"
#include "bibseq/alignment.h"
#include "bibseq/seqToolsUtils.h"
#include "bibseq/readVectorManipulation.h"
#include "bibseq/objects/kmer/kmerCalculator.hpp"
#include "bibseq/objects/seqObjects/Clusters/cluster.hpp"
#include "bibseq/programUtils/seqSetUp.hpp"

namespace bibseq {

class collapser {

public:
	// constructor
	collapser(bool findingBestMatch, uint32_t bestMatchCheck, bool local,
			bool checkKmers, bool kmersByPosition, uint32_t runCutOff,
			uint32_t kLength, bool verbose, bool smallestFirst,
			bool condensedCollapse, bool weighHomopolyer,
			bool skipOnLetterCounterDifference, double fractionDifferenceCutOff,
			bool adjustHomopolyerRuns) :
			opts_(findingBestMatch, bestMatchCheck, local, checkKmers,
					kmersByPosition, runCutOff, kLength, verbose, smallestFirst,
					condensedCollapse, weighHomopolyer, skipOnLetterCounterDifference,
					fractionDifferenceCutOff, adjustHomopolyerRuns) {
	}
	collapser(const collapserOpts & opts) :
			opts_(opts) {
	}

  collapserOpts opts_;
private:
  template<class CLUSTER>
	void findMatch(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
			const std::vector<uint64_t> & positions,
			const runningParameters &runParams, size_t &amountAdded,
			aligner &alignerObj);

	template<class CLUSTER>
	void findMatchOnPerId(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
			const std::vector<uint64_t> & positions,
			const runningParameters &runParams, size_t &amountAdded,
			aligner &alignerObj);


  template <class CLUSTER>
  void collapseWithParameters(std::vector<CLUSTER> &comparingReads,
                                     const runningParameters &runParams,
                                     aligner &alignerObj, bool sort);
  template <class CLUSTER>
  void collapseWithParameters(
      std::vector<CLUSTER> &comparingReads, std::vector<uint64_t> & positions,
      const runningParameters &runParams,
      aligner &alignerObj);

  template <class CLUSTER>
	void collapseWithPerId(
	    std::vector<CLUSTER> &comparingReads, const runningParameters &runParams,
	    aligner &alignerObj, bool sort);

  template <class CLUSTER>
	void collapseWithPerId(
	    std::vector<CLUSTER> &comparingReads,
	    std::vector<uint64_t> & positions,
	    const runningParameters &runParams,
	    aligner &alignerObj);


public:

  template <class CLUSTER>
  std::vector<CLUSTER> runClustering(std::vector<CLUSTER> &currentClusters,
  		std::map<int, std::vector<double>> iteratorMap, aligner &alignerObj);

  template <class CLUSTER>
  void runClustering(std::vector<CLUSTER> &currentClusters,
  		std::vector<uint64_t> & positons,
  		std::map<int, std::vector<double>> iteratorMap, aligner &alignerObj);

  template <class CLUSTER>
  std::vector<CLUSTER> collapseCluster(
      const std::vector<CLUSTER> & clusters,
      std::map<int, std::vector<double>> iteratorMap, const std::string &sortBy,
      aligner &alignerObj);

  template <class CLUSTER>
  std::vector<CLUSTER> collapseCluster(
      std::vector<CLUSTER> clusters,
      std::map<int, runningParameters> iteratorMap, const std::string &sortBy,
      aligner &alignerObj);


  template <class CLUSTER>
	void runClusteringOnId(std::vector<CLUSTER> &currentClusters,
			std::vector<uint64_t> & positions,
			std::map<int, std::vector<double>> iteratorMap,
			aligner &alignerObj);
  template <class CLUSTER>
  std::vector<CLUSTER> collapseClusterOnPerId(
      std::vector<CLUSTER> clusters,
      std::map<int, runningParameters> iteratorMap, const std::string &sortBy,
      aligner &alignerObj);
  template <class CLUSTER>
  std::vector<CLUSTER> collapseClusterOnPerId(
      const std::vector<CLUSTER> & clusters,
      std::map<int, std::vector<double>> iteratorMap, const std::string &sortBy,
      aligner &alignerObj);

  void runFullClustering(std::vector<cluster> & clusters,
  		bool onPerId, std::map<int, std::vector<double>> iteratorMap,
			std::map<int, std::vector<double>> binIteratorMap,
  		bool useNucComp,bool useMinLenNucComp, bool findBestNuc, const std::vector<double> & diffCutOffVec,
  		bool useKmerBinning, uint32_t kCompareLen, double kmerCutOff,
  		aligner & alignerObj, const SeqSetUpPars & setUpPars,
			bool snapShots, const std::string & snapShotsDirName);


  void markChimerasAdvanced(std::vector<cluster> &processedReads,
                                   aligner &alignerObj,
																	 double parentFreqs,
                                   int runCutOff,
                                   const comparison &chiOverlap,
                                   uint32_t overLapSizeCutoff,
                                   uint32_t &chimeraCount,
                                   uint32_t allowableError);

private:
	template<class CLUSTER>
	void findMatchOneOff(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
			const runningParameters &runParams, double fracCutoff,
			size_t &amountAdded, aligner &alignerObj);
	template<class CLUSTER>
	std::vector<CLUSTER> collapseOneOffs(std::vector<CLUSTER> clusters,
			runningParameters runParams, double fracCutOff, aligner &alignerObj);

};

template <class CLUSTER>
void collapser::collapseWithPerId(
    std::vector<CLUSTER> &comparingReads,
    std::vector<uint64_t> & positions,
    const runningParameters &runParams,
    aligner &alignerObj) {
	int32_t sizeOfReadVector = 0;
	for(const auto & pos : positions){
		if(!comparingReads[pos].remove){
			++sizeOfReadVector;
		}
	}
  if (sizeOfReadVector < 2) {
    return;
  }
  if (opts_.verbose_){
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  }
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (opts_.smallestFirst_) {
    for (const auto &reverseReadPos : iter::reverse(positions)) {
    	auto & reverseRead = comparingReads[reverseReadPos];
      if (reverseRead.remove) {
        continue;
      } else {
        ++clusterCounter;
      }
      if (opts_.verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << "\r";
        std::cout.flush();
      }
      findMatchOnPerId(reverseRead, comparingReads, positions,
      		runParams, amountAdded, alignerObj);
    }
  } else {
    for (const auto &forwardReadPos : positions) {
    	auto & forwardRead = comparingReads[forwardReadPos];
      if (forwardRead.remove) {
        continue;
      } else {
        ++clusterCounter;
      }
      if (opts_.verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << "\r";
        std::cout.flush();
      }
      findMatchOnPerId(forwardRead, comparingReads,positions,
      		runParams, amountAdded,alignerObj);
    }
  }
  bib::stopWatch watch;
  watch.setLapName("updatingName");
	for(const auto & pos : positions){
		comparingReads[pos].updateName();
	}

	watch.startNewLap("allCalculateConsensus");
	for(const auto & pos : positions){
		comparingReads[pos].calculateConsensus(alignerObj, true);
	}
	watch.startNewLap("removeLowQualityBases");
	if (opts_.removeLowQualityBases_) {
		for(const auto & pos : positions){
			comparingReads[pos].seqBase_.removeLowQualityBases(opts_.lowQualityBaseTrim_);
		}
	}
	watch.startNewLap("adjustHomopolyerRuns");
	if (opts_.adjustHomopolyerRuns_) {
		for(const auto & pos : positions){
			comparingReads[pos].adjustHomopolyerRunQualities();
		}
	}
	/*
	watch.startNewLap("sortReadVector");
	auto comp = [&comparingReads](const uint64_t & pos1, const uint64_t & pos2){
		return comparingReads[pos1] < comparingReads[pos2];
	};
	std::sort(positions.begin(), positions.end(), comp);
	*/
  if(opts_.debug_){
  	watch.logLapTimes(std::cout, true, 6, true);
  }
  if (opts_.verbose_){
  	if(clusterCounter >= 100){
  		std::cout << std::endl;
  	}
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
  }
}

template<class CLUSTER>
void collapser::collapseWithPerId(std::vector<CLUSTER> &comparingReads,
		const runningParameters &runParams, aligner &alignerObj, bool sort) {
	std::vector<uint64_t> positions(comparingReads.size());
	bib::iota<uint64_t>(positions, 0);
	collapseWithPerId(comparingReads, positions, runParams, alignerObj);
	if(sort){
		bib::stopWatch watch;
		watch.setLapName("sorting vector");
		readVecSorter::sortReadVector(comparingReads, "totalCount");
		if(opts_.debug_){
			watch.logLapTimes(std::cout, true, 6,true);
		}
	}
}


template <class CLUSTER>
void collapser::findMatch(CLUSTER &read,
                                 std::vector<CLUSTER> &comparingReads,
                                 const std::vector<uint64_t> & positions,
                                 const runningParameters &runParams,
                                 size_t &amountAdded,
																 aligner &alignerObj) {
	//std::cout << "FindMatch start" << std::endl;
  int count = -1;
  double bestScore = 0;
  bool foundMatch = false;
  uint32_t bestSearching = 0;
  uint32_t bestClusterPos = std::numeric_limits<uint32_t>::max();
  for (const auto &clusPos : positions) {
    if (comparingReads[clusPos].remove) {
      continue;
    }
  	auto & clus = comparingReads[clusPos];
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      continue;
      //break;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (opts_.skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > opts_.fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 20) {
        continue;
      }
    }

    if (opts_.useReadLen_){
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > opts_.readLenDiff_) {
        continue;
      }
    }

    if (foundMatch) {
      ++bestSearching;
    }
    ++count;
    if (1 + count > runParams.stopCheck_ || bestSearching > opts_.bestMatchCheck_) {
      break;
    }
    bool matching = clus.compare(read, alignerObj, runParams.errors_, opts_);
    if (matching) {
      foundMatch = true;
      if (opts_.findingBestMatch_) {
        if(opts_.noAlign_){
        	alignerObj.noAlignSetAndScore(clus, read);
        }else{
        	alignerObj.alignCacheGlobal(clus, read);
        }
      	double currentScore = 0;
      	if(opts_.eventBased_){
      		alignerObj.profilePrimerAlignment(clus, read, opts_.weighHomopolyer_);
      		currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
      	}else{
      		currentScore = alignerObj.parts_.score_;
      	}
        if (currentScore > bestScore) {
          bestScore = alignerObj.parts_.score_;
          bestClusterPos = clusPos;
        }
      } else {
        ++amountAdded;
        clus.addRead(read);
        read.remove = true;
        break;
      }
    }
  }
  if (foundMatch) {
    if (opts_.findingBestMatch_) {
    	if(bestClusterPos != std::numeric_limits<uint32_t>::max()){
    		comparingReads[bestClusterPos].addRead(read);
        read.remove = true;
        ++amountAdded;
    	}
    }
  }
}

template <class CLUSTER>
void collapser::findMatchOnPerId(CLUSTER &read,
                                 std::vector<CLUSTER> &comparingReads,
                                 const std::vector<uint64_t> & positions,
                                 const runningParameters &runParams,
                                 size_t &amountAdded, aligner &alignerObj) {
  int count = -1;
  double bestScore = 0;
  bool foundMatch = false;
  uint32_t bestSearching = 0;
  uint32_t bestClusterPos = std::numeric_limits<uint32_t>::max();
  for (const auto &clusPos : positions) {
    if (comparingReads[clusPos].remove) {
      continue;
    }
  	auto & clus = comparingReads[clusPos];
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      continue;
      //break;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (opts_.skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > opts_.fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 10) {
        continue;
      }
    }

    if (opts_.useReadLen_){
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > opts_.readLenDiff_) {
        continue;
      }
    }

    if (foundMatch) {
      ++bestSearching;
    }



    ++count;
    if (1 + count > runParams.stopCheck_ || bestSearching > opts_.bestMatchCheck_) {
      break;
    }


    bool matching = clus.compareId(read, alignerObj,runParams.errors_,opts_);
    if (matching) {
      foundMatch = true;
      if (opts_.findingBestMatch_) {
        if(opts_.noAlign_){
        	alignerObj.noAlignSetAndScore(clus, read);
        }else{
        	alignerObj.alignCacheGlobal(clus, read);
        }
      	double currentScore = 0;
      	if(opts_.eventBased_){
      		alignerObj.profilePrimerAlignment(clus, read, opts_.weighHomopolyer_);
      		currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
      	}else{
      		currentScore = alignerObj.parts_.score_;
      	}
        if (currentScore > bestScore) {
          bestScore = alignerObj.parts_.score_;
          bestClusterPos = clusPos;
        }
      } else {
        ++amountAdded;
        clus.addRead(read);
        //clus.allInputClusters.push_back(read);
        read.remove = true;
        break;
      }
    }
  }
  if (foundMatch) {
    if (opts_.findingBestMatch_) {
    	if(bestClusterPos != std::numeric_limits<uint32_t>::max()){
    		comparingReads[bestClusterPos].addRead(read);
    		//comparingReads[bestClusterPos].allInputClusters.emplace_back(read);
        read.remove = true;
        ++amountAdded;
    	}
    }
  }
}



template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseClusterOnPerId(
    std::vector<CLUSTER> clusters,
    std::map<int, runningParameters> iteratorMap, const std::string &sortBy,
    aligner &alignerObj){
  // set kmer maps
  //kmerMaps kMaps = kmerCalculator::indexKmerMpas(clusters, opts_.kLength_, opts_.runCutOff_);
  //alignerObj.setKmerMpas(kMaps);
  // go over the iterations collapsing down and down, sort the reads after each
  // iteration
  // and calculate a consensus after each one
  readVec::allSetLetterCount(clusters);
  for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {
    if (opts_.verbose_) {
    	iteratorMap[i].printIterInfo(std::cout, true, true);
    }
    collapseWithPerId(clusters,iteratorMap[i], alignerObj, true);
  }
  return readVecSplitter::splitVectorOnRemove(clusters).first;
}


template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseClusterOnPerId(
		const std::vector<CLUSTER> & clusters,
    std::map<int, std::vector<double>> iteratorMap, const std::string &sortBy,
    aligner &alignerObj) {
	std::map<int, runningParameters> iteratorMapPars;
	bool onPerId = false;
	for(const auto & iteration : iteratorMap){
		iteratorMapPars[iteration.first] = runningParameters(iteration.second,iteration.first,
				readVec::getTotalReadCount(clusters), onPerId);
	}
  return collapseClusterOnPerId(clusters, iteratorMapPars, sortBy, alignerObj);
}

template <class CLUSTER>
void collapser::runClusteringOnId(std::vector<CLUSTER> &currentClusters,
		std::vector<uint64_t> & positions,
		std::map<int, std::vector<double>> iteratorMap,
		aligner &alignerObj){
	{
		bib::stopWatch watch;
		auto comp = [&currentClusters](const uint64_t & pos1, const uint64_t & pos2){
			return currentClusters[pos1] < currentClusters[pos2];
		};
		watch.setLapName("sortReadVector");
		bib::sort(positions, comp);
		if(opts_.debug_){
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
	for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {
		uint32_t sizeOfReadVector = 0;
		for(const auto & pos : positions){
			if(!currentClusters[pos].remove){
				++sizeOfReadVector;
			}
		}
		bool onPerId = false;
		runningParameters runPars(iteratorMap[i], i, sizeOfReadVector, onPerId);
		if(opts_.verbose_){
			std::cout << std::endl;
			runPars.printIterInfo(std::cout, true, true);
		}
		collapseWithPerId(currentClusters, positions, runPars, alignerObj);
		{
			bib::stopWatch watch;
			auto comp = [&currentClusters](const uint64_t & pos1, const uint64_t & pos2){
				return currentClusters[pos1] < currentClusters[pos2];
			};
			watch.setLapName("sortReadVector");
			bib::sort(positions, comp);
			if(opts_.debug_){
				watch.logLapTimes(std::cout, true, 6, true);
			}
		}
	}
}


template <class CLUSTER>
std::vector<CLUSTER> collapser::runClustering(std::vector<CLUSTER> &currentClusters,
		std::map<int, std::vector<double>> iteratorMap, aligner &alignerObj){
	for (int i = 1; i <= iteratorMap.size(); ++i) {
		uint32_t sizeOfReadVector = readVec::getReadVectorSize(currentClusters);
		bool onPerId = false;
		runningParameters runPars(iteratorMap[i],i, sizeOfReadVector, onPerId);
		if(opts_.verbose_){
			std::cout << std::endl;
			runPars.printIterInfo(std::cout, true, false);
		}
		if (opts_.condensedCollapse_ && runPars.errors_.hqMismatches_ == 0 &&
				runPars.errors_.lqMismatches_ == 0 &&
				runPars.errors_.largeBaseIndel_ == 0) {
			auto byCondensed = readVec::organizeByCondensedSeq(currentClusters);
			for (auto& condensedReads : byCondensed) {
				collapseWithParameters(condensedReads.second,runPars, alignerObj, false);
			}
			currentClusters.clear();
			for (const auto& condensedReads : byCondensed) {
				addOtherVec(currentClusters, condensedReads.second);
			}
			bib::stopWatch watch;
			watch.setLapName("sorting vector");
			readVecSorter::sortReadVector(currentClusters, "totalCount");
			if(opts_.debug_){
				watch.logLapTimes(std::cout, true, 6,true);
			}
		} else {
			collapseWithParameters(currentClusters, runPars, alignerObj, true);
		}
	}
	return readVecSplitter::splitVectorOnRemove(currentClusters).first;
}

template <class CLUSTER>
void collapser::runClustering(std::vector<CLUSTER> &currentClusters,
		std::vector<uint64_t> & positions,
		std::map<int, std::vector<double>> iteratorMap,
		aligner &alignerObj){
	{
		bib::stopWatch watch;
		auto comp = [&currentClusters](const uint64_t & pos1, const uint64_t & pos2){
			return currentClusters[pos1] < currentClusters[pos2];
		};
		watch.setLapName("sortReadVector");
		bib::sort(positions, comp);
		if(opts_.debug_){
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
	for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {
		uint32_t sizeOfReadVector = 0;
		for(const auto & pos : positions){
			if(!currentClusters[pos].remove){
				++sizeOfReadVector;
			}
		}
		bool onPerId = false;
		runningParameters runPars(iteratorMap[i],i, sizeOfReadVector, onPerId);
		if(opts_.verbose_){
			std::cout << std::endl;
			runPars.printIterInfo(std::cout, true, false);
		}
		if (opts_.condensedCollapse_ && runPars.errors_.hqMismatches_ == 0 &&
				runPars.errors_.lqMismatches_ == 0 &&
				runPars.errors_.largeBaseIndel_ == 0) {
			auto byCondensed = readVec::organizeByCondensedSeqPositions(currentClusters, positions);
			for (auto& condensedReads : byCondensed) {

				collapseWithParameters(currentClusters,condensedReads.second,runPars, alignerObj);
			}
		} else {
			collapseWithParameters(currentClusters, positions, runPars,alignerObj);
		}
		bib::stopWatch watch;
		auto comp = [&currentClusters](const uint64_t & pos1, const uint64_t & pos2){
			return currentClusters[pos1] < currentClusters[pos2];
		};
		watch.setLapName("sortReadVector");
		bib::sort(positions, comp);
		if(opts_.debug_){
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
}

template <class CLUSTER>
void collapser::collapseWithParameters(
    std::vector<CLUSTER> &comparingReads, const runningParameters &runParams,
    aligner &alignerObj, bool sort) {
	std::vector<uint64_t> positions(comparingReads.size());
	bib::iota<uint64_t>(positions, 0);
	collapseWithParameters(comparingReads, positions, runParams, alignerObj);
	if(sort){
		bib::stopWatch watch;
		watch.setLapName("sorting vector");
		readVecSorter::sortReadVector(comparingReads, "totalCount");
		if(opts_.debug_){
			watch.logLapTimes(std::cout, true, 6,true);
		}
	}
}

template <class CLUSTER>
void collapser::collapseWithParameters(
    std::vector<CLUSTER> &comparingReads,
    std::vector<uint64_t> & positions,
    const runningParameters &runParams,
    aligner &alignerObj) {
	uint32_t sizeOfReadVector = 0;
	for(const auto & pos : positions){
		if(!comparingReads[pos].remove){
			++sizeOfReadVector;
		}
	}
  if (sizeOfReadVector < 2) {
    return;
  }
  if (opts_.verbose_){
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  }
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (opts_.smallestFirst_) {
    for (const auto &reverseReadPos : iter::reverse(positions)) {
    	auto & reverseRead = comparingReads[reverseReadPos];
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (opts_.verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << "\r";
        std::cout.flush();
      }
      findMatch(reverseRead, comparingReads, positions,
      		runParams, amountAdded, alignerObj);
    }
  } else {
    for (const auto &forwardReadPos : positions) {
    	auto & forwardRead = comparingReads[forwardReadPos];
      if (forwardRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (opts_.verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << "\r";
        std::cout.flush();
      }
      findMatch(forwardRead, comparingReads,positions,
      		runParams, amountAdded,
                       alignerObj);
    }
  }

  bib::stopWatch watch;
  watch.setLapName("updatingName");
	for(const auto & pos : positions){
		comparingReads[pos].updateName();
	}
	watch.startNewLap("allCalculateConsensus");
	for(const auto & pos : positions){
		comparingReads[pos].calculateConsensus(alignerObj, true);
	}
	watch.startNewLap("removeLowQualityBases");
	if (opts_.removeLowQualityBases_) {
		for(const auto & pos : positions){
			comparingReads[pos].seqBase_.removeLowQualityBases(opts_.lowQualityBaseTrim_);
		}
	}
	watch.startNewLap("adjustHomopolyerRuns");
	if (opts_.adjustHomopolyerRuns_) {
		for(const auto & pos : positions){
			comparingReads[pos].adjustHomopolyerRunQualities();
		}
	}
	/*
	watch.startNewLap("sortReadVector");
	auto comp = [&comparingReads](const uint64_t & pos1, const uint64_t & pos2){
		return comparingReads[pos1] < comparingReads[pos2];
	};
	bib::sort(positions, comp);
	 */
  if(opts_.debug_){
  	watch.logLapTimes(std::cout, true, 6, true);
  }

  if (opts_.verbose_){
  	if(clusterCounter >= 100){
  		std::cout << std::endl;
  	}
  	 std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
  	              << " clusters" << std::endl;
  }
}


template<class CLUSTER>
std::vector<CLUSTER> collapser::collapseCluster(std::vector<CLUSTER> clusters,
		std::map<int, runningParameters> iteratorMap, const std::string &sortBy,
		aligner &alignerObj) {
	//set kmer frequncy maps
	bool expandKmerPos = false;
	uint32_t expandKmerPosSize = 5;
	KmerMaps kMaps = indexKmers(clusters, opts_.kLength_, opts_.runCutOff_, opts_.kmersByPosition_, expandKmerPos, expandKmerPosSize);
	alignerObj.setKmerMpas(kMaps);
	// go over the iterations collapsing down and down, sort the reads after each
	// iteration
	// and calculate a consensus after each one
	readVec::allSetLetterCount(clusters);
	for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {
		runningParameters runPars(iteratorMap[i]);
		if (opts_.verbose_) {
			runPars.printIterInfo(std::cout, true, false);
		}
		if (opts_.condensedCollapse_ && runPars.errors_.hqMismatches_ == 0
				&& runPars.errors_.lqMismatches_ == 0
				&& runPars.errors_.largeBaseIndel_ == 0) {
			auto byCondensed = readVec::organizeByCondensedSeq(clusters);
			for (auto &condensedReads : byCondensed) {
				collapseWithParameters(condensedReads.second, runPars, alignerObj, false);
			}
			clusters.clear();
			for (const auto &condensedReads : byCondensed) {
				addOtherVec(clusters, condensedReads.second);
			}
			bib::stopWatch watch;
			watch.setLapName("sorting vector");
			readVecSorter::sortReadVector(clusters, "totalCount");
			if(opts_.debug_){
				watch.logLapTimes(std::cout, true, 6,true);
			}
		} else {
			collapseWithParameters(clusters, runPars, alignerObj, true);
		}
	}
	return readVecSplitter::splitVectorOnRemove(clusters).first;
}

template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseCluster(
    const std::vector<CLUSTER> & clusters,
    std::map<int, std::vector<double>> iteratorMap, const std::string &sortBy,
    aligner &alignerObj) {
	std::map<int, runningParameters> iteratorMapPars;
	bool onPerId = false;
	for(const auto & iteration : iteratorMap){
		iteratorMapPars[iteration.first] = runningParameters(iteration.second,iteration.first,
				readVec::getTotalReadCount(clusters), onPerId);
	}
  return collapseCluster(clusters, iteratorMapPars, sortBy, alignerObj);
}




template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseOneOffs(
    std::vector<CLUSTER> clusters, runningParameters runParams,
    double fracCutOff, aligner &alignerObj) {
	throw std::runtime_error{std::string(__PRETTY_FUNCTION__) + " not yet implemented"};
  runParams.errors_.hqMismatches_ = 1;
	bool expandKmerPos = false;
	uint32_t expandKmerPosSize = 5;
	KmerMaps kMaps = indexKmers(clusters, opts_.kLength_, opts_.runCutOff_, opts_.kmersByPosition_, expandKmerPos, expandKmerPosSize);
  // std::cout <<"cc2" << std::endl;
  alignerObj.setKmerMpas(kMaps);
  // std::cout <<"cc3" << std::endl;
  // go over the iterations collapsing down and down, sort the reads after each
  // iteration
  // and calculate a consensus after each one

  // std::cout <<"cc4" << std::endl;
  int sizeOfReadVector = readVec::getReadVectorSize(clusters);
  if (opts_.verbose_)
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (opts_.smallestFirst_) {
    for (auto &reverseRead : iter::reverse(clusters)) {
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (opts_.verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(reverseRead, clusters, runParams, amountAdded,
                        alignerObj);
    }
  } else {
    for (auto &forwardRead : clusters) {
      if (forwardRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (opts_.verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(forwardRead, clusters, runParams, amountAdded,
                        alignerObj);
    }
  }
  if (opts_.verbose_){
  	if(clusterCounter >= 100){
  		std::cout << std::endl;
  	}
  	 std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
  	              << " clusters" << std::endl;
  }
  // collapser::collapseWithParameters(clusters, runPars, alignerObj);
  readVec::allUpdateName(clusters);
  readVecSorter::sortReadVector(clusters, "totalCount");
  return readVecSplitter::splitVectorOnRemove(clusters).first;
}






template <class CLUSTER>
void collapser::findMatchOneOff(CLUSTER &read,
                                std::vector<CLUSTER> &comparingReads,
                                const runningParameters &runParams,
                                double fracCutoff, size_t &amountAdded,
                                aligner &alignerObj) {
	throw std::runtime_error{std::string(__PRETTY_FUNCTION__) + " not yet implemented"};
	//std::cout << "FindMatch start" << std::endl;
  int count = -1;
  double bestScore = 0;
  int bestSearching = 0;
  bool foundMatch = false;
  uint64_t bestClusterPos = std::numeric_limits<uint64_t>::max();
  std::string bestCluster = "";
  // std::cout <<"fm 1 " << std::endl;
  for (const auto &clusPos : iter::range(comparingReads.size())) {
  	auto & clus = comparingReads[clusPos];
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      //continue;
      break;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (opts_.skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > opts_.fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 10) {
        continue;
      }
    }

    if (opts_.useReadLen_){
    	++count;
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > opts_.readLenDiff_) {
        continue;
      }
    }

    if (foundMatch) {
      ++bestSearching;
    }

    if (clus.remove) {
      continue;
    } else {
      ++count;
    }

    if (1 + count > runParams.stopCheck_ || bestSearching > opts_.bestMatchCheck_) {
      break;
    }
    //bool matching = clus.compare(read, alignerObj, runParams.errors_, opts_);
    bool matching = false;

    if (clus.previousErrorChecks_.find(read.firstReadName_) !=
            clus.previousErrorChecks_.end() &&
        read.previousErrorChecks_.find(clus.firstReadName_) !=
            read.previousErrorChecks_.end()) {
      matching = runParams.errors_.passErrorProfile(
          read.previousErrorChecks_.at(clus.firstReadName_));
      if(opts_.noAlign_){
      	alignerObj.noAlignSetAndScore(clus, read);
      }else{
      	alignerObj.alignCacheGlobal(clus, read);
      }
      //alignerObj.profilePrimerAlignment(clus, read, opts_.weighHomopolyer_);
    } else {
      if(opts_.noAlign_){
      	alignerObj.noAlignSetAndScore(clus, read);
      }else{
      	alignerObj.alignCacheGlobal(clus, read);
      }
      alignerObj.profilePrimerAlignment(clus, read, opts_.weighHomopolyer_);
      comparison currentProfile = alignerObj.compareAlignment(
          clus, read, runParams, opts_.checkKmers_, opts_.kmersByPosition_,
          opts_.weighHomopolyer_);
      read.previousErrorChecks_[clus.firstReadName_] = currentProfile;
      clus.previousErrorChecks_[read.firstReadName_] = currentProfile;
      matching = runParams.errors_.passErrorProfile(currentProfile);
    }

    if (matching) {
      foundMatch = true;
      if (opts_.findingBestMatch_) {
      	alignerObj.alignCacheGlobal(clus, read);
      	double currentScore = 0;
      	if(opts_.eventBased_){
      		alignerObj.profilePrimerAlignment(clus, read, opts_.weighHomopolyer_);
      		currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
      	}else{
      		currentScore = alignerObj.parts_.score_;
      	}
        if (currentScore > bestScore) {
          bestScore = alignerObj.parts_.score_;
          bestClusterPos = clusPos;
        }
      } else {
        ++amountAdded;
        clus.addRead(read);
        clus.allInputClusters.push_back(read);
        read.remove = true;
        break;
      }
    }

  }
  if (foundMatch) {
    if (opts_.findingBestMatch_) {
      read.remove = true;
      comparingReads[bestClusterPos].addRead(read);
     // readVec::getReadByName(comparingReads, bestCluster).addRead(read);
      ++amountAdded;
    }
  }
  //std::cout << "FindMatch stop" << std::endl;
}


}  // namespace bib
#ifndef NOT_HEADER_ONLY
#include "collapser.cpp"
#endif
