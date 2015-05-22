#pragma once
//
//  collapser.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/1/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/alignment.h"
#include "bibseq/seqToolsUtils.h"
#include "bibseq/readVectorManipulation.h"

namespace bibseq {

class collapser {

 public:
  // constructor
  collapser(bool findingBestMatch, int bestMatchCheck, bool local,
            bool checkKmers, bool kmersByPosition, int runCutOff,
            int qualRunCutOff, int kLength, bool verbose, bool smallestFirst,
            bool condensedCollapse, bool weighHomopolyer,
            bool skipOnLetterCounterDifference, double fractionDifferenceCutOff,
            bool regKmer, bool adjustHomopolyerRuns)
      : findingBestMatch_(findingBestMatch),
        bestMatchCheck_(bestMatchCheck),
        local_(local),
        checkKmers_(checkKmers),
        kmersByPosition_(kmersByPosition),
        runCutOff_(runCutOff),
        qualRunCutOff_(qualRunCutOff),
        kLength_(kLength),
        verbose_(verbose),
        smallestFirst_(smallestFirst),
        condensedCollapse_(condensedCollapse),
        weighHomopolyer_(weighHomopolyer),
        skipOnLetterCounterDifference_(skipOnLetterCounterDifference),
        fractionDifferenceCutOff_(fractionDifferenceCutOff),
        regKmer_(regKmer), adjustHomopolyerRuns_(adjustHomopolyerRuns)  {}
  // Members
  // aligner alignerObj;
  bool findingBestMatch_;
  int bestMatchCheck_;
  bool local_;
  bool checkKmers_;
  bool kmersByPosition_;
  int runCutOff_;
  int qualRunCutOff_;
  int kLength_;
  bool verbose_;
  bool smallestFirst_;
  bool condensedCollapse_;
  bool weighHomopolyer_;
  bool skipOnLetterCounterDifference_;
  double fractionDifferenceCutOff_;
  bool regKmer_;
  bool adjustHomopolyerRuns_;
  bool useReadLen_ = false;
  uint32_t readLenDiff_ = 10;
  bool eventBased = true;


  // functions


  template <class CLUSTER>
  std::vector<CLUSTER> runClustering(std::vector<CLUSTER> &currentClusters,
  		std::map<int, std::vector<double>> iteratorMap, aligner &alignerObj);

  template <class CLUSTER>
  void runClustering(std::vector<CLUSTER> &currentClusters,
  		std::vector<uint64_t> & positons,
  		std::map<int, std::vector<double>> iteratorMap, aligner &alignerObj);

  template <class CLUSTER>
  void findMatchRegKmer(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
                        const runningParameters &runParams, size_t &amountAdded,
                        aligner &alignerObj);
  template <class CLUSTER>
  void findMatchRegKmer(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
  											const std::vector<uint64_t> & positions,
                        const runningParameters &runParams, size_t &amountAdded,
                        aligner &alignerObj);



  template <class CLUSTER>
  void findMatchRegKmerThreaded(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
                        const runningParameters &runParams, size_t &amountAdded,
                        aligner &alignerObj);
  template <class CLUSTER>
  void findMatchQualKmer(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
  		const std::vector<uint64_t> & positions,
                         const runningParameters &runParams,
                         size_t &amountAdded, aligner &alignerObj);

  template <class CLUSTER>
  void findMatchQualKmer(CLUSTER &read, std::vector<CLUSTER> &comparingReads,

                         const runningParameters &runParams,
                         size_t &amountAdded, aligner &alignerObj);

  template <class CLUSTER>
  void findMatchOneOff(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
                       const runningParameters &runParams, double fracCutoff,
                       size_t &amountAdded, aligner &alignerObj);
  template <class CLUSTER>
  void collapseWithParametersRegKmer(std::vector<CLUSTER> &comparingReads,
                                     const runningParameters &runParams,
                                     aligner &alignerObj);
  template <class CLUSTER>
  void collapseWithParametersRegKmer(
      std::vector<CLUSTER> &comparingReads,const std::vector<uint64_t> & positions,
      const runningParameters &runParams,
      aligner &alignerObj);

  template <class CLUSTER>
  void collapseWithParametersRegKmerThreaded(std::vector<CLUSTER> &comparingReads,
                                     const runningParameters &runParams,
                                     std::vector<aligner> & aligners, uint32_t numThreads);
  template <class CLUSTER>
  void collapseWithParametersQualKmer(std::vector<CLUSTER> &comparingReads,
                                      const runningParameters &runParams,
                                      aligner &alignerObj);
  template <class CLUSTER>
  void collapseWithParametersQualKmer(std::vector<CLUSTER> &comparingReads,
  																		const std::vector<uint64_t> & positions,
                                      const runningParameters &runParams,
                                      aligner &alignerObj);

  // collapse samples
  template <class CLUSTER>
  void collapseSamples(const std::vector<CLUSTER> &clusters,
                       std::map<int, std::vector<double>> iteratorMap,
                       const std::string &sortBy, const std::string &sampName,
                       std::map<std::string, std::vector<CLUSTER>> &samples,
                       aligner &alignerObj);

  template <class CLUSTER>
  std::vector<CLUSTER> collapseCluster(
      std::vector<CLUSTER> clusters,
      std::map<int, std::vector<double>> iteratorMap, const std::string &sortBy,
      aligner &alignerObj);

  template <class CLUSTER>
  std::vector<CLUSTER> collapseCluster(
      std::vector<CLUSTER> clusters,
      std::map<int, runningParameters> iteratorMap, const std::string &sortBy,
      aligner &alignerObj);

  template <class CLUSTER>
  std::vector<CLUSTER> collapseOneOffsQualKmer(std::vector<CLUSTER> clusters,
                                               runningParameters runParams,
                                               double fracCutOff,
                                               aligner &alignerObj);
  template <class CLUSTER>
  std::vector<CLUSTER> collapseOneOffsRegKmer(std::vector<CLUSTER> clusters,
                                              runningParameters runParams,
                                              double fracCutOff,
                                              aligner &alignerObj);

  /*template <class CLUSTER>
  std::vector<CLUSTER> collapseOnOtu (std::vector<CLUSTER> clusters, double
  percentCutoff,
                                      bool countGaps, )*/
};
template <class CLUSTER>
std::vector<CLUSTER> collapser::runClustering(std::vector<CLUSTER> &currentClusters,
		std::map<int, std::vector<double>> iteratorMap, aligner &alignerObj){
	for (int i = 1; i <= (int)iteratorMap.size(); ++i) {
		uint32_t sizeOfReadVector = readVec::getReadVectorSize(currentClusters);
		runningParameters runPars(iteratorMap[i],i, sizeOfReadVector);
		if(verbose_){
			std::cout << std::endl;
			runPars.printIterInfo(std::cout, true);
		}
		if (condensedCollapse_ && runPars.errors_.hqMismatches_ == 0 &&
				runPars.errors_.lqMismatches_ == 0 &&
				runPars.errors_.largeBaseIndel_ == 0 &&
				runPars.errors_.twoBaseIndel_ == 0) {
			auto byCondensed = readVec::organizeByCondensedSeq(currentClusters);
			for (auto& condensedReads : byCondensed) {
				if (regKmer_) {
					collapseWithParametersRegKmer(condensedReads.second,
																										 runPars, alignerObj);
				} else {
					collapseWithParametersQualKmer(condensedReads.second,
																											runPars, alignerObj);
				}
			}
			currentClusters.clear();
			for (const auto& condensedReads : byCondensed) {
				addOtherVec(currentClusters, condensedReads.second);
			}
		} else {
			if (regKmer_) {
				collapseWithParametersRegKmer(currentClusters, runPars,
																									 alignerObj);
			} else {
				collapseWithParametersQualKmer(currentClusters, runPars,
																										alignerObj);
			}
		}
		readVec::allUpdateName(currentClusters);
		readVecSorter::sortReadVector(currentClusters, "totalCount");
		clusterVec::allCalculateConsensus(currentClusters, alignerObj, true);
		if (adjustHomopolyerRuns_) {
			readVec::allAdjustHomopolymerRunsQualities(currentClusters);
		}

	}
	return readVecSplitter::splitVectorOnRemove(currentClusters).first;
}

template <class CLUSTER>
void collapser::runClustering(std::vector<CLUSTER> &currentClusters,
		std::vector<uint64_t> & positions,
		std::map<int, std::vector<double>> iteratorMap,
		aligner &alignerObj){
	for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {
		uint32_t sizeOfReadVector = 0;
		for(const auto & pos : positions){
			if(!currentClusters[pos].remove){
				++sizeOfReadVector;
			}
		}
		runningParameters runPars(iteratorMap[i],i, sizeOfReadVector);
		if(verbose_){
			std::cout << std::endl;
			runPars.printIterInfo(std::cout, true);
		}
		if (condensedCollapse_ && runPars.errors_.hqMismatches_ == 0 &&
				runPars.errors_.lqMismatches_ == 0 &&
				runPars.errors_.largeBaseIndel_ == 0 &&
				runPars.errors_.twoBaseIndel_ == 0) {
			auto byCondensed = readVec::organizeByCondensedSeqPositions(currentClusters, positions);
			for (auto& condensedReads : byCondensed) {
				if (regKmer_) {
					collapseWithParametersRegKmer(currentClusters,condensedReads.second,
																										 runPars, alignerObj);
				} else {
					collapseWithParametersQualKmer(currentClusters,condensedReads.second,
																											runPars, alignerObj);
				}
			}
		} else {
			if (regKmer_) {
				collapseWithParametersRegKmer(currentClusters,positions, runPars,
																									 alignerObj);
			} else {
				collapseWithParametersQualKmer(currentClusters,positions, runPars,
																										alignerObj);
			}
		}
		for(const auto & pos : positions){
			currentClusters[pos].updateName();
		}
		auto comp = [&currentClusters](const uint64_t & pos1, const uint64_t & pos2){
      if (currentClusters[pos1].seqBase_.cnt_ == currentClusters[pos2].seqBase_.cnt_) {
        return currentClusters[pos1].averageErrorRate < currentClusters[pos2].averageErrorRate;
      } else {
        return currentClusters[pos1].seqBase_.cnt_ > currentClusters[pos2].seqBase_.cnt_;
      }
		};
		sort(positions, comp);
		for(const auto & pos : positions){
			currentClusters[pos].calculateConsensus(alignerObj, true);
		}
		if (adjustHomopolyerRuns_) {
			for(const auto & pos : positions){
				currentClusters[pos].adjustHomopolyerRunQualities();
			}
		}
	}
}

template <class CLUSTER>
void collapser::collapseWithParametersRegKmer(
    std::vector<CLUSTER> &comparingReads, const runningParameters &runParams,
    aligner &alignerObj) {
  int sizeOfReadVector = readVec::getReadVectorSize(comparingReads);
  if (sizeOfReadVector < 2) {
    return;
  }
  if (verbose_){
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  }
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst_) {
    for (auto &reverseRead : iter::reverse(comparingReads)) {
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchRegKmer(reverseRead, comparingReads,
      		runParams, amountAdded, alignerObj);
    }
  } else {
    for (auto &forwardRead : comparingReads) {
      if (forwardRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchRegKmer(forwardRead, comparingReads, runParams, amountAdded,
                       alignerObj);
    }
  }
  if (verbose_)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
}

template <class CLUSTER>
void collapser::collapseWithParametersRegKmer(
    std::vector<CLUSTER> &comparingReads,
    const std::vector<uint64_t> & positions,
    const runningParameters &runParams,
    aligner &alignerObj) {
  //int sizeOfReadVector = readVec::getReadVectorSize(comparingReads);
	//int sizeOfReadVector = positions.size();
	int sizeOfReadVector = 0;
	for(const auto & pos : positions){
		if(!comparingReads[pos].remove){
			++sizeOfReadVector;
		}
	}
  if (sizeOfReadVector < 2) {
    return;
  }
  if (verbose_){
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  }
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst_) {
    for (const auto &reverseReadPos : iter::reverse(positions)) {
    	auto & reverseRead = comparingReads[reverseReadPos];
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchRegKmer(reverseRead, comparingReads, positions,
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
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchRegKmer(forwardRead, comparingReads,positions,
      		runParams, amountAdded,
                       alignerObj);
    }
  }
  if (verbose_)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
}
template <class CLUSTER>
void collapser::collapseWithParametersRegKmerThreaded(
    std::vector<CLUSTER> &comparingReads, const runningParameters &runParams,
    std::vector<aligner> & aligners, uint32_t numThreads) {
  int sizeOfReadVector = readVec::getReadVectorSize(comparingReads);
  if (sizeOfReadVector < 2) {
    return;
  }
  if (verbose_){
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  }
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst_) {
    for (auto &reverseRead : iter::reverse(comparingReads)) {
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchRegKmerThreaded(reverseRead, comparingReads,
      		runParams, amountAdded, aligners.front());
    }
  } else {
    for (auto &forwardRead : comparingReads) {
      if (forwardRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchRegKmerThreaded(forwardRead, comparingReads, runParams, amountAdded,
      		aligners.front());
    }
  }
  if (verbose_)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
}

template <class CLUSTER>
void collapser::collapseWithParametersQualKmer(
    std::vector<CLUSTER> &comparingReads,
    const runningParameters &runParams,
    aligner &alignerObj) {
  int sizeOfReadVector = readVec::getReadVectorSize(comparingReads);
  if (sizeOfReadVector < 2) {
    return;
  }
  if (verbose_)
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst_) {
    for (auto &reverseRead : iter::reverse(comparingReads)) {
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(reverseRead, comparingReads, runParams, amountAdded,
                        alignerObj);
    }
  } else {
    for (auto &forwardRead : comparingReads) {
      if (forwardRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(forwardRead, comparingReads, runParams, amountAdded,
                        alignerObj);
    }
  }
  if (verbose_)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
}

template <class CLUSTER>
void collapser::collapseWithParametersQualKmer(std::vector<CLUSTER> &comparingReads,
		  																		const std::vector<uint64_t> & positions,
		                                      const runningParameters &runParams,
		                                      aligner &alignerObj) {
  int sizeOfReadVector = readVec::getReadVectorSize(comparingReads);
  if (sizeOfReadVector < 2) {
    return;
  }
  if (verbose_)
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst_) {
    for (const auto &reverseReadPos : iter::reverse(positions)) {
    	auto & reverseRead = comparingReads[reverseReadPos];
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(reverseRead, comparingReads,positions, runParams, amountAdded,
                        alignerObj);
    }
  } else {
    for (const auto &forwardReadPos : positions) {
    	auto & forwardRead = comparingReads[forwardReadPos];
      if (forwardRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(forwardRead, comparingReads,positions, runParams, amountAdded,
                        alignerObj);
    }
  }
  if (verbose_)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
}

template <class CLUSTER>
void collapser::collapseSamples(
    const std::vector<CLUSTER> &clusters,
    std::map<int, std::vector<double>> iteratorMap, const std::string &sortBy,
    const std::string &sampName,
    std::map<std::string, std::vector<CLUSTER>> &samples, aligner &alignerObj) {
  // consider putting default runcutoff as zero so it doesn't do kmer matching
  // during sample collapse
  samples[sampName] =
      collapseCluster(clusters, iteratorMap, sortBy, alignerObj);
}

template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseCluster(
    std::vector<CLUSTER> clusters,
    std::map<int, runningParameters> iteratorMap, const std::string &sortBy,
    aligner &alignerObj){
  // std::cout <<"cc1" << std::endl;
  kmerMaps kMaps = kmerCalculator::indexKmerMpas(clusters, kLength_, runCutOff_,
                                                 qualRunCutOff_);
  // std::cout <<"cc2" << std::endl;
  alignerObj.setKmerMpas(kMaps);
  // std::cout <<"cc3" << std::endl;
  // go over the iterations collapsing down and down, sort the reads after each
  // iteration
  // and calculate a consensus after each one
  readVec::allSetLetterCount(clusters);
  // std::cout <<"cc4" << std::endl;
  for (int i = 1; i <= (int)iteratorMap.size(); ++i) {
    runningParameters runPars(iteratorMap[i]);
    // std::cout <<"cc2" << std::endl;
    if (verbose_) {
      runPars.printIterInfo(std::cout, true);
    }
    if (condensedCollapse_ && runPars.errors_.hqMismatches_ == 0 &&
        runPars.errors_.lqMismatches_ == 0 &&
        runPars.errors_.largeBaseIndel_ == 0 &&
        runPars.errors_.twoBaseIndel_ == 0) {
      auto byCondensed = readVec::organizeByCondensedSeq(clusters);
      for (auto &condensedReads : byCondensed) {
        if (regKmer_) {
          collapser::collapseWithParametersRegKmer(condensedReads.second,
                                                   runPars, alignerObj);
        } else {
          collapser::collapseWithParametersQualKmer(condensedReads.second,
                                                    runPars, alignerObj);
        }
      }
      clusters.clear();
      for (const auto &condensedReads : byCondensed) {
        addOtherVec(clusters, condensedReads.second);
      }
    } else {
      if (regKmer_) {
        collapser::collapseWithParametersRegKmer(clusters, runPars, alignerObj);
      } else {
        collapser::collapseWithParametersQualKmer(clusters, runPars,
                                                  alignerObj);
      }
    }
    readVec::allUpdateName(clusters);
    readVecSorter::sortReadVector(clusters, sortBy);
    clusterVec::allCalculateConsensus(clusters, alignerObj, true);
  }
  return readVecSplitter::splitVectorOnRemove(clusters).first;
}

template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseCluster(
    std::vector<CLUSTER> clusters,
    std::map<int, std::vector<double>> iteratorMap, const std::string &sortBy,
    aligner &alignerObj) {
	std::map<int, runningParameters> iteratorMapPars;
	for(const auto & iteration : iteratorMap){
		iteratorMapPars[iteration.first] = runningParameters(iteration.second,iteration.first,
				readVec::getTotalReadCount(clusters));

	}
  return collapseCluster(clusters, iteratorMapPars, sortBy, alignerObj);
}

template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseOneOffsQualKmer(
    std::vector<CLUSTER> clusters, runningParameters runParams,
    double fracCutOff, aligner &alignerObj) {
  runParams.errors_.hqMismatches_ = 1;
  kmerMaps kMaps = kmerCalculator::indexKmerMpas(clusters, kLength_, runCutOff_,
                                                 qualRunCutOff_);
  // std::cout <<"cc2" << std::endl;
  alignerObj.setKmerMpas(kMaps);
  // std::cout <<"cc3" << std::endl;
  // go over the iterations collapsing down and down, sort the reads after each
  // iteration
  // and calculate a consensus after each one

  // std::cout <<"cc4" << std::endl;
  int sizeOfReadVector = readVec::getReadVectorSize(clusters);
  if (verbose_)
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst_) {
    for (auto &reverseRead : iter::reverse(clusters)) {
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
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
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(forwardRead, clusters, runParams, amountAdded,
                        alignerObj);
    }
  }
  if (verbose_)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
  // collapser::collapseWithParameters(clusters, runPars, alignerObj);
  readVec::allUpdateName(clusters);
  readVecSorter::sortReadVector(clusters, "totalCount");
  return readVecSplitter::splitVectorOnRemove(clusters).first;
}

template <class CLUSTER>
std::vector<CLUSTER> collapser::collapseOneOffsRegKmer(
    std::vector<CLUSTER> clusters, runningParameters runParams,
    double fracCutOff, aligner &alignerObj) {
  runParams.errors_.hqMismatches_ = 1;
  kmerMaps kMaps = kmerCalculator::indexKmerMpas(clusters, kLength_, runCutOff_,
                                                 qualRunCutOff_);
  // std::cout <<"cc2" << std::endl;
  alignerObj.setKmerMpas(kMaps);
  // std::cout <<"cc3" << std::endl;
  // go over the iterations collapsing down and down, sort the reads after each
  // iteration
  // and calculate a consensus after each one

  // std::cout <<"cc4" << std::endl;
  int sizeOfReadVector = readVec::getReadVectorSize(clusters);
  if (verbose_)
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst_) {
    for (auto &reverseRead : iter::reverse(clusters)) {
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose_ && clusterCounter % 100 == 0) {
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
      if (verbose_ && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      findMatchQualKmer(forwardRead, clusters, runParams, amountAdded,
                        alignerObj);
    }
  }
  if (verbose_)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
  // collapser::collapseWithParameters(clusters, runPars, alignerObj);
  readVec::allUpdateName(clusters);
  readVecSorter::sortReadVector(clusters, "totalCount");
  return readVecSplitter::splitVectorOnRemove(clusters).first;
}

template <class CLUSTER>
void collapser::findMatchQualKmer(CLUSTER &read,
                                  std::vector<CLUSTER> &comparingReads,
                                  const runningParameters &runParams,
                                  size_t &amountAdded, aligner &alignerObj) {
  int count = -1;
  double bestScore = 0;
  bool foundMatch = false;
  uint32_t bestSearching = 0;
  uint32_t bestClusterPos = std::numeric_limits<uint32_t>::max();
  for (const auto &clusPos : iter::range(comparingReads.size())) {
  	auto & clus = comparingReads[clusPos];
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      //continue;
      break;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 10) {
        continue;
      }
    }

    if (useReadLen_){
    	++count;
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > readLenDiff_) {
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

    if (1 + count > runParams.stopCheck_ || bestSearching > bestMatchCheck_) {
      break;
    }

    bool matching = false;
    if (clus.previousErrorChecks.find(read.firstReadName) !=
            clus.previousErrorChecks.end() &&
        read.previousErrorChecks.find(clus.firstReadName) !=
            read.previousErrorChecks.end()) {
      matching = runParams.errors_.passErrorProfileLowKmer(
          read.previousErrorChecks.at(clus.firstReadName));
      alignerObj.alignVec(clus, read, local_);
      //alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
    } else {
      alignerObj.alignVec(clus, read, local_);
      alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
			if(alignerObj.distances_.queryCoverage_ < 0.50){
				matching = false;
			}else{
				errorProfile currentProfile = alignerObj.checkAlignmentBothQualKmer(
				          clus, read, runParams, checkKmers_, kmersByPosition_,
				          weighHomopolyer_);
				read.previousErrorChecks[clus.firstReadName] = currentProfile;
				clus.previousErrorChecks[read.firstReadName] = currentProfile;
				matching = runParams.errors_.passErrorProfileLowKmer(currentProfile);
			}
		}

    if (matching) {
      foundMatch = true;
      if (findingBestMatch_) {
      	//std::cout << "score: " << alignerObj.score_ << std::endl;
      	double currentScore = 0;
      	if(eventBased){
      		alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
      		currentScore = alignerObj.distances_.eventBasedIdentity_;
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
    if (findingBestMatch_) {
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
void collapser::findMatchQualKmer(CLUSTER &read,
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
  	auto & clus = comparingReads[clusPos];
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      //continue;
      break;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 10) {
        continue;
      }
    }

    if (useReadLen_){
    	++count;
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > readLenDiff_) {
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

    if (1 + count > runParams.stopCheck_ || bestSearching > bestMatchCheck_) {
      break;
    }

    bool matching = false;
    if (clus.previousErrorChecks.find(read.firstReadName) !=
            clus.previousErrorChecks.end() &&
        read.previousErrorChecks.find(clus.firstReadName) !=
            read.previousErrorChecks.end()) {
      matching = runParams.errors_.passErrorProfileLowKmer(
          read.previousErrorChecks.at(clus.firstReadName));
      alignerObj.alignVec(clus, read, local_);
      //alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
    } else {
      alignerObj.alignVec(clus, read, local_);
      alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
			if(alignerObj.distances_.queryCoverage_ < 0.50){
				matching = false;
			}else{
				errorProfile currentProfile = alignerObj.checkAlignmentBothQualKmer(
				          clus, read, runParams, checkKmers_, kmersByPosition_,
				          weighHomopolyer_);
				read.previousErrorChecks[clus.firstReadName] = currentProfile;
				clus.previousErrorChecks[read.firstReadName] = currentProfile;
				matching = runParams.errors_.passErrorProfileLowKmer(currentProfile);
			}
		}

    if (matching) {
      foundMatch = true;
      if (findingBestMatch_) {
      	//std::cout << "score: " << alignerObj.score_ << std::endl;

      	double currentScore = 0;
      	if(eventBased){
      		alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
      		currentScore = alignerObj.distances_.eventBasedIdentity_;
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
    if (findingBestMatch_) {
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
void collapser::findMatchRegKmer(CLUSTER &read,
                                 std::vector<CLUSTER> &comparingReads,
                                 const runningParameters &runParams,
                                 size_t &amountAdded, aligner &alignerObj) {
  int count = -1;
  double bestScore = 0;
  bool foundMatch = false;
  uint32_t bestSearching = 0;
  uint32_t bestClusterPos = std::numeric_limits<uint32_t>::max();
  for (const auto &clusPos : iter::range(comparingReads.size())) {
  	auto & clus = comparingReads[clusPos];
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      //continue;
      break;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 10) {
        continue;
      }
    }

    if (useReadLen_){
    	++count;
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > readLenDiff_) {
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

    if (1 + count > runParams.stopCheck_ || bestSearching > bestMatchCheck_) {
      break;
    }

    bool matching = false;
    if (clus.previousErrorChecks.find(read.firstReadName) !=
            clus.previousErrorChecks.end() &&
        read.previousErrorChecks.find(clus.firstReadName) !=
            read.previousErrorChecks.end()) {
      matching = runParams.errors_.passErrorProfileLowKmer(
          read.previousErrorChecks.at(clus.firstReadName));
      alignerObj.alignVec(clus, read, local_);
      //alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
    } else {
      alignerObj.alignVec(clus, read, local_);
      alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
			if(alignerObj.distances_.queryCoverage_ < 0.50){
				matching = false;
			}else{
				errorProfile currentProfile = alignerObj.checkAlignmentBothRegKmer(
				          clus, read, runParams, checkKmers_, kmersByPosition_,
				          weighHomopolyer_);
				read.previousErrorChecks[clus.firstReadName] = currentProfile;
				clus.previousErrorChecks[read.firstReadName] = currentProfile;
				matching = runParams.errors_.passErrorProfileLowKmer(currentProfile);
			}
		}

    if (matching) {
      foundMatch = true;
      if (findingBestMatch_) {
      	//std::cout << "score: " << alignerObj.score_ << std::endl;

      	double currentScore = 0;
      	if(eventBased){
      		alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
      		currentScore = alignerObj.distances_.eventBasedIdentity_;
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
    if (findingBestMatch_) {
    	if(bestClusterPos != std::numeric_limits<uint32_t>::max()){
    		comparingReads[bestClusterPos].addRead(read);
    		//comparingReads[bestClusterPos].allInputClusters.emplace_back(read);
        read.remove = true;
        ++amountAdded;
    	}
    }
  }
  /*
	if(read.seqBase_.name_ == "UNC20:154:000000000-A5EFK:1:2109:4118:18405 1:N:0:TAAGGCGTAGATCG_t1"){
		std::cout << read.seqBase_.name_ << std::endl;
		std::cout << read.seqBase_.seq_ << std::endl;
		std::cout << vectorToString(read.seqBase_.qual_) << std::endl;
	}*/
}

template <class CLUSTER>
void collapser::findMatchRegKmer(CLUSTER &read,
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
  	auto & clus = comparingReads[clusPos];
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      //continue;
      break;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 10) {
        continue;
      }
    }

    if (useReadLen_){
    	++count;
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > readLenDiff_) {
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

    if (1 + count > runParams.stopCheck_ || bestSearching > bestMatchCheck_) {
      break;
    }

    bool matching = false;
    if (clus.previousErrorChecks.find(read.firstReadName) !=
            clus.previousErrorChecks.end() &&
        read.previousErrorChecks.find(clus.firstReadName) !=
            read.previousErrorChecks.end()) {
      matching = runParams.errors_.passErrorProfileLowKmer(
          read.previousErrorChecks.at(clus.firstReadName));
      alignerObj.alignVec(clus, read, local_);
      //alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
    } else {
      alignerObj.alignVec(clus, read, local_);
      alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
			if(alignerObj.distances_.queryCoverage_ < 0.50){
				matching = false;
			}else{
				errorProfile currentProfile = alignerObj.checkAlignmentBothRegKmer(
				          clus, read, runParams, checkKmers_, kmersByPosition_,
				          weighHomopolyer_);
				read.previousErrorChecks[clus.firstReadName] = currentProfile;
				clus.previousErrorChecks[read.firstReadName] = currentProfile;
				matching = runParams.errors_.passErrorProfileLowKmer(currentProfile);
			}
		}

    if (matching) {
      foundMatch = true;
      if (findingBestMatch_) {

      	double currentScore = 0;
      	if(eventBased){
      		alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
      		currentScore = alignerObj.distances_.eventBasedIdentity_;
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
    if (findingBestMatch_) {
    	if(bestClusterPos != std::numeric_limits<uint32_t>::max()){
    		comparingReads[bestClusterPos].addRead(read);
    		//comparingReads[bestClusterPos].allInputClusters.emplace_back(read);
        read.remove = true;
        ++amountAdded;
    	}
    }
  }
  /*
	if(read.seqBase_.name_ == "UNC20:154:000000000-A5EFK:1:2109:4118:18405 1:N:0:TAAGGCGTAGATCG_t1"){
		std::cout << read.seqBase_.name_ << std::endl;
		std::cout << read.seqBase_.seq_ << std::endl;
		std::cout << vectorToString(read.seqBase_.qual_) << std::endl;
	}*/
}

template <class CLUSTER>
void collapser::findMatchRegKmerThreaded(CLUSTER &read,
                                 std::vector<CLUSTER> &comparingReads,
                                 const runningParameters &runParams,
                                 size_t &amountAdded, aligner &alignerObj) {
	std::lock_guard<std::mutex> lockRead(read.mtx_);
  int count = -1;
  double bestScore = 0;
  int bestSearching = 0;
  bool foundMatch = false;
  uint64_t bestCluster = -1;
  for (auto &clusPos : iter::range(comparingReads.size())) {
  	auto & clus = comparingReads[clusPos];
  	std::lock_guard<std::mutex> lockCompareRead(clus.mtx_);
    if (foundMatch) {
      ++bestSearching;
    }
    if (clus.remove) {
      continue;
    } else {
      ++count;
    }
    if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      continue;
    }
    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }
    if (skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (1 + count > runParams.stopCheck_ || bestSearching > bestMatchCheck_) {
      break;
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(read.seqBase_.seq_.length() -
                   clus.seqBase_.seq_.length()) > 10) {
        continue;
      }
    }
    bool matching = false;
    if (clus.previousErrorChecks.find(read.firstReadName) !=
            clus.previousErrorChecks.end() &&
        read.previousErrorChecks.find(clus.firstReadName) !=
            read.previousErrorChecks.end()) {
      matching = runParams.errors_.passErrorProfileLowKmer(
          read.previousErrorChecks.at(clus.firstReadName));
    } else {
      alignerObj.alignVec(clus, read, local_);
      //alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
			if(alignerObj.distances_.queryCoverage_ < 0.50){
				matching = false;
			}else{
				errorProfile currentProfile = alignerObj.checkAlignmentBothRegKmer(
				          clus, read, runParams, checkKmers_, kmersByPosition_,
				          weighHomopolyer_);
				read.previousErrorChecks[clus.firstReadName] = currentProfile;
				clus.previousErrorChecks[read.firstReadName] = currentProfile;
				matching = runParams.errors_.passErrorProfileLowKmer(currentProfile);
			}
		}
    if (matching) {
      foundMatch = true;
      if (findingBestMatch_) {
        if (alignerObj.parts_.score_ > bestScore) {
          bestScore = alignerObj.parts_.score_;
          bestCluster = clusPos;
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
    if (findingBestMatch_) {
    	std::lock_guard<std::mutex> lockCompareRead(comparingReads[bestCluster].mtx_);
    	comparingReads[bestCluster].addRead(read);
    	comparingReads[bestCluster].allInputClusters.push_back(read);
      read.remove = true;
      ++amountAdded;
    }
  }
}
template <class CLUSTER>
void collapser::findMatchOneOff(CLUSTER &read,
                                std::vector<CLUSTER> &comparingReads,
                                const runningParameters &runParams,
                                double fracCutoff, size_t &amountAdded,
                                aligner &alignerObj) {
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
    if (skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > fractionDifferenceCutOff_) {
        continue;
      }
    }
    if (alignerObj.CountEndGaps()) {
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > 10) {
        continue;
      }
    }

    if (useReadLen_){
    	++count;
      if (std::abs(static_cast<int32_t>(read.seqBase_.seq_.length()) -
      		static_cast<int32_t>(clus.seqBase_.seq_.length())) > readLenDiff_) {
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

    if (1 + count > runParams.stopCheck_ || bestSearching > bestMatchCheck_) {
      break;
    }
    bool matching = false;
    if (clus.previousErrorChecks.find(read.firstReadName) !=
            clus.previousErrorChecks.end() &&
        read.previousErrorChecks.find(clus.firstReadName) !=
            read.previousErrorChecks.end()) {
      matching = runParams.errors_.passErrorProfileLowKmer(
          read.previousErrorChecks.at(clus.firstReadName));
      alignerObj.alignVec(clus, read, local_);
      //alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
    } else {
      alignerObj.alignVec(clus, read, local_);
      alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
      errorProfile currentProfile = alignerObj.checkAlignmentBothQualKmer(
          clus, read, runParams, checkKmers_, kmersByPosition_,
          weighHomopolyer_);
      read.previousErrorChecks[clus.firstReadName] = currentProfile;
      clus.previousErrorChecks[read.firstReadName] = currentProfile;
      matching = runParams.errors_.passErrorProfileLowKmer(currentProfile);
    }

    if (matching) {
      foundMatch = true;
      if (findingBestMatch_) {
      	double currentScore = 0;
      	if(eventBased){
      		alignerObj.profilePrimerAlignment(clus, read, weighHomopolyer_);
      		currentScore = alignerObj.distances_.eventBasedIdentity_;
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
    if (findingBestMatch_) {
      read.remove = true;
      comparingReads[bestClusterPos].addRead(read);
     // readVec::getReadByName(comparingReads, bestCluster).addRead(read);
      ++amountAdded;
    }
  }
}
}  // namespace bib
#ifndef NOT_HEADER_ONLY
#include "collapser.cpp"
#endif