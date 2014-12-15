#pragma once
//
//  clusterCollapser.h
//  sequenceTools
//
//  Created by Nick Hathaway on 11/16/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils/utils.hpp"
#include "bibseq/objects/seqObjects.h"
#include "bibseq/objects/helperObjects/kmer.hpp"
#include "bibseq/alignment.h"
#include "bibseq/helpers/kmerCalculator.hpp"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/IO/readObjectIO.hpp"
#include "bibseq/seqToolsUtils.h"
#include "bibseq/readVectorManipulation.h"

#include <cppitertools/reverse.hpp>

namespace bibseq {

class clusterCollapser {

 public:
  clusterCollapser() {}
  template<typename T>
  static std::vector<identicalCluster> collapseToUniqueReads(
        const std::vector<T> &input, const std::string &repQual);

  // function to collpase to identical reads
  static std::vector<identicalCluster> collapseIdenticalReads(
      const std::vector<readObject> &seqList, const std::string &repQual,
      const std::string &lower);

  template <class CLUSTER>
  static std::vector<uint> findMatchIndex(
      const CLUSTER &read, const std::vector<CLUSTER> &comparingReads,
      const runningParameters &runParams, aligner &alignerObj,
      size_t &amountAdded, bool findingBestMatch, int bestMatchCheck,
      kmerMaps &kMaps, bool kmersByPos, int runCutOff, int kLength,
      bool usingQuality, bool local, bool weighHomopolyer,
      bool skipOnLetterCounterDifference, double fractionDifferenceCutOff);

  template <class CLUSTER>
  static void collapseWithParametersByIndices(
      std::vector<CLUSTER> &comparingReads, aligner &alignerObj,
      const runningParameters &runParams, bool findingBestMatch,
      int bestMatchCheck, bool local, bool usingQuality, kmerMaps &kMaps,
      bool kmersByPosition, int runCutOff, bool verbose, int kLength,
      bool smallestFirst, bool weighHomopolyer,
      bool skipOnLetterCounterDifference, double fractionDifferenceCutOff);

  // mark chimeras
  static std::vector<cluster> markChimeras(std::vector<cluster> &processedReads,
                                           aligner &alignerObj,
                                           double parentFreqs, bool kmersByPos,
                                           int kLength, int runCutOff,
                                           bool local, bool weighHomopolyer);

  static void markChimerasAdvanced(std::vector<cluster> &processedReads,
                                   aligner &alignerObj, double parentFreqs,
                                   int runCutOff, bool local,
                                   const errorProfile &chiOverlap,
                                   uint32_t overLapSizeCutoff,
                                   bool weighHomopolyer, uint32_t &chimeraCount,
                                   uint32_t allowableError);

  // cluster down on tandems
  static void collapseTandems(std::vector<cluster> &processedReads,
                              aligner &alignerObj, int runCutOff, int kLength,
                              bool kMersByPosition, double freqCutoff,
                              bool local, bool weighHomopolyer);
};
template<typename T>
std::vector<identicalCluster> clusterCollapser::collapseToUniqueReads(
        const std::vector<T> &input, const std::string &repQual){
	std::vector<identicalCluster> output;
	std::unordered_map<std::string, std::vector<T>> clusters;
	for(const auto & read : input){
		clusters[read.seqBase_.seq_].emplace_back(read);
	}
	for(const auto & clus : clusters){
		output.emplace_back(identicalCluster(clus.second, repQual));
	}
	return output;
}
template <class CLUSTER>
std::vector<uint> clusterCollapser::findMatchIndex(
    const CLUSTER &read, const std::vector<CLUSTER> &comparingReads,
    const runningParameters &runParams, aligner &alignerObj,
    size_t &amountAdded, bool findingBestMatch, int bestMatchCheck,
    kmerMaps &kMaps, bool kmersByPos, int runCutOff, int kLength,
    bool usingQuality, bool local, bool weighHomopolyer,
    bool skipOnLetterCounterDifference, double fractionDifferenceCutOff) {
  int count = -1;
  uint32_t pos = 0;
  double bestScore = 0;
  int bestSearching = 0;
  bool foundMatch = false;
  bool beginning = true;
  std::vector<uint> indices;
  for (auto &clus : comparingReads) {
    if (!beginning) {
      ++pos;
    } else {
      beginning = false;
    }
    if (foundMatch) {
      ++bestSearching;
    }
    // std::cout<<"mark 1"<<std::endl;
    if (clus.remove) {
      continue;
    } else {
      ++count;
    }
    // std::cout<<"mark 2"<<std::endl;
    if (clus.totalCount <= runParams.smallCheckStop_) {
      continue;
    }
    // std::cout<<"mark 3"<<std::endl;
    if (clus.name == read.name) {
      continue;
    }
    if (skipOnLetterCounterDifference) {
      // std::cout<<"skipping on letter counter difference"<<std::endl;
      double sum = 0;
      for (const auto &count : read.counter.fractions) {
        sum += clus.counter.fractions.at(count.first) - count.second;
      }
      if (sum > fractionDifferenceCutOff) {
        continue;
      }
    }
    // std::cout<<"mark 4"<<std::endl;
    /*std::cout<<"\t"<<count<<std::endl;
     std::cout<<"\t"<<runParams.stopCheck<<std::endl;
     std::cout<<"\t"<<bestSearching<<std::endl;
     std::cout<<"\t"<<bestMatchCheck<<std::endl;*/
    if (1 + count > runParams.stopCheck_ || bestSearching > bestMatchCheck) {
      // std::cout<<"mark breaking on stop check or bestMatchCheck"<<std::endl;
      break;
    }
    // std::cout<<"mark 5"<<std::endl;
    if (alignerObj.CountEndGaps()) {
      if (abs((int)read.seq.length() - (int)clus.seq.length()) > 10) {
        continue;
      }
    }
    bool matching = false;
    alignerObj.alignVec(clus, read, local);
    matching = alignerObj.checkAlignmentBool(clus, read, runParams, kLength,
                                             kmersByPos, usingQuality, false,
                                             weighHomopolyer);
    if (matching) {
      foundMatch = true;
      if (findingBestMatch) {
        if (alignerObj.parts_.score_ > bestScore) {
          bestScore = alignerObj.parts_.score_;
          indices.clear();
          indices.push_back(pos);
        } else if (alignerObj.parts_.score_ == bestScore) {
          indices.push_back(pos);
        }
      } else {
        indices.push_back(pos);
        return indices;
        break;
      }
    }
  }
  return indices;
}
template <class CLUSTER>
void clusterCollapser::collapseWithParametersByIndices(
    std::vector<CLUSTER> &comparingReads, aligner &alignerObj,
    const runningParameters &runParams, bool findingBestMatch,
    int bestMatchCheck, bool local, bool usingQuality, kmerMaps &kMaps,
    bool kmersByPosition, int runCutOff, bool verbose, int kLength,
    bool smallestFirst, bool weighHomopolyer,
    bool skipOnLetterCounterDifference, double fractionDifferenceCutOff) {
  int sizeOfReadVector = getReadVectorSize(comparingReads);
  if (sizeOfReadVector < 2) {
    return;
  }
  if (verbose)
    std::cout << "Starting with " << sizeOfReadVector << " clusters"
              << std::endl;
  int clusterCounter = 0;
  size_t amountAdded = 0;
  if (smallestFirst) {
    for (auto &reverseRead : iter::reverse(comparingReads)) {
      if (reverseRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }
      if (verbose && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }

      std::vector<uint> indices = findMatchIndex(
          reverseRead, comparingReads, runParams, alignerObj, amountAdded,
          findingBestMatch, bestMatchCheck, kMaps, kmersByPosition, runCutOff,
          kLength, usingQuality, local, weighHomopolyer,
          skipOnLetterCounterDifference, fractionDifferenceCutOff);
      if (indices.size() > 0) {
        ++amountAdded;
        comparingReads[indices[0]].addRead(reverseRead);
        comparingReads[indices[0]].allInputClusters.push_back(reverseRead);
        reverseRead.remove = true;
      }
    }
  } else {
    for (auto &forwardRead : comparingReads) {
      if (forwardRead.remove) {
        continue;
      } else {
        clusterCounter++;
      }

      if (verbose && clusterCounter % 100 == 0) {
        std::cout << "Currently on cluster " << clusterCounter << " of "
                  << sizeOfReadVector << std::endl;
      }
      std::vector<uint> indices = findMatchIndex(
          forwardRead, comparingReads, runParams, alignerObj, amountAdded,
          findingBestMatch, bestMatchCheck, kMaps, kmersByPosition, runCutOff,
          kLength, usingQuality, local, weighHomopolyer,
          skipOnLetterCounterDifference, fractionDifferenceCutOff);
      if (indices.size() > 0) {
        ++amountAdded;
        comparingReads[indices[0]].addRead(forwardRead);
        comparingReads[indices[0]].allInputClusters.push_back(forwardRead);
        forwardRead.remove = true;
      }
    }
  }
  if (verbose)
    std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
              << " clusters" << std::endl;
}

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "clusterCollapser.cpp"
#endif
