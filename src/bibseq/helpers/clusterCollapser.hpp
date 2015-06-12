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





  // mark chimeras
  static std::vector<cluster> markChimeras(std::vector<cluster> &processedReads,
                                           aligner &alignerObj,
                                           double parentFreqs, bool kmersByPos,
                                           int kLength, int runCutOff,
                                           bool local, bool weighHomopolyer);

  static void markChimerasAdvanced(std::vector<cluster> &processedReads,
                                   aligner &alignerObj, double parentFreqs,
                                   int runCutOff, bool local,
                                   const comparison &chiOverlap,
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



}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "clusterCollapser.cpp"
#endif
