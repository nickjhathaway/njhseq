#pragma once
//
//  collapseCommon.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/sampleCluster.hpp"
#include "bibseq/helpers/clusterCollapser.hpp"
#include "bibseq/objects/collapseObjects/collapser.hpp"
#include "bibseq/readVectorManipulation.h"
#include "bibseq/seqToolsUtils.h"
#include "bibseq/helpers/profiler.hpp"

namespace bibseq {
namespace collapse {
// update any info
template <typename T>
void updateInfo(const T& read, std::map<std::string, sampInfo>& infos) {
  if (infos.find(read.sampName) == infos.end()) {
    infos.emplace(read.sampName, sampInfo(read));
  } else {
    infos[read.sampName].update(read);
  }
}

class clusterSet {
 public:
  // constructors
  clusterSet() {}
  clusterSet(const std::vector<sampleCluster>& clusters)
      : clusters_(clusters) {}
  // members
  std::vector<sampleCluster> clusters_;
  std::map<std::string, sampInfo> infos_;
  std::unordered_map<std::string, uint32_t> subClustersPositions_;

  double totalReadCount_;
  int numberOfClusters_;
  // functions
  // update subClustersPositions
  void updateSubClustersPositions(bool clearFirst);
  void updateSetInfo(bool clearCurrentInfos);
  void updateCounts();
  template <typename REF>
  void checkAgainstExpected(const std::vector<REF>& refSeqs,
                            aligner& alignerObj, bool local,
                            bool weighHomopolyers) {
    bool eventBased = true;
  	for (auto& clus : clusters_) {

    	std::string bestRef = profiler::getBestRef(refSeqs, clus, alignerObj,
    			local, weighHomopolyers, eventBased , false, ",");
      auto bestRefs = profiler::compareToRefSingle(
          refSeqs, clus, alignerObj, local, false, weighHomopolyers);
      clus.expectsString = bestRefs.front();
    }
  }
  // writing
  void writeClusters(std::string filename, std::string format, bool overWrite,
                     bool exitOnFailure) const;
};

}  // namspace collapse
}  // namspace bib

#ifndef NOT_HEADER_ONLY
#include "collapseCommon.cpp"
#endif
