#pragma once
//
//  populationCollapse.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/collapseObjects/sampleCollapse.hpp"

namespace bibseq {
namespace collapse {
class populationCollapse {

 public:
  populationCollapse(const std::vector<sampleCluster> &inputClusters,
                     const std::string &populationName)
      : input_(clusterSet(inputClusters)), populationName_(populationName) {}
  populationCollapse() {}

  // members
  // the initial clusters of the sample
  clusterSet input_;
  // collapsed clusters
  clusterSet collapsed_;

  std::string populationName_;

  // functions
  void popCluster(collapser &collapserObj,
                  std::map<int, std::vector<double>> iteratorMap,
                  const std::string &sortBy, aligner &alignerObj);

  // update the initial infos
  void updateInitialInfos(bool clearCurrentInfos) {
    input_.updateSetInfo(clearCurrentInfos);
  }
  // update the collapsed infos
  void updateCollapsedInfos(bool clearCurrentInfos) {
    collapsed_.updateSetInfo(clearCurrentInfos);
    collapsed_.updateSubClustersPositions(clearCurrentInfos);
  }
  void renameToOtherPopNames(const std::vector<readObject> &previousPop,
                             aligner &alignerObj);
  void renameClusters(const std::string &sortBy);
  void writeFinal(const std::string &outDirectory, const std::string &outFormat,
                  bool overWrite, bool exitOnFailureToWrite) const;
  void writeFinalInitial(const std::string &outDirectory) const;
};
}  // namspace collpase
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "populationCollapse.cpp"
#endif
