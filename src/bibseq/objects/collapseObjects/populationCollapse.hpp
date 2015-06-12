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
                     const std::string &populationName);

  populationCollapse();
  // members
  // the initial clusters of the sample
  clusterSet input_;
  // collapsed clusters
  clusterSet collapsed_;

  std::string populationName_;

  uint32_t numOfSamps()const;
  // clustering
  void popCluster(collapser &collapserObj,
                  std::map<int, std::vector<double>> iteratorMap,
                  const std::string &sortBy, aligner &alignerObj);
  void popClusterOnId(collapser &collapserObj,
                  std::map<int, std::vector<double>> iteratorMap,
                  const std::string &sortBy, aligner &alignerObj);


  // update the initial infos
  void updateInitialInfos();
  // update the collapsed infos
  void updateCollapsedInfos();

  void renameToOtherPopNames(const std::vector<readObject> &previousPop,
                             aligner &alignerObj);
  void renameClusters();

  void updateInfoWithSampCollapses(const std::map<std::string, sampleCollapse> & sampCollapses);
  void updateInfoWithSampCollapse(const sampleCollapse & sampCollapses);

  //io
  void writeFinal(const std::string &outDirectory, const std::string &outFormat,
                  bool overWrite, bool exitOnFailureToWrite) const;
  void writeFinalInitial(const std::string &outDirectory, const readObjectIOOptions & ioOptions) const;




};
}  // namspace collpase
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "populationCollapse.cpp"
#endif
