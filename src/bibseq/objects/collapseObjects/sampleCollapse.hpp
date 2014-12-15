#pragma once
//
//  sampleCollapse.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "collapseCommon.hpp"

namespace bibseq {
namespace collapse {
class sampleCollapse {

 public:
  // constructor,
  // template<typename T>
  sampleCollapse(const std::vector<std::vector<sampleCluster>> &inputClusters,
                 const std::string &sampName, int sizeCutOff)
      : sampName_(sampName) {
    uint32_t singletCount = 0;
    for (const auto &reads : inputClusters) {
      addOtherVec(input_.clusters_, reads);
      auto tempVec = readVecSplitter::splitVectorOnReadCountAdd(
          reads, sizeCutOff, excluded_.clusters_, singletCount, false);
      clusterVec::allSetFractionClusters(tempVec);
      addOtherVec(collapsed_.clusters_, tempVec);
    }
    updateInitialInfos(false);
  }
  sampleCollapse() {}
  // memebers
  // name
  std::string sampName_;
  // the initial clusters of the sample
  clusterSet input_;
  // the excluded clusters
  clusterSet excluded_;
  // collapsed clusters
  clusterSet collapsed_;
  // functions
  // collapse the input clusters
  void cluster(collapser &collapserObj,
               std::map<int, std::vector<double>> iteratorMap,
               const std::string &sortBy, aligner &alignerObj);
  // excludes
  void excludeChimeras(bool update);
  void excludeFraction(double fractionCutOff, bool update);
  void excludeBySampNum(uint64_t sampsRequired, bool update);
  //
  void renameClusters(const std::string &sortBy);
  // update the exclusioninfos
  void updateExclusionInfos(bool clearCurrentInfos) {
    excluded_.updateSetInfo(clearCurrentInfos);
    sampleCluster::updateAllClusters(excluded_.clusters_, input_.infos_);
  }
  // update the initial infos
  void updateInitialInfos(bool clearCurrentInfos) {
    input_.updateSetInfo(clearCurrentInfos);
  }
  // update the collapsed infos
  void updateCollapsedInfos(bool clearCurrentInfos) {
    collapsed_.updateSetInfo(clearCurrentInfos);
    sampleCluster::updateAllClusters(collapsed_.clusters_, collapsed_.infos_);
  }
  void updateAfterExclustion(bool clearCurrentInfos) {
    updateExclusionInfos(clearCurrentInfos);
    updateCollapsedInfos(clearCurrentInfos);
  }
  // output reads
  std::vector<sampleCluster> createOutput(bool renameFirst,
                                          const std::string &sortBy);
  // write
  void writeExcluded(const std::string &outDirectory,
                     const std::string &outFormat, bool overWrite,
                     bool exitOnFailureToWrite) const;
  void writeInitial(const std::string &outDirectory,
                    const std::string &outFormat, bool overWrite,
                    bool exitOnFailureToWrite) const;
  void writeFinal(const std::string &outDirectory, const std::string &outFormat,
                  bool overWrite, bool exitOnFailureToWrite) const;
  void writeFinalOrginals(const std::string &outDirectory) const;
  // output info
  std::string getSimpleSampInfo(const std::string &delim = "\t") const;
  VecStr getAllInfoVec(bool checkingExpected,
                       const std::string &delim = "\t") const;
  std::map<std::string, std::string, std::greater<std::string>> getAllInfoMap(
      bool checkingExpected, const std::string &delim = "\t") const;
};
}  // namspace collpase
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "sampleCollapse.cpp"
#endif
