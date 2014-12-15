
//
//  sampleCollapse.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/31/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "sampleCollapse.hpp"

namespace bibseq {
namespace collapse {
void sampleCollapse::cluster(collapser &collapserObj,
                             std::map<int, std::vector<double>> iteratorMap,
                             const std::string &sortBy, aligner &alignerObj) {
  // std::cout <<"clus 1 " << std::endl;
  collapsed_.clusters_ = collapserObj.collapseCluster(
      collapsed_.clusters_, iteratorMap, sortBy, alignerObj);
  // std::cout <<"clus 2 " << std::endl;
}
// excludes
void sampleCollapse::excludeChimeras(bool update) {
  for (auto &clus : collapsed_.clusters_) {
    if (clusterVec::isClusterAtLeastHalfChimeric(clus.reads_)) {
      clus.seqBase_.markAsChimeric();
      clus.remove = true;
    }
  }
  uint32_t chimeraNum = 0;
  collapsed_.clusters_ = readVecSplitter::splitVectorOnRemoveAdd(
      collapsed_.clusters_, excluded_.clusters_, chimeraNum, "none", false);
  if (update) {
    updateAfterExclustion(true);
  }
}
void sampleCollapse::excludeFraction(double fractionCutOff, bool update) {
  uint32_t fractionCutOffNum = 0;
  collapsed_.clusters_ = readVecSplitter::splitVectorOnReadFractionAdd(
      collapsed_.clusters_, fractionCutOff, excluded_.clusters_,
      fractionCutOffNum, false);
  if (update) {
    updateAfterExclustion(true);
  }
}
void sampleCollapse::excludeBySampNum(uint64_t sampsRequired, bool update) {
  uint32_t splitCount = 0;
  for (auto &clus : collapsed_.clusters_) {
    if (clus.sampInfos_.size() < sampsRequired) {
      clus.remove = true;
    }
  }
  collapsed_.clusters_ = readVecSplitter::splitVectorOnRemoveAdd(
      collapsed_.clusters_, excluded_.clusters_, splitCount, "none", false);
  if (update) {
    updateAfterExclustion(true);
  }
}
void sampleCollapse::renameClusters(const std::string &sortBy) {
  readVecSorter::sortReadVector(collapsed_.clusters_, sortBy);
  renameReadNames(collapsed_.clusters_, sampName_, true, false);
  for (auto &clus : collapsed_.clusters_) {
    if (clusterVec::isClusterAtLeastHalfChimeric(clus.reads_)) {
      clus.seqBase_.markAsChimeric();
      //clus.remove = true;
    }
  }
  readVec::allUpdateName(collapsed_.clusters_);
}

// output
std::string sampleCollapse::getSimpleSampInfo(const std::string &delim) const {
  // sampName\treadsUsed(%total)\thapsUsed
  return (sampName_ + delim + getPercentageString(collapsed_.totalReadCount_,
                                                  input_.totalReadCount_) +
          delim + std::to_string(collapsed_.numberOfClusters_));
}
VecStr sampleCollapse::getAllInfoVec(bool checkingExpected,
                                     const std::string &delim) const {
  auto mapOfInfos = getAllInfoMap(checkingExpected, delim);
  VecStr ans = getVectorOfMapValues(mapOfInfos);
  return ans;
}
std::map<std::string, std::string, std::greater<std::string>>
sampleCollapse::getAllInfoMap(bool checkingExpected,
                              const std::string &delim) const {
  std::map<std::string, std::string, std::greater<std::string>> ans;
  std::string sampInfoStr = getSimpleSampInfo();
  uint32_t pos = 0;
  uint32_t emptyRepAmount = 8;
  for (const auto &clus : iter::reverse(collapsed_.clusters_)) {
    std::stringstream currentInfo;
    currentInfo << sampInfoStr << delim
                << collapsed_.clusters_.size() - 1 - pos;
    currentInfo << delim
                << clus.getStandardInfo(collapsed_.totalReadCount_, delim);
    for (const auto &info : collapsed_.infos_) {
      if (clus.sampInfos_.find(info.first) == clus.sampInfos_.end()) {
        currentInfo << delim << repeatString(delim, emptyRepAmount);
      } else {
        currentInfo << delim << info.first;
        if (excluded_.infos_.find(info.first) == excluded_.infos_.end()) {
          currentInfo << delim << 0 << delim << 0 << delim << 0 << delim << 0;
        } else {
        	currentInfo << delim << excluded_.infos_.at(info.first).readCnt_;
          currentInfo << delim << excluded_.infos_.at(info.first).readCnt_ /
          														input_.infos_.at(info.first).readCnt_;
          currentInfo << delim << excluded_.infos_.at(info.first).chiReadCnt_;
          currentInfo << delim << excluded_.infos_.at(info.first).chiReadCnt_ /
                                      input_.infos_.at(info.first).readCnt_;
        }
        currentInfo << delim << clus.sampInfos_.at(info.first).getReadInfo();
        currentInfo << delim << info.second.readCnt_;
      }
    }
    if (checkingExpected) {
      currentInfo << delim << clus.expectsString;
    }
    ans.insert({clus.getStubName(true), currentInfo.str()});
    ++pos;
  }
  return ans;
}
std::vector<sampleCluster> sampleCollapse::createOutput(
    bool renameFirst, const std::string &sortBy) {
  if (renameFirst) {
    renameClusters(sortBy);
  }
  std::vector<sampleCluster> output;
  for (const auto &out : collapsed_.clusters_) {
    output.emplace_back(sampleCluster(bibseq::cluster(out.createRead())));
  }
  readVec::allUpdateName(output);
  // need to change subreads names so look up can happen correctly
  for (auto &out : output) {
    for (auto &subRead : out.reads_) {
      subRead.seqBase_.name_ = out.seqBase_.name_;
    }
  }
  return output;
}
// writing
void sampleCollapse::writeExcluded(const std::string &outDirectory,
                                   const std::string &outFormat, bool overWrite,
                                   bool exitOnFailureToWrite) const {
  excluded_.writeClusters(outDirectory + sampName_, outFormat, overWrite,
                          exitOnFailureToWrite);
}
void sampleCollapse::writeInitial(const std::string &outDirectory,
                                  const std::string &outFormat, bool overWrite,
                                  bool exitOnFailureToWrite) const {
  input_.writeClusters(outDirectory + sampName_, outFormat, overWrite,
                       exitOnFailureToWrite);
}
void sampleCollapse::writeFinal(const std::string &outDirectory,
                                const std::string &outFormat, bool overWrite,
                                bool exitOnFailureToWrite) const {
  collapsed_.writeClusters(outDirectory + sampName_, outFormat, overWrite,
                           exitOnFailureToWrite);
}
void sampleCollapse::writeFinalOrginals(const std::string &outDirectory) const {
  for (const auto &clus : collapsed_.clusters_) {
    clus.writeOutClusters(outDirectory);
  }
}
}  // napsace collapse
}  // namespace bib
