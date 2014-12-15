//
//  collapseCommon.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/7/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "collapseCommon.hpp"

namespace bibseq {
namespace collapse {
void clusterSet::updateSubClustersPositions(bool clearFirst) {
  if (clearFirst) {
    subClustersPositions_.clear();
  }
  uint32_t pos = 0;
  for (const auto& clus : clusters_) {
    for (const auto& read : clus.reads_) {
      subClustersPositions_[read.getStubName(true)] = pos;
    }
    ++pos;
  }
}
void clusterSet::updateSetInfo(bool clearCurrentInfos) {
  if (clearCurrentInfos) {
    infos_.clear();
  }
  for (const auto& read : clusters_) {
    for (const auto& subRead : read.reads_) {
      updateInfo(subRead, infos_);
    }
    // updateInfo(read, infos_);
  }
  updateCounts();
}
void clusterSet::updateCounts() {
  totalReadCount_ = 0;
  numberOfClusters_ = 0;
  for (const auto& info : infos_) {
    totalReadCount_ += info.second.readCnt_;
    numberOfClusters_ += info.second.numberOfClusters_;
  }
}
// writing
void clusterSet::writeClusters(std::string filename, std::string format,
                               bool overWrite, bool exitOnFailure) const {
  readObjectIO::write(clusters_, filename, format, overWrite, exitOnFailure);
}

}  // collapse
}  // bib