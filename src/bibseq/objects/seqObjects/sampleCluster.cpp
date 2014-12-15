//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include "sampleCluster.hpp"

namespace bibseq {

void sampleCluster::addRead(const cluster& cr) {
  // add the cluster's reads and update the info
  for (const auto& read : cr.reads_) {
    updateInfoWithRead(read, len(reads_));
    reads_.push_back(read);
  }
  // update the fraction and totalCounts
  cumulativeFraction += cr.cumulativeFraction;
  seqBase_.cnt_ += cr.seqBase_.cnt_;
  seqBase_.frac_ = cumulativeFraction / sampInfos_.size();
  normalizedFraction += cr.normalizedFraction;
  needToCalculateConsensus = true;
  updateName();
}

void sampleCluster::update(const std::map<std::string, sampInfo>& infos) {
  // clear and update it all
  sampleClusters_.clear();
  sampInfos_.clear();
  uint32_t pos = 0;
  for (const auto& read : reads_) {
    updateInfoWithRead(read, pos);
    ++pos;
  }
  for (auto& info : sampInfos_) {
    info.second.updateRunReadCnt(infos.at(info.first).readCnt_);
  }
  updateFractionInfo();
}
void sampleCluster::updateInfoWithRead(const readObject& read, uint32_t pos) {
  // update infos with the read
  if (sampleClusters_.find(read.sampName) == sampleClusters_.end()) {
    sampleClusters_.emplace(read.sampName, std::vector<uint32_t>(1, pos));
    sampInfos_.emplace(read.sampName, sampInfo(read));
  } else {
    sampleClusters_[read.sampName].push_back(pos);
    sampInfos_[read.sampName].update(read);
  }
}
void sampleCluster::updateFractionInfo() {
  // clear the counts and update them
  seqBase_.frac_ = 0;
  cumulativeFraction = 0;
  seqBase_.cnt_ = 0;
  for (const auto& samp : sampInfos_) {
    seqBase_.cnt_ += samp.second.readCnt_;
    cumulativeFraction += samp.second.fraction_;
  }
  /*
  for (const auto & sample : sampleClusters_) {
    for (const auto clus : sample.second) {
      totalCount += reads[clus].totalCount;
      cumulativeFraction += reads[clus].fraction;
    }
  }*/
  seqBase_.frac_ = cumulativeFraction / sampleClusters_.size();
  numberOfRuns_ = len(sampleClusters_);
  updateName();
}

void sampleCluster::updateName() {
  seqBase_.name_ = getStubName(false) + "_f" + std::to_string(seqBase_.frac_);
}
void sampleCluster::setName(const std::string& newName) {
  seqBase_.name_ = newName + "_f" + std::to_string(seqBase_.frac_);
}
std::string sampleCluster::getChimeraInfo(const std::string& delim) const {
  int chiRepCnt = 0;
  int chiReadCnt = 0;
  int chiClusCnt = 0;
  for (const auto& info : sampInfos_) {
    if (info.second.chiReadCnt_ > 0) {
      ++chiRepCnt;
      chiReadCnt += info.second.chiReadCnt_;
      chiClusCnt += info.second.chiNumberOfClusters_;
    }
  }
  return std::to_string(chiReadCnt) + delim + std::to_string(chiClusCnt) +
         delim + std::to_string(chiRepCnt);
}
std::string sampleCluster::getStandardInfo(int sampReadCnt,
                                           const std::string& delim) const {
  std::stringstream outStream;
  outStream << seqBase_.frac_ << delim << seqBase_.cnt_ / sampReadCnt << delim
            << sampInfos_.size() << delim << seqBase_.seq_ << delim
            << reads_.size() << delim << getChimeraInfo(delim) << delim
            << vectorToString(readVec::getNames(reads_), ",");
  return outStream.str();
}
std::string sampleCluster::getPopStandardInfo(double popReadCnt,
                                              uint32_t popClusNum,
                                              uint32_t sampNum, bool addProtein,
                                              const std::string& delim) const {
  std::stringstream outStream;
  outStream << getStubName(false) << delim << popReadCnt << delim << popClusNum
            << delim << cumulativeFraction / sampNum << delim
            << cumulativeFraction << delim << seqBase_.frac_ << delim
            << seqBase_.cnt_ / popReadCnt << delim << sampInfos_.size() << delim
            << (double)sampInfos_.size() / sampNum << delim << seqBase_.cnt_
            << delim << reads_.size() << delim
            << vectorToString(readVec::getNames(reads_), ",") << delim
            << seqBase_.seq_;
  if (addProtein) {
    outStream << delim << getProteinFromcDNA(false);
  }
  return outStream.str();
}
}  // namespace bib
