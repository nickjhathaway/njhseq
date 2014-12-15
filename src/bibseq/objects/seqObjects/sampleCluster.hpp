#pragma once
//
//  sampleCluster.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/10/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/cluster.hpp"

namespace bibseq {

struct sampInfo {
  sampInfo()
      : runName_(""),
        readCnt_(0),
        runReadCnt_(0),
        fraction_(0),
        numberOfClusters_(0),
        chiReadCnt_(0),
        chiNumberOfClusters_(0) {}
  sampInfo(const readObject& cr) {
    if (cr.seqBase_.name_.find("CHI") == std::string::npos) {
      chiReadCnt_ = 0;
      chiNumberOfClusters_ = 0;
    } else {
      chiReadCnt_ = cr.seqBase_.cnt_;
      chiNumberOfClusters_ = 1;
    }
    numberOfClusters_ = 1;
    readCnt_ = cr.seqBase_.cnt_;
    runName_ = cr.getOwnSampName();
  };
  // updates
  void update(const readObject& cr) {
    if (cr.seqBase_.name_.find("CHI") != std::string::npos) {
      chiReadCnt_ += cr.seqBase_.cnt_;
      ++chiNumberOfClusters_;
    }
    ++numberOfClusters_;
    readCnt_ += cr.seqBase_.cnt_;
  }
  void updateRunReadCnt(double runReadCnt) {
    runReadCnt_ = runReadCnt;
    updateFraction();
  }
  void updateFraction() { fraction_ = readCnt_ / runReadCnt_; }
  // samp info
  std::string runName_;
  double readCnt_;
  double runReadCnt_;
  double fraction_;
  int numberOfClusters_;
  // chimeric info
  double chiReadCnt_;
  int chiNumberOfClusters_;

  std::string getReadInfo(const std::string& delim = "\t") const {
    return std::to_string(fraction_) + delim + std::to_string(readCnt_) +
           delim + std::to_string(numberOfClusters_);
  }
  std::string getReadInfo(int cnt, const std::string& delim = "\t") const {
    return std::to_string(readCnt_ / cnt) + delim + std::to_string(readCnt_) +
           delim + std::to_string(numberOfClusters_);
  }
  std::string getChimeraInfo(const std::string& delim = "\t") const {
    return std::to_string(chiReadCnt_ / runReadCnt_) + delim +
           std::to_string(chiReadCnt_) + delim +
           std::to_string(chiNumberOfClusters_);
  }
  std::string getChimeraInfo(int cnt, const std::string& delim = "\t") const {
    return std::to_string(chiReadCnt_ / cnt) + delim +
           std::to_string(chiReadCnt_) + delim +
           std::to_string(chiNumberOfClusters_);
  }
};

class sampleCluster : public cluster {

 public:
  template <typename T>
  sampleCluster(const T& firstCluster) {
    // std::cout<<"mark con1"<<std::endl;
    seqBase_.name_ = firstCluster.seqBase_.name_;
    seqBase_.seq_ = firstCluster.seqBase_.seq_;
    seqBase_.qual_ = firstCluster.seqBase_.qual_;
    firstReadName = firstCluster.seqBase_.name_;
    firstReadCount = firstCluster.seqBase_.cnt_;
    // normalizedFraction = 0;
    seqBase_.frac_ = firstCluster.seqBase_.frac_;
    cumulativeFraction = firstCluster.seqBase_.frac_;
    seqBase_.cnt_ = firstCluster.seqBase_.cnt_;
    normalizedFraction = firstCluster.normalizedFraction;
    sampName = getOwnSampName();
    // create tempObject to add to reads
    readObject tempObject(seqInfo(firstCluster.seqBase_.name_,
                                  firstCluster.seqBase_.seq_,
                                  firstCluster.seqBase_.qual_));
    tempObject.seqBase_.frac_ = firstCluster.seqBase_.frac_;
    tempObject.cumulativeFraction = firstCluster.cumulativeFraction;
    tempObject.seqBase_.cnt_ = firstCluster.seqBase_.cnt_;
    reads_.emplace_back(tempObject);
    updateInfoWithRead(tempObject, 0);
    // std::cout<<"sampName: "<<sampName<<std::endl;
    // reads.clear();
    // addRead(firstCluster);
    remove = false;
    // normalizedFraction=firstCluster.normalizedFraction;
    // std::cout<<"mark con2"<<std::endl;
  }
  sampleCluster() {}

  virtual void updateName();
  virtual void setName(const std::string& newName);
  void addRead(const cluster& cr);

  void update(const std::map<std::string, sampInfo>& infos);
  void update(const cluster& cr);
  void updateInfoWithRead(const readObject& read, uint32_t pos);
  void updateFractionInfo();
  // info
  std::string getChimeraInfo(const std::string& delim = "\t") const;
  std::string getStandardInfo(int sampReadCnt,
                              const std::string& delim = "\t") const;
  std::string getPopStandardInfo(double popReadCnt, uint32_t popClusNum,
                                 uint32_t sampNum, bool addProtein,
                                 const std::string& delim = "\t") const;
  // members
  uint32_t numberOfRuns_;
  // std::map<std::string, std::vector<readObject>> sampleClusters_;
  std::map<std::string, std::vector<uint32_t>> sampleClusters_;
  std::map<std::string, sampInfo> sampInfos_;

  template <typename T>
  static void updateAllClusters(std::vector<T>& clusters,
                                const std::map<std::string, sampInfo>& infos) {
    for (auto& clus : clusters) {
      clus.update(infos);
    }
  }
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "sampleCluster.cpp"
#endif
