#pragma once
//
//  sampInfo.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 05/03/15.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//


#include "bibseq/common/allSystemIncludes.h"
#include "bibseq/objects/seqObjects/readObject.hpp"

namespace bibseq {
class sampInfo {
public:

  sampInfo();
  sampInfo(const readObject& cr);
  sampInfo(const std::string & runName, double totalRunCnt);
  // updates
  void update(const readObject& cr);
  void updateRunReadCnt(double runReadCnt);
  void updateFraction();

  void resetBasicInfo();
  // samp info

  std::string runName_;
  double runReadCnt_; //total number of reads associated with runName_(MID rep) in all clusters

  double readCnt_; //amount of reads associated with runName_ (MID rep)
  double fraction_; //fraction of rep reads in this cluster,
  uint32_t numberOfClusters_;  //total number of clusters
  // chimeric info
  double chiReadCnt_;
  uint32_t chiNumberOfClusters_;

  std::string getReadInfo(const std::string& delim = "\t") const;
  std::string getReadInfo(uint32_t cnt, const std::string& delim = "\t") const;
  std::string getChimeraInfo(const std::string& delim = "\t") const;
  std::string getChimeraInfo(uint32_t cnt, const std::string& delim = "\t") const;
};


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "sampInfo.cpp"
#endif
