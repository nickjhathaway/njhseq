
#include "sampInfo.hpp"
#include "bibseq/utils.h"

namespace bibseq {


sampInfo::sampInfo()
    : runName_(""),
      runReadCnt_(0),
			readCnt_(0),
      fraction_(0),
      numberOfClusters_(0),
      chiReadCnt_(0),
      chiNumberOfClusters_(0) {}

sampInfo::sampInfo(const readObject& cr) {
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
  fraction_ = 1;
  runReadCnt_ = readCnt_;
};

sampInfo::sampInfo(const std::string & runName, double totalRunCnt) :
		runName_(runName), runReadCnt_(totalRunCnt), readCnt_(0),  fraction_(0), numberOfClusters_(
				0), chiReadCnt_(0), chiNumberOfClusters_(0) {
}

void sampInfo::resetBasicInfo(){
  readCnt_ = 0;
  fraction_ = 0;
  numberOfClusters_ = 0;
  chiReadCnt_ = 0;
  chiNumberOfClusters_ = 0;
}

// updates
void sampInfo::update(const readObject& cr) {
  if (cr.seqBase_.name_.find("CHI") != std::string::npos) {
    chiReadCnt_ += cr.seqBase_.cnt_;
    ++chiNumberOfClusters_;
  }
  ++numberOfClusters_;
  readCnt_ += cr.seqBase_.cnt_;
}

void sampInfo::updateRunReadCnt(double runReadCnt) {
  runReadCnt_ = runReadCnt;
  updateFraction();
}

void sampInfo::updateFraction() {
	fraction_ = readCnt_ / runReadCnt_;
}

std::string sampInfo::getReadInfo(const std::string& delim) const {
  return estd::to_string(fraction_) + delim + estd::to_string(readCnt_) +
         delim + estd::to_string(numberOfClusters_);
}
std::string sampInfo::getReadInfo(uint32_t cnt, const std::string& delim ) const {
  return estd::to_string(readCnt_ / cnt) + delim + estd::to_string(readCnt_) +
         delim + estd::to_string(numberOfClusters_);
}
std::string sampInfo::getChimeraInfo(const std::string& delim ) const {
  return estd::to_string(chiReadCnt_ / runReadCnt_) + delim +
         estd::to_string(chiReadCnt_) + delim +
         estd::to_string(chiNumberOfClusters_);
}
std::string sampInfo::getChimeraInfo(uint32_t cnt, const std::string& delim ) const {
  return estd::to_string(chiReadCnt_ / cnt) + delim +
         estd::to_string(chiReadCnt_) + delim +
         estd::to_string(chiNumberOfClusters_);
}

}  // namespace bibseq
