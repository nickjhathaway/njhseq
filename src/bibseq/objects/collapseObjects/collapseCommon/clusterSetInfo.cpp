/*
 * clusterSetInfo.cpp
 *
 *  Created on: Aug 7, 2016
 *      Author: nick
 */


#include "clusterSetInfo.hpp"

namespace bibseq {
namespace collapse {


clusterSetInfo::clusterSetInfo() {
}

void clusterSetInfo::updateInfo(const seqInfo & read) {
	auto search = infos_.find(read.getOwnSampName());
	if (search == infos_.end()) {
		infos_.emplace(read.getOwnSampName(), sampInfo(read));
	} else {
		search->second.update(read);
	}
	totalReadCount_ += read.cnt_;
	++numberOfClusters_;
}

void clusterSetInfo::clear() {
	std::map<std::string, sampInfo>().swap(infos_);
	totalReadCount_ = 0;
	numberOfClusters_ = 0;
}

void clusterSetInfo::resetRunReadCnt() {
	for (auto & i : infos_) {
		i.second.runReadCnt_ = i.second.readCnt_;
	}
}

void clusterSetInfo::updateCoi(uint32_t coi) {
	cois_.emplace_back(coi);
}

double clusterSetInfo::meanCoi() const {
	return vectorMean(cois_);
}
double clusterSetInfo::medianCoi() const {
	return vectorMedianCopy(cois_);
}
uint32_t clusterSetInfo::maxCoi() const {
	return vectorMaximum(cois_);
}

uint32_t clusterSetInfo::minCoi() const {
	return vectorMinimum(cois_);
}

std::string clusterSetInfo::getCoiInfoHeader(const std::string & delim) {
	return vectorToString(getCoiInfoHeaderVec(), delim);
}

std::string clusterSetInfo::getCoiInfo(const std::string & delim) const {
	return vectorToString(getCoiInfoVec(), delim);
}

VecStr clusterSetInfo::getCoiInfoHeaderVec() {
	return VecStr { "p_meanCoi", "p_medianCoi", "p_minCoi", "p_maxCoi" };
}

VecStr clusterSetInfo::getCoiInfoVec() const {
	return toVecStr(meanCoi(), medianCoi(), minCoi(), maxCoi());
}

} // namespace collapse
} // namespace bibseq
