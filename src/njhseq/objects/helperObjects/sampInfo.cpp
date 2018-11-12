#include "sampInfo.hpp"
#include "njhseq/utils.h"
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
namespace njhseq {

sampInfo::sampInfo() :
		runName_(""), runReadCnt_(0), readCnt_(0), fraction_(0), numberOfClusters_(
				0), chiReadCnt_(0), chiNumberOfClusters_(0) {
}

sampInfo::sampInfo(const seqInfo& cr) :
		runName_(cr.getOwnSampName()), runReadCnt_(cr.cnt_), readCnt_(0), fraction_(
				0), numberOfClusters_(0), chiReadCnt_(0), chiNumberOfClusters_(0) {
	update(cr);
	updateFraction();
}

sampInfo::sampInfo(const std::string & runName, double totalRunCnt) :
		runName_(runName), runReadCnt_(totalRunCnt), readCnt_(0), fraction_(0), numberOfClusters_(
				0), chiReadCnt_(0), chiNumberOfClusters_(0) {
}

void sampInfo::resetBasicInfo() {
	readCnt_ = 0;
	fraction_ = 0;
	numberOfClusters_ = 0;
	chiReadCnt_ = 0;
	chiNumberOfClusters_ = 0;
}

// updates
void sampInfo::update(const seqInfo& cr) {
	if (cr.isChimeric()) {
		chiReadCnt_ += cr.cnt_;
		++chiNumberOfClusters_;
	}
	++numberOfClusters_;
	readCnt_ += cr.cnt_;
}

void sampInfo::updateRunReadCnt(double runReadCnt) {
	runReadCnt_ = runReadCnt;
	updateFraction();
}

void sampInfo::updateFraction() {
	fraction_ = readCnt_ / runReadCnt_;
}

std::string sampInfo::getReadInfo(const std::string& delim) const {
	return estd::to_string(fraction_) + delim + estd::to_string(readCnt_) + delim
			+ estd::to_string(numberOfClusters_);
}
std::string sampInfo::getReadInfo(uint32_t cnt,
		const std::string& delim) const {
	return estd::to_string(readCnt_ / cnt) + delim + estd::to_string(readCnt_)
			+ delim + estd::to_string(numberOfClusters_);
}
std::string sampInfo::getChimeraInfo(const std::string& delim) const {
	return estd::to_string(chiReadCnt_ / runReadCnt_) + delim
			+ estd::to_string(chiReadCnt_) + delim
			+ estd::to_string(chiNumberOfClusters_);
}
std::string sampInfo::getChimeraInfo(uint32_t cnt,
		const std::string& delim) const {
	return estd::to_string(chiReadCnt_ / cnt) + delim
			+ estd::to_string(chiReadCnt_) + delim
			+ estd::to_string(chiNumberOfClusters_);
}

}  // namespace njhseq
