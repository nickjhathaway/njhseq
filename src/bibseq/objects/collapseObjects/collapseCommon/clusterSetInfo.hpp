#pragma once
/*
 * clusterSetInfo.hpp
 *
 *  Created on: Aug 7, 2016
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "bibseq/objects/helperObjects/sampInfo.hpp"


namespace bibseq {
namespace collapse {

class clusterSetInfo {
public:
	clusterSetInfo();

	template<typename T>
	clusterSetInfo(const std::vector<T> & reads) {
		setInfo(reads);
	}

	std::map<std::string, sampInfo> infos_;
	double totalReadCount_ = 0;
	uint32_t numberOfClusters_ = 0;
	std::vector<uint32_t> cois_;

	template<typename T>
	void setInfo(const std::vector<T> & reads) {
		clear();
		for (const auto & read : reads) {
			updateInfo(getSeqBase(read));
		}
		//for cluster set info set the runReadCnt equal to the readCnt count
		resetRunReadCnt();
	}

	template<typename T>
	void updateInfo(const std::vector<T> & reads) {
		for (const auto & read : reads) {
			updateInfo(getSeqBase(read));
		}
		//for cluster set info set the runReadCnt equal to the readCnt count
		resetRunReadCnt();
	}

	void updateInfo(const seqInfo & read);

	void clear();

	void resetRunReadCnt();

	void updateCoi(uint32_t coi);

	double meanCoi() const;
	double medianCoi() const;
	uint32_t maxCoi() const;

	uint32_t minCoi() const;

	static std::string getCoiInfoHeader(const std::string & delim);

	std::string getCoiInfo(const std::string & delim) const;

	static VecStr getCoiInfoHeaderVec();
	VecStr getCoiInfoVec() const;

};

} // namespace collapse
} // namespace bibseq
