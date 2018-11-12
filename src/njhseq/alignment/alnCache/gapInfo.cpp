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
/*
 * gapInfo.cpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */


#include "gapInfo.hpp"

namespace njhseq {


gapInfo::gapInfo() :
		pos_(0), size_(0), gapInA_(false) {
}
gapInfo::gapInfo(uint32_t pos, uint32_t size, bool gapInA) :
		pos_(pos), size_(size), gapInA_(gapInA) {
}


void gapInfo::writeInfo(std::ostream& out) const {
	out << pos_ << "," << size_ << "," << gapInA_ << std::endl;
}
void gapInfo::writeInfoNoEnd(std::ostream& out) const {
	out << pos_ << "," << size_ << "," << gapInA_;
}
bool gapInfo::operator>(const gapInfo& otherInfo) const {
	return pos_ > otherInfo.pos_;
}
bool gapInfo::operator<(const gapInfo& otherInfo) const {
	return pos_ < otherInfo.pos_;
}
bool gapInfo::operator==(const gapInfo& otherInfo) const {
	return (pos_ == otherInfo.pos_ && size_ == otherInfo.size_
			&& gapInA_ == otherInfo.gapInA_);
}
bool gapInfo::operator<=(const gapInfo& otherInfo) const {
	return pos_ <= otherInfo.pos_;
}
bool gapInfo::operator>=(const gapInfo& otherInfo) const {
	return pos_ >= otherInfo.pos_;
}

Json::Value gapInfo::toJson() const {
	Json::Value ret;
	ret["class"] = "njhseq::gapInfo";
	ret["pos_"] = njh::json::toJson(pos_);
	ret["size_"] = njh::json::toJson(size_);
	ret["gapInA_"] = njh::json::toJson(gapInA_);
	return ret;
}

}


