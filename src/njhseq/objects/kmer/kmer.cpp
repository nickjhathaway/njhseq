#include "kmer.hpp"
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

// constructors
kmer::kmer() :
		k_(""), count_(0), readCnt_(0) {
}

kmer::kmer(const std::string & kSeq) :
		k_(kSeq), count_(0), readCnt_(0) {
}

kmer::kmer(const std::string& kSeq, uint32_t firstPos) :
		k_(kSeq), count_(1), positions_( { firstPos }), readCnt_(1) {
}
kmer::kmer(const std::string& kSeq, uint32_t firstPos,
		const std::string& firstName, uint32_t numReads) :
		k_(kSeq), count_(1), positions_( { firstPos }), readCnt_(numReads) {
	++names_[firstName];
}

Json::Value kmer::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson("njhseq::kmer");
	ret["k_"] = njh::json::toJson(k_);
	ret["count_"] = njh::json::toJson(count_);
	ret["names_"] = njh::json::toJson(names_);
	ret["positions_"] = njh::json::toJson(positions_);
	ret["readCnt_"] = njh::json::toJson(readCnt_);
	return ret;
}

void kmer::addPosition(uint32_t pos) {
	positions_.push_back(pos);
	++count_;
}

void kmer::addPosition(uint32_t pos, const std::string& name,
		uint32_t numReads) {
	positions_.push_back(pos);
	count_ += numReads;
	readCnt_ += numReads;
	++names_[name];
}

void kmer::increaseCnt(uint32_t readCount){
	count_ += readCount;
	readCnt_ += readCount;
}




// comparisons
bool kmer::operator<(const kmer& otherKmer) const {
	if (positions_[0] < otherKmer.positions_[0]) {
		return true;
	} else {
		return false;
	}
}

bool kmer::operator>(const kmer& otherKmer) const {
	if (positions_[0] > otherKmer.positions_[0]) {
		return true;
	} else {
		return false;
	}
}

bool kmer::operator==(const kmer& otherKmer) const {
	if (positions_[0] == otherKmer.positions_[0]) {
		return true;
	} else {
		return false;
	}
}
bool kmer::operator>=(const kmer& otherKmer) const {
	if (positions_[0] >= otherKmer.positions_[0]) {
		return true;
	} else {
		return false;
	}
}

bool kmer::operator<=(const kmer& otherKmer) const {
	if (positions_[0] <= otherKmer.positions_[0]) {
		return true;
	} else {
		return false;
	}
}

bool kmer::operator==(const std::string& kString) const {
	return (k_ == kString);
}



}  // namespace njh
