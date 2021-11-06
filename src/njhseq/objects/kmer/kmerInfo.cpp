#include "kmerInfo.hpp"
#include "njhseq/helpers/seqUtil.hpp"
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


double kmerInfo::DetailedKmerDist::getDistTotalShared() const {
	return static_cast<double>(totalShared_) / (totalKmersIn1_ + totalKmersIn2_);
}

double kmerInfo::DetailedKmerDist::getDistTotalSharedLenAdjusted() const {
	return static_cast<double>(totalShared_)
			/ (std::min(totalKmersIn1_, totalKmersIn2_));
}

double kmerInfo::DetailedKmerDist::getDistUniqueShared() const {
	return static_cast<double>(totalUniqShared_) / totalUniqBetween_;
}

double kmerInfo::DetailedKmerDist::getDistUniqueSharedLenAdjusted() const {
	return static_cast<double>(totalUniqShared_)
			/ (std::min(totalUniqKmersIn1_, totalUniqKmersIn2_));
}




kmerInfo::kmerInfo() :
		kLen_(1), seqLen_(0) {
}

kmerInfo::kmerInfo(const std::string & seq, uint32_t kLength, bool setReverse) :
		kLen_(kLength), seqLen_(seq.size()), infoSet_(true) {
	setKmers(seq, kLength, setReverse);
}

kmerInfo::kmerInfo(const std::string & seq, uint32_t kLength,
			size_t pos, uint32_t len, bool setReverse){
	setKmersFromPortion(seq, kLength, pos, len, setReverse);
}

void kmerInfo::setKmersFromPortion(const std::string & seq, uint32_t kLength,
		size_t pos, uint32_t len, bool setReverse){

	//safety checks
	//check pos
	if (pos > seq.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": pos: " << pos
				<< " is greater than seq size: " << seq.size() << std::endl;
		throw std::runtime_error { ss.str() };
	}
	//check len
	if (pos + len > seq.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": pos: " << pos
				<< " plus len: " << len << " is greater than seq size: " << seq.size() << std::endl;
		throw std::runtime_error { ss.str() };
	}
	//check to see if it will be possible to get a kmer
	if(pos + kLength > seq.size()
			|| kLength > len
			){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": with pos: " << pos
				<< ", kmer size of " << kLength
				<< ",  a len of: " << len
				<< ", and a seq size of " << seq.size()
				<< " it is not possible to get kmers" << std::endl;
		throw std::runtime_error { ss.str() };
	}

	//reset information
	infoSet_ = true;
	kLen_ = kLength;
	seqLen_ = len;
	seqPos_ = pos;
	kmers_.clear();
	kmersRevComp_.clear();
	//set kmer information in the current direction
	/**@todo check if this needs a plus 1 */
	for (const auto seqPos : iter::range(pos, pos + len - kLen_)) {
		auto currentK = seq.substr(seqPos, kLen_);
		auto k = kmers_.find(currentK);
		if (k != kmers_.end()) {
			k->second.addPosition(seqPos);
		} else {
			kmers_.emplace(currentK, kmer(currentK, seqPos));
		}
	}

	//if needed set kmer information for the reverse direction
	if (setReverse) {
		std::string reverseComplement = seqUtil::reverseComplement(seq.substr(pos, len), "DNA");
		for (const auto seqPos : iter::range(reverseComplement.size() + 1 - kLen_)) {
			auto currentK = reverseComplement.substr(seqPos, kLen_);
			auto k = kmersRevComp_.find(currentK);
			if (k != kmersRevComp_.end()) {
				k->second.addPosition((seq.size() - pos - len) + seqPos);
			} else {
				kmersRevComp_.emplace(currentK, kmer(currentK, (seq.size() - pos - len) + seqPos));
			}
		}
	}
}


void kmerInfo::updateKmers(const std::string & seq, bool setReverse) {

	if(kLen_ <=seq.length()){
		infoSet_ = true;
		seqLen_ += seq.length();
		//set kmer information in the current direction
		for (const auto pos : iter::range(seq.size() + 1 - kLen_)) {
			auto currentK = seq.substr(pos, kLen_);
			auto k = kmers_.find(currentK);
			if (k != kmers_.end()) {
				k->second.addPosition(pos);
			} else {
				kmers_.emplace(currentK, kmer(currentK, pos));
			}
		}
		//if needed set kmer information for the reverse direction
		if (setReverse) {
			std::string reverseComplement = seqUtil::reverseComplement(seq, "DNA");
			for (const auto pos : iter::range(reverseComplement.size() + 1 - kLen_)) {
				auto currentK = reverseComplement.substr(pos, kLen_);
				auto k = kmersRevComp_.find(currentK);
				if (k != kmersRevComp_.end()) {
					k->second.addPosition(pos);
				} else {
					kmersRevComp_.emplace(currentK, kmer(currentK, pos));
				}
			}
		}
	}
}

void kmerInfo::setKmers(const std::string & seq, uint32_t kLength, bool setReverse) {
	//reset information
	infoSet_ = true;
	kLen_ = kLength;
	seqLen_ = seq.size();
	kmers_.clear();
	kmersRevComp_.clear();
	if(kLength <=seq.length()){
		//set kmer information in the current direction
		for (const auto pos : iter::range(seq.size() + 1 - kLen_)) {
			auto currentK = seq.substr(pos, kLen_);
			auto k = kmers_.find(currentK);
			if (k != kmers_.end()) {
				k->second.addPosition(pos);
			} else {
				kmers_.emplace(currentK, kmer(currentK, pos));
			}
		}
		//if needed set kmer information for the reverse direction
		if (setReverse) {
			std::string reverseComplement = seqUtil::reverseComplement(seq, "DNA");
			for (const auto pos : iter::range(reverseComplement.size() + 1 - kLen_)) {
				auto currentK = reverseComplement.substr(pos, kLen_);
				auto k = kmersRevComp_.find(currentK);
				if (k != kmersRevComp_.end()) {
					k->second.addPosition(pos);
				} else {
					kmersRevComp_.emplace(currentK, kmer(currentK, pos));
				}
			}
		}
	}
}

std::pair<uint32_t, double> kmerInfo::compareKmers(
		const kmerInfo & info) const {
	uint32_t kShared = 0;
	for (const auto & k : kmers_) {
		auto otherK = info.kmers_.find(k.first);
		if (otherK != info.kmers_.end()) {
			kShared += std::min(otherK->second.count_, k.second.count_);
		}
	}
	return {kShared,
		kShared/static_cast<double>(std::min(info.seqLen_, seqLen_) + 1 - kLen_)};
}

kmerInfo::DetailedKmerDist kmerInfo::compareKmersDetailed(const kmerInfo & info) const{
	kmerInfo::DetailedKmerDist ret;
	if(info.kmers_.size() > kmers_.size()){
		for (const auto & k : kmers_) {
			auto otherK = info.kmers_.find(k.first);
			if (otherK != info.kmers_.end()) {
				++ret.totalUniqShared_;
				ret.totalShared_ += std::min(otherK->second.count_, k.second.count_);
			}
		}
	}else{
		for (const auto & k : info.kmers_) {
			auto otherK = kmers_.find(k.first);
			if (otherK != kmers_.end()) {
				++ret.totalUniqShared_;
				ret.totalShared_ += std::min(otherK->second.count_, k.second.count_);
			}
		}
	}

	{
		std::unordered_set<std::string> uniqKmersBetween = njh::vecToUOSet(getVectorOfMapKeys(kmers_));
		njh::addVecToUOSet(getVectorOfMapKeys(info.kmers_), uniqKmersBetween);
		ret.totalUniqBetween_ = uniqKmersBetween.size();
	}

	ret.totalKmersIn1_ = seqLen_ + 1 - kLen_;
	ret.totalKmersIn2_ = info.seqLen_ + 1 - kLen_;

	ret.totalUniqKmersIn1_ = kmers_.size();
	ret.totalUniqKmersIn2_ = info.kmers_.size();

	return ret;
}



std::pair<uint32_t, double> kmerInfo::compareKmersRevComp(
		const kmerInfo & info) const {
	uint32_t kShared = 0;
	for (const auto & k : kmers_) {
		auto otherK = info.kmersRevComp_.find(k.first);
		if (otherK != info.kmersRevComp_.end()) {
			kShared += std::min(otherK->second.count_, k.second.count_);
		}
	}
	return {kShared,
		kShared/static_cast<double>(std::min(info.seqLen_,
						seqLen_) + 1 - kLen_)};
}

kmerInfo::DetailedKmerDist kmerInfo::compareKmersRevCompDetailed(const kmerInfo & info) const{
	kmerInfo::DetailedKmerDist ret;
	if(info.kmersRevComp_.size() > kmers_.size()){
		for (const auto & k : kmers_) {
			auto otherK = info.kmersRevComp_.find(k.first);
			if (otherK != info.kmersRevComp_.end()) {
				++ret.totalUniqShared_;
				ret.totalShared_ += std::min(otherK->second.count_, k.second.count_);
			}
		}
	}else{
		for (const auto & k : info.kmersRevComp_) {
			auto otherK = kmers_.find(k.first);
			if (otherK != kmers_.end()) {
				++ret.totalUniqShared_;
				ret.totalShared_ += std::min(otherK->second.count_, k.second.count_);
			}
		}
	}

	{
		std::unordered_set<std::string> uniqKmersBetween = njh::vecToUOSet(getVectorOfMapKeys(kmers_));
		njh::addVecToUOSet(getVectorOfMapKeys(info.kmersRevComp_), uniqKmersBetween);
		ret.totalUniqBetween_ = uniqKmersBetween.size();
	}

	ret.totalKmersIn1_ = seqLen_ + 1 - kLen_;
	ret.totalKmersIn2_ = info.seqLen_ + 1 - kLen_;

	ret.totalUniqKmersIn1_ = kmers_.size();
	ret.totalUniqKmersIn2_ = info.kmersRevComp_.size();

	return ret;
}



std::pair<uint32_t, double> kmerInfo::compareKmers(const kmerInfo & info,
		uint32_t startPos, uint32_t windowSize) const {
	/**@todo Consider putting in a check for windowSize vs kmersize and for startPos though it would be better to just check once outside of function*/
	uint32_t kmersShared = 0;
	//calculate the maximum number of kmers possibly shared in a window this size
	double maxKMers = windowSize - kLen_ + 1;
	//get kmers at this position
	for (const auto & k : kmers_) {
		uint32_t kmersInWindow = 0;
		for (const auto & pos : k.second.positions_) {
			if (pos >= startPos && pos <= (startPos + windowSize - kLen_)) {
				++kmersInWindow;
			}
			if (pos > (startPos + windowSize - kLen_)) {
				break;
			}
		}
		//if this kmer isn't found at this position go on
		if(kmersInWindow == 0){
			continue;
		}
		uint32_t kmersInWindowOther = 0;
		auto otherK = info.kmers_.find(k.first);
		if (otherK != info.kmers_.end()) {
			for (const auto & pos : otherK->second.positions_) {
				if (pos >= startPos && pos <= (startPos + windowSize - kLen_)) {
					++kmersInWindowOther;
				}
				if (pos > (startPos + windowSize - kLen_)) {
					break;
				}
			}
		}
		kmersShared += std::min(kmersInWindow, kmersInWindowOther);
	}
	return {kmersShared, kmersShared/maxKMers};
}

std::pair<uint32_t, double> kmerInfo::compareKmers(const kmerInfo & info,
		uint32_t startPos,uint32_t otherStartPos, uint32_t windowSize) const{
	/**@todo Consider putting in a check for windowSize vs kmersize and for startPos though it would be better to just check once outside of function*/
	uint32_t kmersShared = 0;
	//calculate the maximum number of kmers possibly shared in a window this size
	double maxKMers = windowSize - kLen_ + 1;
	//get kmers at this position
	for (const auto & k : kmers_) {
		uint32_t kmersInWindow = 0;
		for (const auto & pos : k.second.positions_) {
			if (pos >= startPos && pos <= (startPos + windowSize - kLen_)) {
				++kmersInWindow;
			}
			if (pos > (startPos + windowSize - kLen_)) {
				break;
			}
		}
		//if this kmer isn't found at this position go on
		if(kmersInWindow == 0){
			continue;
		}
		uint32_t kmersInWindowOther = 0;
		auto otherK = info.kmers_.find(k.first);
		if (otherK != info.kmers_.end()) {
			for (const auto & pos : otherK->second.positions_) {
				if (pos >= otherStartPos && pos <= (otherStartPos + windowSize - kLen_)) {
					++kmersInWindowOther;
				}
				if (pos > (otherStartPos + windowSize - kLen_)) {
					break;
				}
			}
		}
		kmersShared += std::min(kmersInWindow, kmersInWindowOther);
	}
	return {kmersShared, kmersShared/maxKMers};
}

std::pair<uint32_t, double> kmerInfo::compareSubKmersToFull(
		const kmerInfo & info, uint32_t startPos, uint32_t windowSize) const {
	uint32_t kmersShared = 0;
	double maxKMers = windowSize - kLen_ + 1;
	for (const auto & k : kmers_) {
		uint32_t kmersInWindow = 0;
		for (const auto & pos : k.second.positions_) {
			if (pos >= startPos && pos <= (startPos + windowSize - kLen_)) {
				++kmersInWindow;
			}
			if (pos > (startPos + windowSize - kLen_)) {
				break;
			}
		}
		uint32_t kmersInWindowOther = 0;
		auto otherK = info.kmers_.find(k.first);
		if (otherK != info.kmers_.end()) {
			kmersInWindowOther = otherK->second.count_;
		}
		kmersShared += std::min(kmersInWindow, kmersInWindowOther);
	}
	return {kmersShared, kmersShared/maxKMers};
}

std::unordered_map<size_t, std::pair<uint32_t, double>> kmerInfo::slideCompareKmers(
		const kmerInfo & info, uint32_t windowSize, uint32_t windowStepSize) const {
	std::unordered_map<size_t, std::pair<uint32_t, double>> ret;
	uint64_t minLen = std::min(seqLen_, info.seqLen_);
	for (const auto pos : iter::range<uint32_t>(0, minLen - windowSize + 1,
			windowStepSize)) {
		ret.emplace(pos,compareKmers(info, pos, windowSize));
	}
	return ret;
}

std::unordered_map<size_t, std::pair<uint32_t, double>> kmerInfo::slideCompareSubKmersToFull(
		const kmerInfo & info, uint32_t windowSize, uint32_t windowStepSize) const {
	std::unordered_map<size_t, std::pair<uint32_t, double>> ret;
	for (const auto pos : iter::range<uint32_t>(0, seqLen_ - windowSize + 1,
			windowStepSize)) {
		ret.emplace(pos,compareSubKmersToFull(info, pos, windowSize));
	}
	return ret;
}

std::unordered_map<size_t,std::unordered_map<size_t,std::pair<uint32_t, double>>> kmerInfo::slideCompareSubKmersToSubKmers(
		const kmerInfo & info, uint32_t windowSize,
		uint32_t windowStepSize) const{
	std::unordered_map<size_t,std::unordered_map<size_t,std::pair<uint32_t, double>>> ret;
	for (const auto pos : iter::range<uint32_t>(0, seqLen_ - windowSize + 1,
			windowStepSize)) {
		for (const auto otherPos : iter::range<uint32_t>(0, info.seqLen_ - windowSize + 1,
				windowStepSize)) {
			ret[pos][otherPos] = compareKmers(info, pos, otherPos, windowSize);
		}
	}
	return ret;
}

uint32_t kmerInfo::getMinimumNonRedundant(const std::string & seq){
	if(seq.size() <= 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error seq size should at least be " << 3 << "\n";
		throw std::runtime_error{ss.str()};
	}
  uint32_t klen = 2;
  bool foundLength = false;
  while(klen < seq.size()  && !foundLength){
		kmerInfo kinfo(seq, klen, false);
		foundLength = true;
		for(const auto & k : kinfo.kmers_){
			if(k.second.count_ > 1){
				foundLength = false;
				break;
			}
		}
		if(!foundLength){
			++klen;
		}
  }
  return klen;
}


double kmerInfo::computeKmerEntropy() const{
	uint32_t totalCount = 0;
	uint32_t totalKmers = std::pow(4, kLen_);
	double div = std::sqrt(totalKmers);
	for(const auto & k : kmers_){
		totalCount += k.second.count_;
	}
	double sum = 0;
	for(const auto & k : kmers_){
		if(k.second.count_ != totalCount){
			double frac = k.second.count_/static_cast<double>(totalCount);
			sum += frac * std::log(frac)/std::log(div);
		}
	}
	return -1 * sum;
}

}  // namespace njh
