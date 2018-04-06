/*
 * kmerMap.cpp
 *
 *  Created on: Dec 2, 2015
 *      Author: nick
 */

//
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

#include "kmerMap.hpp"

namespace bibseq {

KmerMapBase::KmerMapBase(size_t kLength, uint32_t runCutOff) :
		kLength_(kLength), runCutOff_(runCutOff) {
}
KmerMapBase::KmerMapBase(size_t kLength) :
		kLength_(kLength), runCutOff_(0) {
}
KmerMapBase::KmerMapBase() :
		kLength_(1), runCutOff_(0) {
}

Json::Value KmerMapBase::toJson()const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::KmerMapBase");
	ret["k_"] = bib::json::toJson(kLength_);
	ret["runCutOff_"] = bib::json::toJson(runCutOff_);
	return ret;
}


KmerMap::KmerMap(size_t kLength, uint32_t runCutOff,
		const std::unordered_map<std::string, kmer> & kmers) :
			KmerMapBase(kLength, runCutOff), kmers_(kmers) {
}
KmerMap::KmerMap(size_t kLength, uint32_t runCutOff) : KmerMapBase(kLength, runCutOff){}
KmerMap::KmerMap(size_t kLength) : KmerMapBase(kLength){}
KmerMap::KmerMap() : KmerMapBase(){}


void KmerMap::increaseCounts(const std::string & str) {
	increaseCounts(str, 1);
}

void KmerMap::increaseCounts(const std::string & str, uint32_t readCount) {
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto currentKmer = str.substr(cursor, kLength_);
		auto search = kmers_.find(currentKmer);
		if (search == kmers_.end()) {
			auto currentKmerObj = kmer(currentKmer);
			currentKmerObj.increaseCnt(readCount);
			kmers_.emplace(currentKmer, currentKmerObj);
		} else {
			search->second.increaseCnt(readCount);
		}
		++cursor;
	}
}

void KmerMap::increaseCountsExtra(const std::string & str) {
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto currentKmer = str.substr(cursor, kLength_);
		auto search = kmers_.find(currentKmer);
		if (search == kmers_.end()) {
			kmers_.emplace(currentKmer, kmer(currentKmer, cursor));
		} else {
			search->second.addPosition(cursor);
		}
		++cursor;
	}
}

void KmerMap::increaseCountsExtra(const std::string & str, const std::string & uid,
		uint32_t readCount) {
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto currentKmer = str.substr(cursor, kLength_);
		auto search = kmers_.find(currentKmer);
		if (search == kmers_.end()) {
			kmers_.emplace(currentKmer, kmer(currentKmer, cursor, uid, readCount));
		} else {
			search->second.addPosition(cursor, uid, readCount);
		}
		++cursor;
	}
}


bool KmerMap::isKmerLowFreq(const kmerWithPos & k) const {
	return isKmerLowFreq(k, runCutOff_);
}
bool KmerMap::isKmerLowFreq(const kmerWithPos & k, uint32_t runCutOff) const {
	auto search = kmers_.find(k.kmer_);
	return (search == kmers_.end() || search->second.readCnt_ <= runCutOff);
}

uint32_t KmerMap::getKmerFreq(const kmerWithPos & k) const {
	auto search = kmers_.find(k.kmer_);
	if (search == kmers_.end()) {
		return 0;
	} else {
		return search->second.readCnt_;
	}
}

Json::Value KmerMap::toJson()const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::KmerMap");
	ret["super"] = KmerMapBase::toJson();
	ret["kmers_"] = bib::json::toJson(kmers_);
	return ret;
}

KmerMap::~KmerMap() {
}



KmerMapByPos::KmerMapByPos(size_t kLength, uint32_t runCutOff,
		const std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & kmers) :
			KmerMapBase(kLength, runCutOff), kmers_(kmers) {
}
KmerMapByPos::KmerMapByPos(size_t kLength, uint32_t runCutOff) : KmerMapBase(kLength, runCutOff){}
KmerMapByPos::KmerMapByPos(size_t kLength) : KmerMapBase(kLength){}
KmerMapByPos::KmerMapByPos() : KmerMapBase(){}



std::pair<uint32_t, uint32_t> KmerMapByPos::determineExpandPos(uint32_t expandSize,
		uint32_t pos, const std::string & seq) {
	uint32_t startPos = 0;
	uint32_t stopPos = 0;
	if (expandSize >= pos) {
		startPos = 0;
	} else {
		startPos = pos - expandSize;
	}
	if ((expandSize + pos + 1) > seq.size()) {
		stopPos = seq.size();
	} else {
		stopPos = pos + expandSize + 1;
	}
	return {startPos, stopPos};
}

void KmerMapByPos::increaseCounts(const std::string & str){
	increaseCounts(str,1);
}
void KmerMapByPos::increaseCounts(const std::string & str, uint32_t readCount){
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto currentKmer = str.substr(cursor, kLength_);
		auto search = kmers_[cursor].find(currentKmer);
		if (search == kmers_[cursor].end()) {
			auto currentKmerObj = kmer(currentKmer);
			currentKmerObj.increaseCnt(readCount);
			kmers_[cursor].emplace(currentKmer,currentKmerObj);
		} else {
			search->second.increaseCnt(readCount);
		}
		++cursor;
	}
}

void KmerMapByPos::increaseCountsExtra(const std::string & str){
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto currentKmer = str.substr(cursor, kLength_);
		auto search = kmers_[cursor].find(currentKmer);
		if (search == kmers_[cursor].end()) {
			kmers_[cursor].emplace(currentKmer, kmer(currentKmer, cursor));
		} else {
			search->second.addPosition(cursor);
		}
		++cursor;
	}
}

void KmerMapByPos::increaseCountsExtra(const std::string & str, const std::string & uid,
		uint32_t readCount){
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto currentKmer = str.substr(cursor, kLength_);
		auto search = kmers_[cursor].find(currentKmer);
		if (search == kmers_[cursor].end()) {
			kmers_[cursor].emplace(currentKmer, kmer(currentKmer, cursor, uid, readCount));
		} else {
			search->second.addPosition(cursor, uid, readCount);
		}
		++cursor;
	}
}

void KmerMapByPos::increaseCountsExp(const std::string & str, uint32_t expandSize){
	increaseCountsExp(str, expandSize,1);
}
void KmerMapByPos::increaseCountsExp(const std::string & str, uint32_t expandSize,
		uint32_t readCount){
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto positions = determineExpandPos(expandSize, cursor, str);
		auto currentKmer = str.substr(cursor, kLength_);
		for(auto pos : iter::range(positions.first, positions.second)){
			auto search = kmers_[pos].find(currentKmer);
			if (search == kmers_[pos].end()) {
				auto currentKmerObj = kmer(currentKmer);
				currentKmerObj.increaseCnt(readCount);
				kmers_[pos].emplace(currentKmer,currentKmerObj);
			} else {
				search->second.increaseCnt(readCount);
			}
		}
		++cursor;
	}
}

void KmerMapByPos::increaseCountsExtraExp(const std::string & str, uint32_t expandSize){
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto positions = determineExpandPos(expandSize, cursor, str);
		auto currentKmer = str.substr(cursor, kLength_);
		for(auto pos : iter::range(positions.first, positions.second)){
			auto search = kmers_[pos].find(currentKmer);
			if (search == kmers_[pos].end()) {
				kmers_[pos].emplace(currentKmer, kmer(currentKmer, pos));
			} else {
				search->second.addPosition(pos);
			}
		}
		++cursor;
	}
}
void KmerMapByPos::increaseCountsExtraExp(const std::string & str, uint32_t expandSize,
		const std::string & uid, uint32_t readCount){
	uint32_t cursor = 0;
	while (cursor + kLength_ <= static_cast<uint32_t>(str.size())) {
		auto positions = determineExpandPos(expandSize, cursor, str);
		auto currentKmer = str.substr(cursor, kLength_);
		for(auto pos : iter::range(positions.first, positions.second)){
			auto search = kmers_[pos].find(currentKmer);
			if (search == kmers_[pos].end()) {
				kmers_[pos].emplace(currentKmer, kmer(currentKmer, pos,uid, readCount));
			} else {
				search->second.addPosition(pos,uid, readCount);
			}
		}
		++cursor;
	}
}


bool KmerMapByPos::isKmerLowFreq(const kmerWithPos & k) const {
	return isKmerLowFreq(k, runCutOff_);
}

bool KmerMapByPos::isKmerLowFreq(const kmerWithPos & k, uint32_t runCutOff) const {
	auto posSearch = kmers_.find(k.kPos_);
	if(posSearch == kmers_.end()){
		return false;
	}else{
		auto kSearch = posSearch->second.find(k.kmer_);
		return (kSearch == posSearch->second.end() || kSearch->second.readCnt_ <= runCutOff);
	}
}

uint32_t KmerMapByPos::getKmerFreq(const kmerWithPos & k) const {
	auto posSearch = kmers_.find(k.kPos_);
	if (posSearch == kmers_.end()) {
		return 0;
	} else {
		auto kSearch = posSearch->second.find(k.kmer_);
		if (kSearch == posSearch->second.end()) {
			return 0;
		} else {
			return kSearch->second.readCnt_;
		}
	}
}

Json::Value KmerMapByPos::toJson() const {
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::KmerMapByPos");
	ret["super"] = KmerMapBase::toJson();
	ret["kmers_"] = bib::json::toJson(kmers_);
	return ret;
}


KmerMapByPos::~KmerMapByPos() {
}



KmerMaps::KmerMaps(size_t kLength, uint32_t runCutOff) :
		kLength_(kLength), runCutOff_(runCutOff), kmersNoPos_(
				std::make_shared<KmerMap>(kLength, runCutOff)), kmersByPos_(
				std::make_shared<KmerMapByPos>(kLength, runCutOff)) {
}
KmerMaps::KmerMaps(size_t kLength) :
		kLength_(kLength), runCutOff_(0), kmersNoPos_(
				std::make_shared<KmerMap>(kLength, 0)), kmersByPos_(
				std::make_shared<KmerMapByPos>(kLength, 0)) {
}
KmerMaps::KmerMaps() :
		kLength_(1), runCutOff_(0), kmersNoPos_(std::make_shared<KmerMap>(1, 0)), kmersByPos_(
				std::make_shared<KmerMapByPos>(1, 0)) {
	kmers_ = kmersByPos_;
	kmersByPosition_ = true;
}

void KmerMaps::setKmers(bool kmersByPosition){
	kmersByPosition_ = kmersByPosition;
	kmerSet_ = true;
	if(kmersByPosition){
		kmers_ = kmersByPos_;
	}else{
		kmers_ = kmersNoPos_;
	}
}

bool KmerMaps::isKmerLowFreq(const kmerWithPos & k) const {
	return kmers_->isKmerLowFreq(k, runCutOff_);
}

bool KmerMaps::isKmerLowFreq(const kmerWithPos & k, uint32_t runCutOff) const {
	return kmers_->isKmerLowFreq(k, runCutOff);
}

void KmerMaps::setAllRunCutOffs(uint32_t runCutOff) {
	runCutOff_ = runCutOff;
	if (kmersByPos_) {
		kmersByPos_->runCutOff_ = runCutOff;
	}
	if (kmersNoPos_) {
		kmersNoPos_->runCutOff_ = runCutOff;
	}
}


void KmerMaps::resetMaps(){
	kmersNoPos_ = std::make_shared<KmerMap>(kLength_, runCutOff_);
	kmersByPos_ = std::make_shared<KmerMapByPos>(kLength_, runCutOff_);
	if(kmersByPosition_){
		kmers_ = kmersByPos_;
	}else{
		kmers_ = kmersNoPos_;
	}
}

Json::Value KmerMaps::toJson() const {
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::KmerMaps");
	ret["kmersByPos_"] = bib::json::toJson(*kmersByPos_);
	ret["kmersNoPos_"] = bib::json::toJson(*kmersNoPos_);
	ret["kmersByPos_"] = bib::json::toJson(kmersByPosition_);
	ret["kLength_"] = bib::json::toJson(kLength_);
	ret["runCutOff_"] = bib::json::toJson(runCutOff_);
	ret["kLength_"] = bib::json::toJson(kLength_);

	return ret;
}





}  // namespace bibseq
