#include "readChecker.hpp"
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

ReadChecker::ReadChecker(std::string markWith, bool mark) :
		markWith_(std::move(markWith)), mark_(mark) {
}

void ReadChecker::markName(seqInfo & info) const {
	if (mark_) {
		info.name_.append(markWith_);
	}
}

bool ReadChecker::checkRead(seqInfo & info) const {
	throw std::runtime_error { std::string(__PRETTY_FUNCTION__)
			+ ": ReadChecker base class should not be used it self" };
}

bool ReadChecker::checkRead(PairedRead & seq) const {
	if(checkRead(seq.seqBase_)){
		if(!checkRead(seq.mateSeqBase_)){
			seq.seqBase_.on_ = false;
			markName(seq.seqBase_);
		}else{
			return true;
		}
	}
	return false;
}

ReadChecker::~ReadChecker() {
}

ReadCheckerLenWithin::ReadCheckerLenWithin(uint32_t basesWithin,
		double givenLen, bool mark) :
		ReadChecker(
				std::string(
						"_lenMoreThan" + std::to_string(basesWithin) + "_From"
								+ std::to_string(givenLen)), mark), basesWithin_(basesWithin), givenLen_(
				givenLen) {
}

bool ReadCheckerLenWithin::checkRead(seqInfo & info) const {
	if (uAbsdiff(givenLen_, info.seq_.length()) > basesWithin_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

ReadCheckerLenWithin::~ReadCheckerLenWithin() {
}

ReadCheckerLenBelow::ReadCheckerLenBelow(uint32_t maxLen, bool mark) :
		ReadChecker(std::string("_length>" + std::to_string(maxLen)), mark), maxLen_(
				maxLen) {
}

bool ReadCheckerLenBelow::checkRead(seqInfo & info) const {
	if (info.seq_.length() > maxLen_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

ReadCheckerLenBelow::~ReadCheckerLenBelow() {
}

ReadCheckerLenAbove::ReadCheckerLenAbove(uint32_t minLen, bool mark) :
		ReadChecker(std::string("_length<" + std::to_string(minLen)), mark), minLen_(
				minLen) {
}

bool ReadCheckerLenAbove::checkRead(seqInfo & info) const {
	if (len(info) < minLen_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}




ReadCheckerLenAbove::~ReadCheckerLenAbove() {
}

ReadCheckerLenBetween::ReadCheckerLenBetween(uint32_t maxLen, uint32_t minLen,
		bool mark) :
		ReadChecker(std::string("readBetweenLens"), mark), minLen_(minLen), maxLen_(
				maxLen) {
}

bool ReadCheckerLenBetween::checkRead(seqInfo & info) const {
	std::string markWithMin = "_length<" + std::to_string(minLen_);
	std::string markWithMax = "_length>" + std::to_string(maxLen_);
	if (info.seq_.length() < minLen_) {
		info.on_ = false;
		if (mark_)
			info.name_.append(markWithMin);
		return false;
	} else if (info.seq_.length() > maxLen_) {
		info.on_ = false;
		if (mark_)
			info.name_.append(markWithMax);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

ReadCheckerLenBetween::~ReadCheckerLenBetween() {
}

ReadCheckerQualCheck::ReadCheckerQualCheck(uint32_t qualCutOff,
		double qualFracCutOff, bool mark) :
		ReadChecker(
				std::string(
						"_q" + estd::to_string(qualCutOff) + "<"
								+ std::to_string(qualFracCutOff)), mark), qualCutOff_(
				qualCutOff), qualFracCutOff_(qualFracCutOff) {
}

bool ReadCheckerQualCheck::checkRead(PairedRead & info) const{
	if (info.getQualCheck(qualCutOff_) < qualFracCutOff_) {
		info.seqBase_.on_ = false;
		info.mateSeqBase_.on_ = false;
		markName(info.seqBase_);
		markName(info.mateSeqBase_);
		return false;
	} else {
		info.seqBase_.on_ = true;
		info.mateSeqBase_.on_ = true;
		return true;
	}
}


bool ReadCheckerQualCheck::checkRead(seqInfo & info) const {
	if (info.getQualCheck(qualCutOff_) < qualFracCutOff_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}



ReadCheckerQualCheck::~ReadCheckerQualCheck() {
}

ReadCheckerOnCount::ReadCheckerOnCount(double countCutOff, bool mark) :
		ReadChecker(std::string("_cnt<=" + estd::to_string(countCutOff)), mark), countCutOff_(
				countCutOff) {
}

bool ReadCheckerOnCount::checkRead(seqInfo & info) const {
	if (info.cnt_ <= countCutOff_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool ReadCheckerOnCount::checkRead(PairedRead & seq) const {
	return checkRead(seq.seqBase_);
}



ReadCheckerOnCount::~ReadCheckerOnCount() {
}

ReadCheckerOnFrac::ReadCheckerOnFrac(double fracCutOff, bool mark) :
		ReadChecker(std::string("_frac<=" + estd::to_string(fracCutOff)), mark), fracCutOff_(
				fracCutOff) {
}

bool ReadCheckerOnFrac::checkRead(seqInfo & info) const {
	if (info.frac_ <= fracCutOff_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool ReadCheckerOnFrac::checkRead(PairedRead & seq) const {
	return checkRead(seq.seqBase_);
}

ReadCheckerOnFrac::~ReadCheckerOnFrac() {
}

ReadCheckerOnNucComp::ReadCheckerOnNucComp(charCounter counter, double fracDiff,
		bool mark) :
		ReadChecker(std::string("_nucCompFrac>" + estd::to_string(fracDiff)), mark), counter_(
				std::move(counter)), fracDiff_(fracDiff) {

}

bool ReadCheckerOnNucComp::checkRead(seqInfo & info) const {
	charCounter count(counter_.alphabet_);
	count.increaseCountByString(info.seq_);
	count.setFractions();
	if (counter_.getFracDifference(count, counter_.alphabet_) > fracDiff_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool ReadCheckerOnNucComp::checkRead(PairedRead & seq) const {
	charCounter count(counter_.alphabet_);
	count.increaseCountByString(seq.seqBase_.seq_);
	if(seq.mateRComplemented_){
		count.increaseCountByString(seq.mateSeqBase_.seq_);
	}else{
		count.increaseCountByString(seqUtil::reverseComplement(seq.mateSeqBase_.seq_, "DNA"));
	}
	count.setFractions();
	if (counter_.getFracDifference(count, counter_.alphabet_) > fracDiff_) {
		seq.seqBase_.on_ = false;
		markName(seq.seqBase_);
		return false;
	} else {
		seq.seqBase_.on_ = true;
		return true;
	}
}


ReadCheckerOnNucComp::~ReadCheckerOnNucComp() {

}

ReadCheckerOnSeqContaining::ReadCheckerOnSeqContaining(std::string str,
		uint32_t occurences, bool mark) :
		ReadChecker(std::string("_seqContains" + estd::to_string(str)), mark), str_(
				std::move(str)), occurences_(occurences) {
}

bool ReadCheckerOnSeqContaining::checkRead(seqInfo & info) const {

	std::string lowerCaseSearch = stringToLowerReturn(str_);
	if (countOccurences(info.seq_, str_)
			+ countOccurences(info.seq_, lowerCaseSearch) >= occurences_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool ReadCheckerOnSeqContaining::checkRead(PairedRead & seq) const {
	if(checkRead(seq.seqBase_)){
		if(!checkRead(seq.mateSeqBase_)){
			seq.seqBase_.on_ = false;
			markName(seq.seqBase_);
		}else{
			return true;
		}
	}
	return false;
}

ReadCheckerOnSeqContaining::~ReadCheckerOnSeqContaining() {
}

ReadCheckerOnNameContaining::ReadCheckerOnNameContaining(std::string str,
		bool mark) :
		ReadChecker(std::string("_nameContains" + str), mark), str_(std::move(str)) {
}

bool ReadCheckerOnNameContaining::checkRead(seqInfo & info) const {
	if (std::string::npos != info.name_.find(str_)) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool ReadCheckerOnNameContaining::checkRead(PairedRead & seq) const {
	return checkRead(seq.seqBase_);
}



ReadCheckerOnNameContaining::~ReadCheckerOnNameContaining() {
}

ReadCheckerOnQualityWindow::ReadCheckerOnQualityWindow(
		uint32_t qualityWindowLength, uint32_t qualityWindowStep,
		uint32_t qualityWindowThres, bool mark) :
		ReadChecker(
				njh::err::F() << "_failedWindowOf_wl:" << qualityWindowLength << ",ws:"
						<< qualityWindowStep << ",wt:" << qualityWindowThres, mark), qualityWindowLength_(
				qualityWindowLength), qualityWindowStep_(qualityWindowStep), qualityWindowThres_(
				qualityWindowThres) {
}

bool ReadCheckerOnQualityWindow::checkRead(seqInfo & info) const {
	if (!seqUtil::checkQualityWindow(qualityWindowLength_, qualityWindowThres_,
			qualityWindowStep_, info.qual_)) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

ReadCheckerOnQualityWindow::~ReadCheckerOnQualityWindow() {
}

ReadCheckerOnQualityWindowTrim::ReadCheckerOnQualityWindowTrim(
		uint32_t qualityWindowLength, uint32_t qualityWindowStep,
		uint32_t qualityWindowThres, uint32_t minLen, bool mark) :
		ReadChecker(
				njh::err::F() << "_failedWindowOf_wl:" << qualityWindowLength << ",ws:"
						<< qualityWindowStep << ",wt:" << qualityWindowThres, mark), qualityWindowLength_(
				qualityWindowLength), qualityWindowStep_(qualityWindowStep), qualityWindowThres_(
				qualityWindowThres), minLen_(minLen) {
}

bool ReadCheckerOnQualityWindowTrim::checkRead(seqInfo & info) const {
	auto pos = seqUtil::checkQualityWindowPos(qualityWindowLength_,
			qualityWindowThres_, qualityWindowStep_, info.qual_);
	if (pos <= minLen_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.setClip(0, pos - 1);
		info.on_ = true;
		return true;
	}
}

ReadCheckerOnQualityWindowTrim::~ReadCheckerOnQualityWindowTrim() {
}

ReadCheckerOnNs::ReadCheckerOnNs(bool mark) :
		ReadCheckerOnSeqContaining("N", 1, mark) {
}

ReadCheckerOnNs::~ReadCheckerOnNs() {
}

ReadCheckerOnKmerComp::ReadCheckerOnKmerComp(kmerInfo compareInfo,
		uint32_t kLength, double kmerCutoff, bool mark) :
		ReadChecker(
				njh::err::F() << "_failedKmerCutOff:" << kmerCutoff << "_kLen_"
						<< kLength, mark), compareInfo_(std::move(compareInfo)), kLength_(
				kLength), kmerCutoff_(kmerCutoff) {
}

bool ReadCheckerOnKmerComp::checkRead(seqInfo & info) const {
	kmerInfo currentInfo(info.seq_, kLength_, false);
	auto dist = compareInfo_.compareKmers(currentInfo);
	if (dist.second < kmerCutoff_) {
		info.on_ = false;
		markName(info);
		return false;
	} else {
		info.on_ = true;
		return true;
	}
}

bool ReadCheckerOnKmerComp::checkRead(PairedRead & seq) const {
	/**@todo this might not be the best way to do this, doing some positional checks or something would be a little better
	 *
	 */
	if(checkRead(seq.seqBase_)){
		std::unique_ptr<kmerInfo> kInfo;
		if(seq.mateRComplemented_){
			kInfo = std::make_unique<kmerInfo>(seq.mateSeqBase_.seq_, kLength_, false);
		}else{
			kInfo = std::make_unique<kmerInfo>(seqUtil::reverseComplement(seq.mateSeqBase_.seq_, "DNA"), kLength_, false);
		}
		auto dist = compareInfo_.compareKmers(*kInfo);
		if(dist.second < kmerCutoff_){
			seq.seqBase_.on_ = false;
			markName(seq.seqBase_);
		}else{
			return true;
		}
	}
	return false;
}

//

ReadCheckerOnKmerComp::~ReadCheckerOnKmerComp() {
}

}  // namespace njhseq
