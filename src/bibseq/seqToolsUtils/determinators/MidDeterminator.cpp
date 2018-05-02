/*
 * MidDeterminator.cpp
 *
 *  Created on: Jun 10, 2015
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
#include "MidDeterminator.hpp"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/IO/SeqIO.h"

namespace bibseq {

MidDeterminator::midPos::midPos() :
		midName_("unrecognized"), midPos_(midPos::npos), barcodeSize_(0), barcodeScore_(
				0) {
}

MidDeterminator::midPos::midPos(const std::string & midName, uint64_t midPos,
		uint64_t barcodeSize, uint32_t barcodeScore) :
		midName_(midName), midPos_(midPos), barcodeSize_(barcodeSize), barcodeScore_(
				barcodeScore) {
}

Json::Value MidDeterminator::midPos::toJson() const {
	Json::Value ret;
	ret["midName_"] = bib::json::toJson(midName_);
	ret["midPos_"] = bib::json::toJson(midPos_);
	ret["barcodeSize_"] = bib::json::toJson(barcodeSize_);
	ret["barcodeScore_"] = bib::json::toJson(barcodeScore_);
	ret["inRevComp_"] = bib::json::toJson(inRevComp_);
	return ret;
}

MidDeterminator::midPos::operator bool() const {
	return npos != midPos_;
}



std::string MidDeterminator::midPos::getFailureCaseName(FailureCase fCase){
	std::string fCaseStr = "unhandledCase";
	switch (fCase) {
	case MidDeterminator::midPos::FailureCase::MISMATCHING_DIRECTION:
		fCaseStr = "MISMATCHING_DIRECTION";
		break;
	case MidDeterminator::midPos::FailureCase::MISMATCHING_MIDS:
		fCaseStr = "MISMATCHING_MIDS";
		break;
	case MidDeterminator::midPos::FailureCase::NONE:
		fCaseStr = "NONE";
		break;
	case MidDeterminator::midPos::FailureCase::NO_MATCHING:
		fCaseStr = "NO_MATCHING";
		break;
	case MidDeterminator::midPos::FailureCase::TOO_MANY_MATCHING:
		fCaseStr = "TOO_MANY_MATCHING";
		break;
	default:
		fCaseStr = "unhandledCase";
		break;
	}
	return fCaseStr;
}

VecStr MidDeterminator::midPos::getFailureCaseNames() {
	return VecStr { "NONE", "NO_MATCHING", "TOO_MANY_MATCHING",
			"MISMATCHING_MIDS", "MISMATCHING_DIRECTION", "unhandledCase" };
}


double MidDeterminator::midPos::normalizeScoreByLen() const{
	return static_cast<double>(barcodeScore_)/barcodeSize_;
}

MidDeterminator::MidInfo::MidInfo(const std::string & midName,
		const std::string & barcode) :
		midName_(midName), bar_(std::make_unique<motif>(barcode)), rcompBar_(
				std::make_unique<motif>(seqUtil::reverseComplement(barcode, "DNA"))) {
	if (barcode.size() > 1) {
		shortenFrontBar_ = std::make_unique<motif>(barcode.substr(1));
		shortenBackBar_ = std::make_unique<motif>(barcode.substr(0, barcode.size() - 1));
		auto rcBar = seqUtil::reverseComplement(barcode, "DNA");
		shortenFrontRCompBar_ = std::make_unique<motif>(rcBar.substr(1));
		shortenBackRCompBar_ = std::make_unique<motif>(rcBar.substr(0, rcBar.size() - 1));
	}
	rCompSame_ = bar_->motifOriginal_ == rcompBar_->motifOriginal_;
}

MidDeterminator::MidDeterminator(const table & mids) {
	if (!bib::in(std::string("id"), mids.columnNames_)
			|| !bib::in(std::string("barcode"), mids.columnNames_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error in creating MidDeterminator, need to have at "
						"least two columns, id and barcode, only have "
						+ bib::conToStr(mids.columnNames_, ",");
		throw std::runtime_error {ss.str()};
	}
	for (const auto & row : mids.content_) {
		addBarcode(row[mids.getColPos("id")],row[mids.getColPos("barcode")]);
	}
}

MidDeterminator::MidDeterminator(const std::unordered_map<std::string, MidDeterminator::MidInfo> & mids){
	for (const auto & mid : mids) {
		addBarcode(mid.second.midName_,mid.second.bar_->motifOriginal_);
	}
}

MidDeterminator::MidDeterminator() {

}


void MidDeterminator::setAllowableMismatches(uint32_t allowableMismatches) {
	allowableMismatches_ = allowableMismatches;
}

void MidDeterminator::setMidEndsRevComp(bool midEndsRevComp){
	midEndsRevComp_ = midEndsRevComp;
}



bool MidDeterminator::containsMidByName(const std::string & name) const {
	return bib::in(name, mids_);
}

bool MidDeterminator::containsMidByBarcode(const std::string & barcode)const{
	for(const auto & mid : mids_){
		if(mid.second.bar_->motifOriginal_ == barcode){
			return true;
		}
	}
	return false;
}

std::string MidDeterminator::getMidName(const std::string & barcode)const{
	for(const auto & mid : mids_){
		if(mid.second.bar_->motifOriginal_ == barcode){
			return mid.first;
		}
	}
	return "no_name_for_barcode:" + barcode;
}

void MidDeterminator::addBarcode(const std::string & name, const std::string & barcode) {
	if (containsMidByName(name)) {
		std::stringstream ss;
		ss << "Error, MidDeterminator already contains mid by name: " << name
				<< "\n";
		ss << "original barcode: " << mids_.at(name).bar_->motifOriginal_
				<< ", adding barcode: " << barcode << "\n";
		throw std::runtime_error { bib::bashCT::boldRed(ss.str()) };
	}
	if (containsMidByBarcode(barcode)) {
		std::stringstream ss;
		ss << "Error, MidDeterminator already contains mid by barcode: "
				<< barcode << "\n";
		ss << "original name: " << getMidName(barcode)
				<< ", adding name: " << name << "\n";
		throw std::runtime_error { bib::bashCT::boldRed(ss.str()) };
	}
	mids_.emplace(name,
			MidInfo { name, barcode});
}







std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPos(
		const std::string & seq, uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.bar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStop = std::min<uint32_t>(within + mid.second.bar_->size(),
				seq.size());
		auto positions = mid.second.bar_->findPositionsFull(seq,
				allowableMismatches_, 0, searchStop);
		if (!positions.empty()) {
			auto currentScore = mid.second.bar_->scoreMotif(
					seq.substr(positions.front(),
							mid.second.bar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.front(),
					mid.second.bar_->motifOriginal_.size(), currentScore);
		}
	}
	return possible;
}



std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPosBack(
		const std::string & seq, uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.bar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStart = seq.size() - mid.second.bar_->size();
		if (within > searchStart) {
			searchStart = 0;
		} else {
			searchStart -= within;
		}
		auto positions = mid.second.bar_->findPositionsFull(seq,
				allowableMismatches_, searchStart, seq.size());
		if (!positions.empty()) {
			auto currentScore = mid.second.bar_->scoreMotif(
					seq.substr(positions.back(),
							mid.second.bar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.back(),
					mid.second.bar_->motifOriginal_.size(), currentScore);
		}
	}
	return possible;
}

std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPosComp(const std::string & seq,
		uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.rcompBar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStop = std::min<uint32_t>(within + mid.second.rcompBar_->size(), seq.size());
		auto positions = mid.second.rcompBar_->findPositionsFull(seq, allowableMismatches_, 0,
				searchStop);
		if (!positions.empty()) {
			auto currentScore = mid.second.rcompBar_->scoreMotif(
					seq.substr(positions.front(), mid.second.rcompBar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.front(),
					mid.second.rcompBar_->motifOriginal_.size(), currentScore);
			possible.back().inRevComp_ = true;
		}
	}
	return possible;
}



std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPosCompBack(
		const std::string & seq, uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.rcompBar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStart = seq.size() - mid.second.rcompBar_->size();
		if (within > searchStart) {
			searchStart = 0;
		} else {
			searchStart -= within;
		}
		auto positions = mid.second.rcompBar_->findPositionsFull(seq,
				allowableMismatches_, searchStart, seq.size());
		if (!positions.empty()) {
			auto currentScore = mid.second.rcompBar_->scoreMotif(
					seq.substr(positions.back(),
							mid.second.rcompBar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.back(),
								mid.second.rcompBar_->motifOriginal_.size(),
								currentScore);
			possible.back().inRevComp_ = true;
		}
	}
	return possible;
}

std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPosShorten(
		const std::string & seq, uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.shortenFrontBar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStop = std::min<uint32_t>(within + mid.second.shortenFrontBar_->size(),
				seq.size());
		auto positions = mid.second.shortenFrontBar_->findPositionsFull(seq,
				allowableMismatches_, 0, searchStop);
		if (!positions.empty()) {
			auto currentScore = mid.second.shortenFrontBar_->scoreMotif(
					seq.substr(positions.front(),
							mid.second.shortenFrontBar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.front(),
					mid.second.shortenFrontBar_->motifOriginal_.size(), currentScore);
		}
	}
	return possible;
}



std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPosBackShorten(
		const std::string & seq, uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.shortenBackBar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStart = seq.size() - mid.second.shortenBackBar_->size();
		if (within > searchStart) {
			searchStart = 0;
		} else {
			searchStart -= within;
		}
		auto positions = mid.second.shortenBackBar_->findPositionsFull(seq,
				allowableMismatches_, searchStart, seq.size());
		if (!positions.empty()) {
			auto currentScore = mid.second.shortenBackBar_->scoreMotif(
					seq.substr(positions.back(),
							mid.second.shortenBackBar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.back(),
					mid.second.shortenBackBar_->motifOriginal_.size(), currentScore);
		}
	}
	return possible;
}

std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPosCompShorten(const std::string & seq,
		uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.shortenFrontRCompBar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStop = std::min<uint32_t>(within + mid.second.shortenFrontRCompBar_->size(), seq.size());
		auto positions = mid.second.shortenFrontRCompBar_->findPositionsFull(seq, allowableMismatches_, 0,
				searchStop);
		if (!positions.empty()) {
			auto currentScore = mid.second.shortenFrontRCompBar_->scoreMotif(
					seq.substr(positions.front(), mid.second.shortenFrontRCompBar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.front(),
					mid.second.shortenFrontRCompBar_->motifOriginal_.size(), currentScore);
			possible.back().inRevComp_ = true;
		}
	}
	return possible;
}



std::vector<MidDeterminator::midPos> MidDeterminator::determinePossibleMidPosCompBackShorten(
		const std::string & seq, uint32_t within) {
	std::vector<midPos> possible;
	for (auto & mid : mids_) {
		if (seq.size() < mid.second.shortenBackRCompBar_->motifOriginal_.size()
				|| allowableMismatches_ >= seq.size()) {
			continue;
		}
		uint32_t searchStart = seq.size() - mid.second.shortenBackRCompBar_->size();
		if (within > searchStart) {
			searchStart = 0;
		} else {
			searchStart -= within;
		}
		auto positions = mid.second.shortenBackRCompBar_->findPositionsFull(seq,
				allowableMismatches_, searchStart, seq.size());
		if (!positions.empty()) {
			auto currentScore = mid.second.shortenBackRCompBar_->scoreMotif(
					seq.substr(positions.back(),
							mid.second.shortenBackRCompBar_->motifOriginal_.size()));
			possible.emplace_back(mid.first, positions.back(),
								mid.second.shortenBackRCompBar_->motifOriginal_.size(),
								currentScore);
			possible.back().inRevComp_ = true;
		}
	}
	return possible;
}


void MidDeterminator::processInfoWithMidPos(seqInfo & info,
		const MidDeterminator::midPos & frontPos,
		const MidDeterminator::midPos & backPos) {
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//	std::cout << "backPos.midName_: " << backPos.midName_ << std::endl;
//	std::cout << "frontPos.midPos_: " << frontPos.midName_ << std::endl;
//	std::cout << "backPos.midPos_: " << backPos.midPos_ << std::endl;
//	std::cout << "frontPos.midPos_: " << frontPos.midPos_ << std::endl;
//	std::cout << "frontPos.midPos_ + frontPos.barcodeSize_: " << frontPos.midPos_ + frontPos.barcodeSize_ << std::endl;
//	std::cout << " len(info): " << len(info) << std::endl;
	info.trimBack(backPos.midPos_);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	info.trimFront(frontPos.midPos_ + frontPos.barcodeSize_);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if(frontPos.inRevComp_){
		info.reverseComplementRead(true, true);
	}
}

void MidDeterminator::processInfoWithMidPos(seqInfo & info, const MidDeterminator::midPos & pos,
		MidDeterminator::MidDeterminePars pars) {
	/**@todo remove the same size other side, should add a better check for the other side barcode
	 */

	if (pars.variableStop_ > pos.midPos_ || (0 == pars.variableStop_ && 0 == pos.midPos_)) {
		info.trimFront(pos.midPos_ + pos.barcodeSize_);
		if (pars.barcodesBothEnds_) {
			info.trimBack(len(info) - mids_.at(pos.midName_).bar_->size());
		}
	} else {
		info.trimBack(pos.midPos_);
		if (pars.barcodesBothEnds_) {
			info.trimFront(mids_.at(pos.midName_).bar_->size());
		}
	}
	if(pos.inRevComp_){
		info.reverseComplementRead(true, true);
	}
}



void MidDeterminator::processInfoWithMidPos(PairedRead & seq,
		const MidDeterminator::midPos & firstPos,
		const MidDeterminator::midPos & matePos) {
	if(seq.mateRComplemented_){
		seq.mateSeqBase_.trimBack(matePos.midPos_);
	}else{
		seq.mateSeqBase_.trimFront(matePos.midPos_ + matePos.barcodeSize_);
	}
	seq.seqBase_.trimFront(firstPos.midPos_ + firstPos.barcodeSize_);
	if(firstPos.inRevComp_){
		seq.seqBase_.reverseComplementRead(true, true);
		seq.mateSeqBase_.reverseComplementRead(true, true);
	}
}

void MidDeterminator::processInfoWithMidPos(PairedRead & seq,
		const MidDeterminator::midPos & pos,
		MidDeterminator::MidDeterminePars pars) {
	/**@todo remove the same size other side, should add a better check for the other side barcode
	 */


}





std::pair<MidDeterminator::midPos, MidDeterminator::midPos>  MidDeterminator::fullDetermine(PairedRead & seq, MidDeterminePars pars){

//	if(seq.mateRComplemented_){
//
//	}
//	std::vector<midPos> frontPossible = determinePossibleMidPos(seq.seqBase_.seq_, pars.variableStop_);
//	std::vector<midPos> backPossible = determinePossibleMidPosComp(seq.mateSeqBase_.seq_, pars.variableStop_);
//
//
//
//	if (pars.barcodesBothEnds_) {
//		addOtherVec(possible,
//				determinePossibleMidPosBack(info.seq_, pars.variableStop_));
//	}
//
//	if (pars.checkComplement_) {
//		addOtherVec(possible,
//				determinePossibleMidPosCompBack(info.seq_, pars.variableStop_));
//	}
//
//	if (pars.checkComplement_ && pars.barcodesBothEnds_) {
//		addOtherVec(possible,
//				determinePossibleMidPosComp(info.seq_, pars.variableStop_));
//	}
//
//	if (pars.checkForShorten_) {
//		if (pars.barcodesBothEnds_) {
//			addOtherVec(possible, determinePossibleMidPosShorten(info.seq_, 0));
//		}
//
//		if (pars.barcodesBothEnds_) {
//			addOtherVec(possible, determinePossibleMidPosBackShorten(info.seq_, 0));
//		}
//
//		if (pars.checkComplement_) {
//			addOtherVec(possible,
//					determinePossibleMidPosCompBackShorten(info.seq_, 0));
//		}
//
//		if (pars.checkComplement_ && pars.barcodesBothEnds_) {
//			addOtherVec(possible, determinePossibleMidPosCompShorten(info.seq_, 0));
//		}
//	}
//
//



	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	/**@todo expensive way to do this for now */
	uint32_t spacerSize = 20;
	seqInfo tempSeq = seq.seqBase_;
	tempSeq.append(std::string(spacerSize, 'N'));
	if(seq.mateRComplemented_){
		tempSeq.append(seq.mateSeqBase_.seq_);
	}else{
		tempSeq.append(seqUtil::reverseComplement(seq.mateSeqBase_.seq_, "DNA"));
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	size_t lengthOfTempSeq = len(tempSeq);
	auto positions = fullDetermine(tempSeq, pars);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  if(positions.first){

		if (positions.second) {
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			if(seq.mateRComplemented_){
				positions.second.midPos_ = positions.second.midPos_ - len(seq.seqBase_) - spacerSize;
			}else{
				positions.second.midPos_ = lengthOfTempSeq - positions.second.midPos_ - positions.second.barcodeSize_;
			}

			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			processInfoWithMidPos(seq, positions.first, positions.second);
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		} else {
			//processInfoWithMidPos(seq, positions.first, pars);
			if (pars.variableStop_ > positions.first.midPos_ || (0 == pars.variableStop_ && 0 == positions.first.midPos_)) {
				seq.seqBase_.trimFront(positions.first.midPos_ + positions.first.barcodeSize_);
				if (pars.barcodesBothEnds_) {
					if(seq.mateRComplemented_){
						seq.mateSeqBase_.trimBack(len(seq.mateSeqBase_) - mids_.at(positions.first.midName_).bar_->size());
					}else{
						seq.mateSeqBase_.trimFront(mids_.at(positions.first.midName_).bar_->size());
					}
				}
			} else {
				if (seq.mateRComplemented_) {
					positions.first.midPos_ = positions.first.midPos_
							- len(seq.seqBase_) - spacerSize;
					seq.mateSeqBase_.trimBack(positions.first.midPos_);
				} else {
					positions.first.midPos_ = lengthOfTempSeq - positions.first.midPos_ - positions.first.barcodeSize_;
					seq.mateSeqBase_.trimFront(positions.first.midPos_ + positions.first.barcodeSize_);
				}
				if (pars.barcodesBothEnds_) {
					/**@todo look into to see if this actually needed */
					seq.seqBase_.trimFront(mids_.at(positions.first.midName_).bar_->size());
				}
			}
			if(positions.first.inRevComp_){
				//seq.seqBase_.reverseComplementRead(true, true);
				//seq.mateSeqBase_.reverseComplementRead(true, true);
				auto tempSeq = seq.seqBase_;
				seq.seqBase_ = seq.mateSeqBase_;
				seq.mateSeqBase_ = tempSeq;
				seq.seqBase_.name_ += "_Comp";
				seq.mateSeqBase_.name_ += "_Comp";
			}
		}
  }
  return positions;
}

std::pair<MidDeterminator::midPos, MidDeterminator::midPos> MidDeterminator::fullDetermine(
		seqInfo & info, MidDeterminePars pars) {
	midPos ret;
	midPos backPos;
	//determine possible mid positions;
	std::vector<midPos> possible = determinePossibleMidPos(info.seq_, pars.variableStop_);
	if (pars.barcodesBothEnds_) {
		addOtherVec(possible,
				determinePossibleMidPosBack(info.seq_, pars.variableStop_));
	}
	if (pars.checkComplement_) {
		addOtherVec(possible,
				determinePossibleMidPosCompBack(info.seq_, pars.variableStop_));
	}
	if (pars.checkComplement_ && pars.barcodesBothEnds_) {
		addOtherVec(possible,
				determinePossibleMidPosComp(info.seq_, pars.variableStop_));
	}
	if (pars.checkForShorten_) {
		if (pars.barcodesBothEnds_) {
			addOtherVec(possible, determinePossibleMidPosShorten(info.seq_, 0));
		}

		if (pars.barcodesBothEnds_) {
			addOtherVec(possible, determinePossibleMidPosBackShorten(info.seq_, 0));
		}

		if (pars.checkComplement_) {
			addOtherVec(possible,
					determinePossibleMidPosCompBackShorten(info.seq_, 0));
		}

		if (pars.checkComplement_ && pars.barcodesBothEnds_) {
			addOtherVec(possible, determinePossibleMidPosCompShorten(info.seq_, 0));
		}
	}

	if (1 == possible.size()) {
		ret = possible.front();
	} else if (!possible.empty()) {
		//determine best
		std::vector<midPos> best;
		double bestScore = 0;
		for (const auto & poss : possible) {
			if (poss.normalizeScoreByLen() > bestScore) {
				bestScore = poss.normalizeScoreByLen();
				best.clear();
				best.emplace_back(poss);
			} else if (poss.normalizeScoreByLen() == bestScore && 0 != bestScore) {
				best.emplace_back(poss);
			}
		}
		if (1 == best.size()) {
			ret = best.front();
		} else if (2 == best.size()) {
			if (best[0].midName_ == best[1].midName_) {
				/**@todo odd way to allow for rev comp on one end and regular on the other if poeple's library's are like that*/
				/**@todo also this doesn't handle if both aren't the best, like if one had a mismatch but the other didn't */
				if (mids_.at(best[0].midName_).rCompSame_) {
					if (best[0].midPos_ == best[1].midPos_) {
						ret = best[0];
					} else {
						if (best[0].midPos_ < best[1].midPos_) {
							ret = best[0];
							backPos = best[1];
						} else {
							ret = best[1];
							backPos = best[0];
						}
					}
				} else {
					if (midEndsRevComp_) {
						if (best[0].inRevComp_ != best[1].inRevComp_) {
							if (best[0].midPos_ < best[1].midPos_) {
								ret = best[0];
								backPos = best[1];
							} else {
								ret = best[1];
								backPos = best[0];
							}
						} else {
							ret.fCase_ = midPos::FailureCase::MISMATCHING_DIRECTION;
						}
					} else {
						if (best[0].inRevComp_ == best[1].inRevComp_) {
							if (best[0].midPos_ < best[1].midPos_) {
								ret = best[0];
								backPos = best[1];
							} else {
								ret = best[1];
								backPos = best[0];
							}
						} else {
							ret.fCase_ = midPos::FailureCase::MISMATCHING_DIRECTION;
						}
					}
				}
			} else {
				//found two mids with different names
				ret.fCase_ = midPos::FailureCase::MISMATCHING_MIDS;
			}
		} else if (best.size() > 2) {
			if (std::all_of(best.begin() + 1, best.end(),
					[&best](const midPos & pos) {return pos.midName_ == best.front().midName_;})
					&& mids_.at(best[0].midName_).rCompSame_) {
				std::sort(best.begin(), best.end(),
						[](const midPos & first, const midPos & second) {
							if(first.midPos_ == second.midPos_) {
								return first.inRevComp_ > second.inRevComp_;
							} else {
								return first.midPos_ < second.midPos_;
							}
						});
				if (best.size() == 3 || best.size() == 4) {
					ret = best[0];
					backPos = best[3];
				} else {
					//found more than one 4 mids with all the same name
					ret.fCase_ = midPos::FailureCase::TOO_MANY_MATCHING;
				}
			} else {
				//found more than two mids and they aren't the same name and their rev comp aren't the same as their forward
				ret.fCase_ = midPos::FailureCase::TOO_MANY_MATCHING;
			}
		} else {
			//best is empty... i think, if i'm following logic correctly
			ret.fCase_ = midPos::FailureCase::NO_MATCHING;
		}
	} else {
		ret.fCase_ = midPos::FailureCase::NO_MATCHING;
	}

	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if (ret) {
		if (backPos) {
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			if(ret.midPos_ < pars.variableStop_ && backPos.midPos_ < pars.variableStop_){
				processInfoWithMidPos(info, backPos, pars);
				ret = backPos;
				backPos = midPos();
			}else if(ret.midPos_ + pars.variableStop_ >= len(info) &&  backPos.midPos_ + pars.variableStop_ >= len(info)){
				processInfoWithMidPos(info, ret, pars);
				backPos = midPos();
			}else{
				processInfoWithMidPos(info, ret, backPos);
			}
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		} else {
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			processInfoWithMidPos(info, ret, pars);
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	return {ret, backPos};
}

void MidDeterminator::increaseFailedBarcodeCounts(
		const MidDeterminator::midPos & pos,
		std::unordered_map<std::string, uint32_t> & counts) {
	std::string fCaseStr = midPos::getFailureCaseName(pos.fCase_);
	counts[fCaseStr] += 1;
}

} /* namespace bibseq */
