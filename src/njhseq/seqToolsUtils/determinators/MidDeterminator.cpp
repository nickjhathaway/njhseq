/*
 * MidDeterminator.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */
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
#include "MidDeterminator.hpp"
#include "njhseq/helpers/seqUtil.hpp"
#include "njhseq/IO/SeqIO.h"

namespace njhseq {

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
	ret["midName_"] = njh::json::toJson(midName_);
	ret["midPos_"] = njh::json::toJson(midPos_);
	ret["barcodeSize_"] = njh::json::toJson(barcodeSize_);
	ret["barcodeScore_"] = njh::json::toJson(barcodeScore_);
	ret["inRevComp_"] = njh::json::toJson(inRevComp_);

	return ret;
}

MidDeterminator::midPos::operator bool() const {
	return npos != midPos_;
}




std::string MidDeterminator::ProcessedRes::getProcessedCaseName(PROCESSED_CASE pCase){
	std::string pCaseStr = "unhandledCase";
	switch (pCase) {
	case MidDeterminator::ProcessedRes::PROCESSED_CASE::MISMATCHING_DIRECTION:
		pCaseStr = "MISMATCHING_DIRECTION";
		break;
	case MidDeterminator::ProcessedRes::PROCESSED_CASE::MISMATCHING_MIDS:
		pCaseStr = "MISMATCHING_MIDS";
		break;
	case MidDeterminator::ProcessedRes::PROCESSED_CASE::NONE:
		pCaseStr = "NONE";
		break;
	case MidDeterminator::ProcessedRes::PROCESSED_CASE::NO_MATCHING:
		pCaseStr = "NO_MATCHING";
		break;
	case MidDeterminator::ProcessedRes::PROCESSED_CASE::TOO_MANY_MATCHING:
		pCaseStr = "TOO_MANY_MATCHING";
		break;
	case MidDeterminator::ProcessedRes::PROCESSED_CASE::PARTIALDUAL:
		pCaseStr = "PARTIALDUAL";
		break;
	case MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH:
		pCaseStr = "MATCH";
		break;
	default:
		pCaseStr = "unhandledCase";
		break;
	}
	return pCaseStr;
}

VecStr MidDeterminator::ProcessedRes::getProcessedCaseNames() {
	return VecStr { "unhandledCase", "MISMATCHING_DIRECTION", "MISMATCHING_MIDS",
			"NONE", "NO_MATCHING", "TOO_MANY_MATCHING", "PARTIALDUAL", "MATCH" };
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
	case MidDeterminator::midPos::FailureCase::PARTIAL:
		fCaseStr = "PARTIAL";
		break;
	default:
		fCaseStr = "unhandledCase";
		break;
	}
	return fCaseStr;
}

VecStr MidDeterminator::midPos::getFailureCaseNames() {
	return VecStr { "NONE", "NO_MATCHING", "TOO_MANY_MATCHING",
			"MISMATCHING_MIDS", "MISMATCHING_DIRECTION", "PARTIAL", "unhandledCase" };
}


double MidDeterminator::midPos::normalizeScoreByLen() const{
	return static_cast<double>(barcodeScore_)/barcodeSize_;
}

MidDeterminator::MidInfo::MidInfo(const std::string & midName,
		const std::string & barcode) :
		midName_(midName),
		bar_(std::make_unique<motif>(barcode)),
		rcompBar_(std::make_unique<motif>(seqUtil::reverseComplement(barcode, "DNA")) ) {
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if(barcode.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, MID: " << midName << " can't be an empty string" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if (barcode.size() > 1) {
		shortenFrontBar_ = std::make_unique<motif>(barcode.substr(1));
		shortenBackBar_ = std::make_unique<motif>(barcode.substr(0, barcode.size() - 1));
		auto rcBar = seqUtil::reverseComplement(barcode, "DNA");
		shortenFrontRCompBar_ = std::make_unique<motif>(rcBar.substr(1));
		shortenBackRCompBar_ = std::make_unique<motif>(rcBar.substr(0, rcBar.size() - 1));
	}
	rCompSame_ = bar_->motifOriginal_ == rcompBar_->motifOriginal_;

}

MidDeterminator::MID::MID(const std::string & name):
	name_(name) {

}

MidDeterminator::MID::MID(const std::string & name, const std::string & forwardBar):
	name_(name),
		forwardBar_(std::make_unique<MidInfo>(name, forwardBar)) {

}

MidDeterminator::MID::MID(const std::string & name,
		const std::string & forwardBar,
		const std::string & reverseBar):
	name_(name),
	reverseBar_(std::make_unique<MidInfo>(name, reverseBar)) {
	if("" != forwardBar){
		forwardBar_ = std::make_unique<MidInfo>(name, forwardBar);
		forSameAsRev_ = forwardBar == reverseBar;
		forSameAsRevShorten_ = forwardBar_->shortenFrontBar_->motifOriginal_ == reverseBar_->shortenFrontBar_->motifOriginal_;
	}
}

bool MidDeterminator::MID::dualBarcoded() const {
	return nullptr != forwardBar_ && nullptr != reverseBar_;
}


MidDeterminator::MidDeterminator(const table & mids, const MidDeterminePars & searchPars):
		searchPars_(searchPars),
		shortenSearchPars_(searchPars) {
	shortenSearchPars_.searchStart_ = 0;
	shortenSearchPars_.searchStop_ = 0;
	mids.checkForColumnsThrow(VecStr { "id", "barcode" }, __PRETTY_FUNCTION__);
	bool containsSecondBarCodeColumn = mids.containsColumn("barcode2");
	if (containsSecondBarCodeColumn) {
		for (const auto & row : mids.content_) {
			addForwardBarcode(
					row[mids.getColPos("id")],
					row[mids.getColPos("barcode")]);
		}
	} else {
		for (const auto & row : mids.content_) {
			if("" == row[mids.getColPos("barcode2")]){
				addForwardBarcode(
						row[mids.getColPos("id")],
						row[mids.getColPos("barcode")]);
			}else{
				addForwardReverseBarcode(
						row[mids.getColPos("id")],
						row[mids.getColPos("barcode")],
						row[mids.getColPos("barcode2")]);
			}
		}
	}
}

MidDeterminator::MidDeterminator(const bfs::path & idFileFnp,
		const MidDeterminePars & searchPars) :
		searchPars_(searchPars), shortenSearchPars_(searchPars) {
	shortenSearchPars_.searchStart_ = 0;
	shortenSearchPars_.searchStop_ = 0;
	if (searchPars.searchStop_ < searchPars.searchStart_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": error, barcode search start should be less than barcode search stop; searchStop: "
				<< searchPars.searchStop_ << ", searchStart: "
				<< searchPars.searchStart_ << "\n";
		throw std::runtime_error { ss.str() };
	}

	if (!bfs::exists(idFileFnp)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, "
				<< njh::bashCT::boldRed(idFileFnp.string()) << " doesn't exist\n";
		throw std::runtime_error { ss.str() };
	}
	auto firstLine = njh::files::getFirstLine(idFileFnp);
	njh::strToLower(firstLine);
	if (!njh::beginsWith(firstLine, "gene")
			&& !njh::beginsWith(firstLine, "target")
			&& !njh::beginsWith(firstLine, "id")) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error the id file "
				<< njh::bashCT::boldRed(idFileFnp.string())
				<< " should start with either target, gene, or id (case insensitive)\n";
		ss << "line: " << firstLine << "\n";
		throw std::runtime_error { ss.str() };
	}
	bool readingPrimers = false;
	bool readingMids = false;
	InputStream idFile { InOptions { idFileFnp } };
	std::string line = "";

	while (njh::files::crossPlatGetline(idFile, line)) {
		if (njh::beginsWith(line, "#") || njh::allWhiteSpaceStr(line)) {
			continue;
		}
		auto lowerLine = njh::strToLowerRet(line);
		auto lowerLineToks = njh::tokenizeString(lowerLine, "whitespace");

		if (lowerLineToks.size() > 1
				&& ("gene" == lowerLineToks[0] || "target" == lowerLineToks[0])) {
			readingPrimers = true;
			readingMids = false;
		} else if (lowerLineToks.size() > 1 && "id" == lowerLineToks[0]) {
			readingPrimers = false;
			readingMids = true;
		} else if (readingPrimers) {
			auto toks = njh::tokenizeString(line, "whitespace");
			if (3 != toks.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error in id file "
						<< njh::bashCT::boldRed(idFileFnp.string())
						<< " primer line should contain 3 items not " << toks.size() << "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error { ss.str() };
			}
			//addTarget(toks[0], toks[1], toks[2]);
		} else if (readingMids) {
			auto toks = njh::tokenizeString(line, "whitespace");
			if (2 == toks.size()) {
				addForwardBarcode(toks[0], toks[1]);
			} else if (3 == toks.size()) {
				addForwardReverseBarcode(toks[0], toks[1], toks[2]);
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": error in id file "
						<< njh::bashCT::boldRed(idFileFnp.string())
						<< " barcode line should contain 2 items not " << toks.size()
						<< "\n";
				ss << "line: " << line << "\n";
				throw std::runtime_error { ss.str() };
			}
		}
	}
}






void MidDeterminator::addForwardReverseBarcode(const std::string & name,
		const std::string & forward, const std::string & reverse) {
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	containsMidByNameThrow(name,__PRETTY_FUNCTION__);
	if("" == forward){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error for MID: " << name << " forward can't be blank" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if("" == reverse){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error for MID: " << name << " reverse can't be blank" << "\n";
		throw std::runtime_error{ss.str()};
	}
	for (const auto & otherMid : mids_) {
		//check forward barcode only mids if they have this forward mid
		if ((nullptr != otherMid.second.forwardBar_ && nullptr == otherMid.second.reverseBar_)
				&& forward == otherMid.second.forwardBar_->bar_->motifOriginal_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error for MID: " << name
					<< " forward barcode, " << forward << " already found for MID: "
					<< otherMid.second.forwardBar_->midName_ << "\n";
			throw std::runtime_error { ss.str() };
		}

		//check other pairings
		if ((nullptr != otherMid.second.forwardBar_
				&& forward == otherMid.second.forwardBar_->bar_->motifOriginal_)
				&& (nullptr != otherMid.second.reverseBar_
						&& reverse == otherMid.second.reverseBar_->bar_->motifOriginal_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error for MID: " << name
					<< " pairing, forward: " << forward << " reverse: " << reverse
					<< " already found with name: " << otherMid.first << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	mids_.emplace(name, MID(name, forward, reverse));
}

void MidDeterminator::addForwardBarcode(const std::string & name,
		const std::string & forward) {
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	containsMidByNameThrow(name, __PRETTY_FUNCTION__);
	if ("" == forward) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error for MID: " << name
				<< " forward can't be blank" << "\n";
		throw std::runtime_error { ss.str() };
	}
	for (const auto & otherMid : mids_) {
		if (nullptr != otherMid.second.forwardBar_
				&& forward == otherMid.second.forwardBar_->bar_->motifOriginal_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error for MID: " << name
					<< " forward barcode, " << forward << " already found for MID: "
					<< otherMid.second.forwardBar_->midName_ << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	mids_.emplace(name, MID(name, forward));
}

void MidDeterminator::addReverseBarcode(const std::string & name,
		const std::string & reverse) {
	containsMidByNameThrow(name,__PRETTY_FUNCTION__);
	if("" == reverse){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error for MID: " << name << " reverse can't be blank" << "\n";
		throw std::runtime_error{ss.str()};
	}
	mids_.emplace(name, MID(name, "", reverse));
}


bool MidDeterminator::containsMidByName(const std::string & name) const {
	return njh::in(name, mids_);
}
void MidDeterminator::containsMidByNameThrow(const std::string & name, const std::string & funcName) const {
	if(containsMidByName(name)){
		std::stringstream ss;
		ss << funcName << ", error already contain MID name: " << name << "\n";
		throw std::runtime_error{ss.str()};
	}
}




bool MidDeterminator::containsMidByBarcode(const std::string & barcode)const{
	for(const auto & mid : mids_){
		if(nullptr != mid.second.forwardBar_ && mid.second.forwardBar_->bar_->motifOriginal_ == barcode){
			return true;
		}
		if(nullptr != mid.second.reverseBar_ && mid.second.reverseBar_->bar_->motifOriginal_ == barcode){
			return true;
		}
	}
	return false;
}

std::string MidDeterminator::getMidName(const std::string & barcode)const{
	for(const auto & mid : mids_){
		if(nullptr != mid.second.forwardBar_ && mid.second.forwardBar_->bar_->motifOriginal_ == barcode){
			return mid.first;
		}
		if(nullptr != mid.second.reverseBar_ && mid.second.reverseBar_->bar_->motifOriginal_ == barcode){
			return mid.first;
		}
	}
	return "no_name_for_barcode:" + barcode;
}



std::vector<MidDeterminator::midPos> MidDeterminator::frontDeterminePosMIDPos(
		const std::string & seq,
		const motif & bar,
		const std::string & midName,
		const MidDeterminator::MidDeterminePars & mPars) {
	std::vector<midPos> ret;
	if (seq.size() > bar.size()
			&& mPars.allowableErrors_ + 1 <  bar.size()
			&& mPars.searchStart_ + bar.size() < seq.size()) {

		uint32_t searchStop  = std::min<uint32_t>(mPars.searchStop_ + bar.size(), seq.size() );
		uint32_t searchStart = mPars.searchStart_;


		if(0 == mPars.searchStop_ ){
			searchStop = mPars.searchStart_ + bar.size();
		}

		auto positions = bar.findPositionsFull(
				seq,
				mPars.allowableErrors_,
				searchStart,
				searchStop);
		if (!positions.empty()) {
			auto currentScore = bar.scoreMotif(
					seq.substr(positions.front(), bar.size()));
			ret.emplace_back(midName, positions.front(), bar.size(), currentScore);
		}
	}
	return ret;
}

std::vector<MidDeterminator::midPos> MidDeterminator::backDeterminePosMIDPos(
		const std::string & seq, const motif & bar, const std::string & midName,
		const MidDeterminator::MidDeterminePars & mPars) {
	std::vector<midPos> ret;
	if (seq.size() > bar.size()
			&& mPars.allowableErrors_ + 1 <  bar.size()
			&& mPars.searchStart_ + bar.size() < seq.size()) {

		uint32_t searchStart = seq.size() - bar.size();
		uint32_t searchStop = seq.size() - mPars.searchStart_;


		if (mPars.searchStop_ > searchStart) {
			searchStart = 0;
		} else {
			searchStart -= mPars.searchStop_;
		}
		if(0 == mPars.searchStop_ ){
			searchStart = seq.size() - mPars.searchStart_ - bar.size();
		}
		auto positions = bar.findPositionsFull(seq,
				mPars.allowableErrors_,
				searchStart, searchStop);
		if (!positions.empty()) {
			auto currentScore = bar.scoreMotif(
					seq.substr(positions.back(), bar.size()));
			ret.emplace_back(midName, positions.back(), bar.size(), currentScore);
		}
	}
	return ret;
}

MidDeterminator::ProcessedRes MidDeterminator::processSearchPairedEndRead(
		 PairedRead & seq, const MidSearchRes & res) const {

	MidDeterminator::ProcessedRes ret;
	if(res.forward_.empty()){
		//currently requiring forward
		if(1 == res.reverse_.size() && mids_.at(res.reverse_.front().midName_).dualBarcoded()){
			//partial match, failed to find forward barcode of a dual barcoded sample
			ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::PARTIALDUAL;
			MetaDataInName meta;
			meta.addMeta("PartialReverseMID", njh::pasteAsStr(res.reverse_.front().midName_, ":", res.reverse_.front().midPos_, ":", (res.reverse_.front().inRevComp_ ? "InRComp": "InFor")));
			seq.seqBase_.name_.append(meta.createMetaName());
			seq.mateSeqBase_.name_.append(meta.createMetaName());
		}else{
			//no match
			ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::NO_MATCHING;
		}
	}else{
		if(res.forward_.size() == 1 && res.reverse_.size() < 2){
			if(res.reverse_.size() == 1){
				if(res.forward_.front().midName_ == res.reverse_.front().midName_){
					if(res.forward_.front().inRevComp_ == res.reverse_.front().inRevComp_){
						//match
						ret.midName_ = res.forward_.front().midName_;
						ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH;
						ret.rcomplement_ = res.forward_.front().inRevComp_;
						if(res.forward_.front().inRevComp_){
							//process
							seq.seqBase_.trimFront(res.reverse_.front().midPos_ + res.reverse_.front().barcodeSize_);
							seq.mateSeqBase_.trimFront(res.forward_.front().midPos_ + res.forward_.front().barcodeSize_);
							//complement
							seq.seqBase_.name_.append("_Comp");
							seq.mateSeqBase_.name_.append("_Comp");
							seqInfo tempMate = seq.seqBase_;
							seq.seqBase_ = seq.mateSeqBase_;
							seq.mateSeqBase_ = tempMate;
						}else{
							//process
							seq.seqBase_.trimFront(res.forward_.front().midPos_ + res.forward_.front().barcodeSize_);
							seq.mateSeqBase_.trimFront(res.reverse_.front().midPos_ + res.reverse_.front().barcodeSize_);
						}
					}else{
						//mismatching directions
						ret.midName_ = res.forward_.front().midName_;
						ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::MISMATCHING_DIRECTION;
						ret.rcomplement_ = res.forward_.front().inRevComp_;
						MetaDataInName meta;
						meta.addMeta("ReverseMID", njh::pasteAsStr(res.reverse_.front().midName_, ":", res.reverse_.front().midPos_, ":", (res.reverse_.front().inRevComp_ ? "InRComp": "InFor")));
						meta.addMeta("ForwardMID", njh::pasteAsStr(res.forward_.front().midName_, ":", res.forward_.front().midPos_, ":", (res.forward_.front().inRevComp_ ? "InRComp": "InFor")));
						seq.seqBase_.name_.append(meta.createMetaName());
						seq.mateSeqBase_.name_.append(meta.createMetaName());
					}
				}else{
					//mismatch
					ret.midName_ = res.forward_.front().midName_;
					ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::MISMATCHING_MIDS;
					ret.rcomplement_ = res.forward_.front().inRevComp_;
					MetaDataInName meta;
					meta.addMeta("ReverseMID", njh::pasteAsStr(res.reverse_.front().midName_, ":", res.reverse_.front().midPos_, ":", (res.reverse_.front().inRevComp_ ? "InRComp": "InFor")));
					meta.addMeta("ForwardMID", njh::pasteAsStr(res.forward_.front().midName_, ":", res.forward_.front().midPos_, ":", (res.forward_.front().inRevComp_ ? "InRComp": "InFor")));
					seq.seqBase_.name_.append(meta.createMetaName());
					seq.mateSeqBase_.name_.append(meta.createMetaName());
				}
			}else{
				if(mids_.at(res.forward_.front().midName_).dualBarcoded()){
					//partial match, failed to find reverse barcode of a dual barcoded sample

					ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::PARTIALDUAL;
					MetaDataInName meta;
					meta.addMeta("PartialForwardMID", njh::pasteAsStr(res.forward_.front().midName_, ":", res.forward_.front().midPos_, ":", (res.forward_.front().inRevComp_ ? "InRComp": "InFor")));
					seq.seqBase_.name_.append(meta.createMetaName());
					seq.mateSeqBase_.name_.append(meta.createMetaName());

				}else{
					//match
					ret.midName_ = res.forward_.front().midName_;
					ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH;
					ret.rcomplement_ = res.forward_.front().inRevComp_;
					if(res.forward_.front().inRevComp_){
						//process
						seq.mateSeqBase_.trimFront(res.forward_.front().midPos_ + res.forward_.front().barcodeSize_);
						//complement
						seq.seqBase_.name_.append("_Comp");
						seq.mateSeqBase_.name_.append("_Comp");
						seqInfo tempMate = seq.seqBase_;
						seq.seqBase_ = seq.mateSeqBase_;
						seq.mateSeqBase_ = tempMate;
					}else{
						//process
						seq.seqBase_.trimFront(res.forward_.front().midPos_ + res.forward_.front().barcodeSize_);
					}
				}
			}
		}else{
			//indeterminate
			ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::TOO_MANY_MATCHING;
			MetaDataInName meta;
			std::string forwardMids = "";
			std::string reverseMids = "";
			for(const auto & forMid : res.forward_){
				if("" != forwardMids){
					forwardMids.append(",");
				}
				forwardMids.append(njh::pasteAsStr(forMid.midName_, ":", forMid.midPos_, ":", (forMid.inRevComp_ ? "InRComp": "InFor")));
			}

			for(const auto & revMid : res.reverse_){
				if("" != forwardMids){
					reverseMids.append(",");
				}
				reverseMids.append(njh::pasteAsStr(revMid.midName_, ":", revMid.midPos_, ":", (revMid.inRevComp_ ? "InRComp": "InFor")));
			}

			if ("" != forwardMids) {
				meta.addMeta("ForwardMIDs", forwardMids);
			}
			if ("" != reverseMids) {
				meta.addMeta("ReverseMIDs", reverseMids);
			}
			seq.seqBase_.name_.append(meta.createMetaName());
			seq.mateSeqBase_.name_.append(meta.createMetaName());
		}
	}

	return ret;
}


MidDeterminator::MidSearchRes MidDeterminator::searchRead(
		const seqInfo & seq) const {
	std::vector<MidDeterminator::midPos> forwardBarMatches;
	std::vector<MidDeterminator::midPos> reverseBarMatches;

	for (const auto & mid : mids_) {
		if (nullptr != mid.second.forwardBar_) {
			auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
					seq.seq_, *(mid.second.forwardBar_->bar_), mid.second.name_,
					searchPars_);
			addOtherVec(forwardBarMatches, midPositions);
		}
		if (nullptr != mid.second.reverseBar_) {
			auto midPositions = MidDeterminator::backDeterminePosMIDPos(
					seq.seq_, *(mid.second.reverseBar_->rcompBar_),
					mid.second.name_, searchPars_);
			addOtherVec(reverseBarMatches, midPositions);
		}
	}

	if (searchPars_.checkForShorten_) {
		for (const auto & mid : mids_) {
			if (nullptr != mid.second.forwardBar_) {
				auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
						seq.seq_, *(mid.second.forwardBar_->shortenFrontBar_),
						mid.second.name_, shortenSearchPars_);
				addOtherVec(forwardBarMatches, midPositions);
			}
			if (nullptr != mid.second.reverseBar_) {
				////std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				auto midPositions = MidDeterminator::backDeterminePosMIDPos(
						seq.seq_, *(mid.second.reverseBar_->shortenBackRCompBar_),
						mid.second.name_, shortenSearchPars_);
				addOtherVec(reverseBarMatches, midPositions);
			}
		}
	}

	if (searchPars_.checkComplement_) {
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

		for (const auto & mid : mids_) {
			if (mid.second.forSameAsRev_) {
				continue;
			}
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

			if (nullptr != mid.second.forwardBar_) {
				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				auto midPositions = MidDeterminator::backDeterminePosMIDPos(
						seq.seq_, *(mid.second.forwardBar_->rcompBar_),
						mid.second.name_, searchPars_);
				for (auto & m : midPositions) {
					m.inRevComp_ = true;
				}
				addOtherVec(forwardBarMatches, midPositions);
			}
			if (nullptr != mid.second.reverseBar_) {
				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
						seq.seq_, *(mid.second.reverseBar_->bar_),
						mid.second.name_, searchPars_);
				for (auto & m : midPositions) {
					m.inRevComp_ = true;
				}
				addOtherVec(reverseBarMatches, midPositions);
			}
		}
		if (searchPars_.checkForShorten_) {
			for (const auto & mid : mids_) {
				if (mid.second.forSameAsRevShorten_) {
					continue;
				}
				if (nullptr != mid.second.forwardBar_) {
					auto midPositions = MidDeterminator::backDeterminePosMIDPos(
							seq.seq_,
							*(mid.second.forwardBar_->shortenBackRCompBar_), mid.second.name_,
							shortenSearchPars_);
					for (auto & m : midPositions) {
						m.inRevComp_ = true;
					}
					addOtherVec(forwardBarMatches, midPositions);
				}
				if (nullptr != mid.second.reverseBar_) {
					////std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

					auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
							seq.seq_, *(mid.second.reverseBar_->shortenFrontBar_),
							mid.second.name_, shortenSearchPars_);
					for (auto & m : midPositions) {
						m.inRevComp_ = true;
					}
					addOtherVec(reverseBarMatches, midPositions);
				}
			}
		}
	}
	//
	//if(0 == forwardBarMatches.size() && 0 ==reverseBarMatches.size()){
	//right now requiring that the forward barcode is found
	MidDeterminator::MidSearchRes res;
	std::vector<MidDeterminator::midPos> bestForwardBarMatches;
	std::vector<MidDeterminator::midPos> bestReverseBarMatches;
	//determine best matches;
	//forward
//	if("[Mixture=Mixture5;Sample=Sample1-FrontBarcode;backEndBluntEndArtifact=false;complement=true;forwardBarcodePosition=0;forwardPrimerPosition=10;frontEndBluntEndArtifact=false;readNumber=1;refName=PfKH02;reversePrimerPosition=432;targetSeqLength=457]" == seq.name_ ){
//		std::cout << njh::bashCT::boldGreen("forwardBarMatches.size(): ")<< forwardBarMatches.size() << std::endl;
//	}

	if (1 == forwardBarMatches.size()) {
		bestForwardBarMatches = forwardBarMatches;
	} else {
		uint32_t bestBarcodeScore = 0;
		for (const auto & forwardMatch : forwardBarMatches) {
//			if("[Mixture=Mixture5;Sample=Sample1-FrontBarcode;backEndBluntEndArtifact=false;complement=true;forwardBarcodePosition=0;forwardPrimerPosition=10;frontEndBluntEndArtifact=false;readNumber=1;refName=PfKH02;reversePrimerPosition=432;targetSeqLength=457]" == seq.name_ ){
//				std::cout << njh::json::writeAsOneLine(forwardMatch.toJson()) << std::endl;
//			}
			if (forwardMatch.barcodeScore_ > bestBarcodeScore) {
				bestForwardBarMatches.clear();
				bestForwardBarMatches.emplace_back(forwardMatch);
				bestBarcodeScore = forwardMatch.barcodeScore_;
			} else if (forwardMatch.barcodeScore_ == bestBarcodeScore) {
				bestForwardBarMatches.emplace_back(forwardMatch);
			}
		}
	}
	//std::cout << njh::bashCT::boldGreen("reverseBarMatches.size(): ")<< reverseBarMatches.size() << std::endl;
	//reverse
	if (1 == reverseBarMatches.size()) {
		bestReverseBarMatches = reverseBarMatches;
	} else {
		uint32_t bestBarcodeScore = 0;
		for (const auto & reverseMath : reverseBarMatches) {
			if (reverseMath.barcodeScore_ > bestBarcodeScore) {
				bestReverseBarMatches.clear();
				bestReverseBarMatches.emplace_back(reverseMath);
				bestBarcodeScore = reverseMath.barcodeScore_;
			} else if (reverseMath.barcodeScore_ == bestBarcodeScore) {
				bestReverseBarMatches.emplace_back(reverseMath);
			}
		}
	}

	//utilize other size if multiple best matches
	if (bestForwardBarMatches.size() > 1 && bestReverseBarMatches.size() == 1) {
		bool foundMatch = false;
		MidDeterminator::midPos frontMatch;
		for (const auto & bestFor : bestForwardBarMatches) {
			if (bestFor.midName_ == bestReverseBarMatches.front().midName_) {
				frontMatch = bestFor;
				foundMatch = true;
				break;
			}
		}
		if (foundMatch) {
			bestForwardBarMatches.clear();
			bestForwardBarMatches.emplace_back(frontMatch);
		}
	} else if (bestForwardBarMatches.size() == 1
			&& bestReverseBarMatches.size() > 1) {
		bool foundMatch = false;
		MidDeterminator::midPos backMatch;
		for (const auto & bestBack : bestReverseBarMatches) {
			if (bestBack.midName_ == bestForwardBarMatches.front().midName_) {
				backMatch = bestBack;
				foundMatch = true;
				break;
			}
		}
		if (foundMatch) {
			bestReverseBarMatches.clear();
			bestReverseBarMatches.emplace_back(backMatch);
		}
	} else if (bestForwardBarMatches.size() > 1 && bestReverseBarMatches.size() > 1) {
		std::set<std::string> forwards;
		std::set<std::string> reverse;
		for(const auto & bestFor : bestForwardBarMatches){
			forwards.emplace(bestFor.midName_);
		}
		for(const auto & bestRev : bestReverseBarMatches){
			reverse.emplace(bestRev.midName_);
		}

		VecStr inBoth;
    std::set_intersection(forwards.begin(), forwards.end(),
    		reverse.begin(), reverse.end(),
                          std::back_inserter(inBoth));
		if (1 == inBoth.size()) {
			std::vector<MidDeterminator::midPos> bestForwardBarMatchesSel;
			std::vector<MidDeterminator::midPos> bestReverseBarMatchesSel;
			for (const auto & bestFor : bestForwardBarMatches) {
				if (inBoth.front() == bestFor.midName_) {
					bestForwardBarMatchesSel.emplace_back(bestFor);
					break;
				}
			}
			for (const auto & bestRev : bestReverseBarMatches) {
				if (inBoth.front() == bestRev.midName_) {
					bestReverseBarMatchesSel.emplace_back(bestRev);
					break;
				}
			}
			bestForwardBarMatches = bestForwardBarMatchesSel;
			bestReverseBarMatches = bestReverseBarMatchesSel;
		}
	}

	res.forward_ = bestForwardBarMatches;
	res.reverse_ = bestReverseBarMatches;

	return res;
}
MidDeterminator::ProcessedRes MidDeterminator::processSearchRead(seqInfo & seq,
		const MidDeterminator::MidSearchRes & res) const {
	MidDeterminator::ProcessedRes ret;
	if (res.forward_.empty()) {
		//currently requiring forward
		if (1 == res.reverse_.size()
				&& mids_.at(res.reverse_.front().midName_).dualBarcoded()) {
			//partial match, failed to find forward barcode of a dual barcoded sample
			ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::PARTIALDUAL;
			MetaDataInName meta;
			meta.addMeta("PartialReverseMID",
					njh::pasteAsStr(res.reverse_.front().midName_, ":",
							res.reverse_.front().midPos_, ":",
							(res.reverse_.front().inRevComp_ ? "InRComp" : "InFor")));
			seq.name_.append(meta.createMetaName());
		} else {
			//no match
			ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::NO_MATCHING;
		}
	} else {
		if (res.forward_.size() == 1 && res.reverse_.size() < 2) {
			if (res.reverse_.size() == 1) {
				if (res.forward_.front().midName_ == res.reverse_.front().midName_) {
					if (res.forward_.front().inRevComp_
							== res.reverse_.front().inRevComp_) {
						//match
						ret.midName_ = res.forward_.front().midName_;
						ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH;
						ret.rcomplement_ = res.forward_.front().inRevComp_;
						if (res.forward_.front().inRevComp_) {
							//process
							seq.trimBack(
									res.forward_.front().midPos_);
							seq.trimFront(
									res.reverse_.front().midPos_
											+ res.reverse_.front().barcodeSize_);

							//complement
							seq.name_.append("_Comp");
							seq.reverseComplementRead(false, true);
						} else {
							//process
							seq.trimBack(
									res.reverse_.front().midPos_);
							seq.trimFront(
									res.forward_.front().midPos_
											+ res.forward_.front().barcodeSize_);
						}
					} else {
						//mismatching directions
						ret.midName_ = res.forward_.front().midName_;
						ret.case_ =
								MidDeterminator::ProcessedRes::PROCESSED_CASE::MISMATCHING_DIRECTION;
						ret.rcomplement_ = res.forward_.front().inRevComp_;
						MetaDataInName meta;
						meta.addMeta("ReverseMID",
								njh::pasteAsStr(res.reverse_.front().midName_, ":",
										res.reverse_.front().midPos_, ":",
										(res.reverse_.front().inRevComp_ ? "InRComp" : "InFor")));
						meta.addMeta("ForwardMID",
								njh::pasteAsStr(res.forward_.front().midName_, ":",
										res.forward_.front().midPos_, ":",
										(res.forward_.front().inRevComp_ ? "InRComp" : "InFor")));
						seq.name_.append(meta.createMetaName());
					}
				} else {
					//mismatch
					ret.midName_ = res.forward_.front().midName_;
					ret.case_ =
							MidDeterminator::ProcessedRes::PROCESSED_CASE::MISMATCHING_MIDS;
					ret.rcomplement_ = res.forward_.front().inRevComp_;
					MetaDataInName meta;
					meta.addMeta("ReverseMID",
							njh::pasteAsStr(res.reverse_.front().midName_, ":",
									res.reverse_.front().midPos_, ":",
									(res.reverse_.front().inRevComp_ ? "InRComp" : "InFor")));
					meta.addMeta("ForwardMID",
							njh::pasteAsStr(res.forward_.front().midName_, ":",
									res.forward_.front().midPos_, ":",
									(res.forward_.front().inRevComp_ ? "InRComp" : "InFor")));
					seq.name_.append(meta.createMetaName());
				}
			} else {
				if (mids_.at(res.forward_.front().midName_).dualBarcoded()) {
					//partial match, failed to find reverse barcode of a dual barcoded sample

					ret.case_ =
							MidDeterminator::ProcessedRes::PROCESSED_CASE::PARTIALDUAL;
					MetaDataInName meta;
					meta.addMeta("PartialForwardMID",
							njh::pasteAsStr(res.forward_.front().midName_, ":",
									res.forward_.front().midPos_, ":",
									(res.forward_.front().inRevComp_ ? "InRComp" : "InFor")));
					seq.name_.append(meta.createMetaName());

				} else {
//					if("[Mixture=Mixture5;Sample=Sample1-FrontBarcode;backEndBluntEndArtifact=false;complement=true;forwardBarcodePosition=0;forwardPrimerPosition=10;frontEndBluntEndArtifact=false;readNumber=1;refName=PfKH02;reversePrimerPosition=432;targetSeqLength=457]" == seq.name_ ){
//						std::cout << res.forward_.front().toJson() << std::endl;
//					}
					//match
					ret.midName_ = res.forward_.front().midName_;
					ret.case_ = MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH;
					ret.rcomplement_ = res.forward_.front().inRevComp_;
					if (res.forward_.front().inRevComp_) {
						//process
						seq.trimBack(
								res.forward_.front().midPos_);
						//complement
						seq.name_.append("_Comp");
						seq.reverseComplementRead(false, true);
					} else {
						//process
						seq.trimFront(
								res.forward_.front().midPos_
										+ res.forward_.front().barcodeSize_);
					}
				}
			}
		} else {
			//indeterminate
			ret.case_ =
					MidDeterminator::ProcessedRes::PROCESSED_CASE::TOO_MANY_MATCHING;
			MetaDataInName meta;
			std::string forwardMids = "";
			std::string reverseMids = "";
			for (const auto & forMid : res.forward_) {
				if ("" != forwardMids) {
					forwardMids.append(",");
				}
				forwardMids.append(
						njh::pasteAsStr(forMid.midName_, ":", forMid.midPos_, ":",
								(forMid.inRevComp_ ? "InRComp" : "InFor")));
			}

			for (const auto & revMid : res.reverse_) {
				if ("" != forwardMids) {
					reverseMids.append(",");
				}
				reverseMids.append(
						njh::pasteAsStr(revMid.midName_, ":", revMid.midPos_, ":",
								(revMid.inRevComp_ ? "InRComp" : "InFor")));
			}

			if ("" != forwardMids) {
				meta.addMeta("ForwardMIDs", forwardMids);
			}
			if ("" != reverseMids) {
				meta.addMeta("ReverseMIDs", reverseMids);
			}
			seq.name_.append(meta.createMetaName());
		}
	}
	return ret;
}

MidDeterminator::MidSearchRes MidDeterminator::searchPairedEndRead(const PairedRead & seq) const{
	std::vector<MidDeterminator::midPos> forwardBarMatches;
	std::vector<MidDeterminator::midPos> reverseBarMatches;

	for (const auto & mid : mids_) {
		if (nullptr != mid.second.forwardBar_) {
			auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
					seq.seqBase_.seq_, *(mid.second.forwardBar_->bar_),
					mid.second.name_, searchPars_);
			addOtherVec(forwardBarMatches, midPositions);
		}
		if (nullptr != mid.second.reverseBar_) {
			auto midPositions =
					MidDeterminator::frontDeterminePosMIDPos(
					seq.mateSeqBase_.seq_, *(mid.second.reverseBar_->bar_),
					mid.second.name_, searchPars_);
			addOtherVec(reverseBarMatches, midPositions);
		}
	}

	if(searchPars_.checkForShorten_){
		for (const auto & mid : mids_) {
			if (nullptr != mid.second.forwardBar_) {
				auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
						seq.seqBase_.seq_, *(mid.second.forwardBar_->shortenFrontBar_),
						mid.second.name_, shortenSearchPars_);
				addOtherVec(forwardBarMatches, midPositions);
			}
			if (nullptr != mid.second.reverseBar_) {
				////std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				auto midPositions =
						MidDeterminator::frontDeterminePosMIDPos(
						seq.mateSeqBase_.seq_, *(mid.second.reverseBar_->shortenFrontBar_),
						mid.second.name_, shortenSearchPars_);
				addOtherVec(reverseBarMatches, midPositions);
			}
		}
	}

	if(searchPars_.checkComplement_){
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

		for (const auto & mid : mids_) {
			if (mid.second.forSameAsRev_){
				continue;
			}
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

			if (nullptr != mid.second.forwardBar_) {
				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
						seq.mateSeqBase_.seq_, *(mid.second.forwardBar_->bar_),
						mid.second.name_, searchPars_);
				for(auto & m : midPositions){
					m.inRevComp_ = true;
				}
				addOtherVec(forwardBarMatches, midPositions);
			}
			if (nullptr != mid.second.reverseBar_) {
				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				auto midPositions =
						MidDeterminator::frontDeterminePosMIDPos(
						seq.seqBase_.seq_, *(mid.second.reverseBar_->bar_),
						mid.second.name_, searchPars_);
				for(auto & m : midPositions){
					m.inRevComp_ = true;
				}
				addOtherVec(reverseBarMatches, midPositions);
			}
		}
		if(searchPars_.checkForShorten_){
			for (const auto & mid : mids_) {
				if (mid.second.forSameAsRevShorten_){
					continue;
				}
				if (nullptr != mid.second.forwardBar_) {
					auto midPositions = MidDeterminator::frontDeterminePosMIDPos(
							seq.mateSeqBase_.seq_, *(mid.second.forwardBar_->shortenFrontBar_),
							mid.second.name_, shortenSearchPars_);
					for(auto & m : midPositions){
						m.inRevComp_ = true;
					}
					addOtherVec(forwardBarMatches, midPositions);
				}
				if (nullptr != mid.second.reverseBar_) {
					////std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

					auto midPositions =
							MidDeterminator::frontDeterminePosMIDPos(
							seq.seqBase_.seq_, *(mid.second.reverseBar_->shortenFrontBar_),
							mid.second.name_, shortenSearchPars_);
					for(auto & m : midPositions){
						m.inRevComp_ = true;
					}
					addOtherVec(reverseBarMatches, midPositions);
				}
			}
		}
	}
	//
	//if(0 == forwardBarMatches.size() && 0 ==reverseBarMatches.size()){
	//right now requiring that the forward barcode is found
	MidDeterminator::MidSearchRes res;
	std::vector<MidDeterminator::midPos> bestForwardBarMatches;
	std::vector<MidDeterminator::midPos> bestReverseBarMatches;
//	if("M04367:213:000000000-CKMCR:1:1102:15808:1331 1:N:0:TAGCTT" == seq.seqBase_.name_){
//		std::cout << "first 8 bases for: " << seq.seqBase_.seq_.substr(0,8) << std::endl;
//		std::cout << "first 8 bases rev: " << seq.mateSeqBase_.seq_.substr(0,8) << std::endl;
//
//		std::cout << "Forward matches: ";
//		for(const auto & f : forwardBarMatches){
//			std::cout << f.toJson() << std::endl << std::endl;
//		}
//		std::cout << "Reverse matches: ";
//		for(const auto & f : reverseBarMatches){
//			std::cout << f.toJson() << std::endl << std::endl;
//		}
//	}
	//determine best matches;

	//forward
	if(1 == forwardBarMatches.size()){
		bestForwardBarMatches = forwardBarMatches;
	}else{
		uint32_t bestBarcodeScore = 0;
		for(const auto & forwardMatch : forwardBarMatches){
			if(forwardMatch.barcodeScore_ > bestBarcodeScore){
				bestForwardBarMatches.clear();
				bestForwardBarMatches.emplace_back(forwardMatch);
				bestBarcodeScore = forwardMatch.barcodeScore_;
			}else if(forwardMatch.barcodeScore_ == bestBarcodeScore){
				bestForwardBarMatches.emplace_back(forwardMatch);
			}
		}
	}
	//std::cout << njh::bashCT::boldGreen("reverseBarMatches.size(): ")<< reverseBarMatches.size() << std::endl;
	//reverse
	if(1 == reverseBarMatches.size()){
		bestReverseBarMatches = reverseBarMatches;
	}else{
		uint32_t bestBarcodeScore = 0;
		for(const auto & reverseMath : reverseBarMatches){
			if(reverseMath.barcodeScore_ > bestBarcodeScore){
				bestReverseBarMatches.clear();
				bestReverseBarMatches.emplace_back(reverseMath);
				bestBarcodeScore = reverseMath.barcodeScore_;
			}else if(reverseMath.barcodeScore_ == bestBarcodeScore){
				bestReverseBarMatches.emplace_back(reverseMath);
			}
		}
	}
//
//	if("M04367:213:000000000-CKMCR:1:1102:15808:1331 1:N:0:TAGCTT" == seq.seqBase_.name_){
//		std::cout << "first 8 bases for: " << seq.seqBase_.seq_.substr(0,8) << std::endl;
//		std::cout << "first 8 bases rev: " << seq.mateSeqBase_.seq_.substr(0,8) << std::endl;
//
//		std::cout << "Best Forward matches: ";
//		for(const auto & f : bestForwardBarMatches){
//			std::cout << f.toJson() << std::endl << std::endl;
//		}
//		std::cout << "Best Reverse matches: ";
//		for(const auto & f : bestReverseBarMatches){
//			std::cout << f.toJson() << std::endl << std::endl;
//		}
//	}


	//utilize other size if multiple best matches
	if(bestForwardBarMatches.size() > 1 && bestReverseBarMatches.size() == 1 ){
		bool foundMatch = false;
		MidDeterminator::midPos frontMatch;
		for(const auto & bestFor : bestForwardBarMatches){
			if(bestFor.midName_ == bestReverseBarMatches.front().midName_){
				frontMatch = bestFor;
				foundMatch = true;
				break;
			}
		}
		if(foundMatch){
			bestForwardBarMatches.clear();
			bestForwardBarMatches.emplace_back(frontMatch);
		}
	}else if(bestForwardBarMatches.size() == 1 && bestReverseBarMatches.size() > 1){
		bool foundMatch = false;
		MidDeterminator::midPos backMatch;
		for(const auto & bestBack : bestReverseBarMatches){
			if(bestBack.midName_ == bestForwardBarMatches.front().midName_){
				backMatch = bestBack;
				foundMatch = true;
				break;
			}
		}
		if(foundMatch){
			bestReverseBarMatches.clear();
			bestReverseBarMatches.emplace_back(backMatch);
		}
	} else if (bestForwardBarMatches.size() > 1 && bestReverseBarMatches.size() > 1) {
		std::set<std::string> forwards;
		std::set<std::string> reverse;
		for(const auto & bestFor : bestForwardBarMatches){
			forwards.emplace(bestFor.midName_);
		}
		for(const auto & bestRev : bestReverseBarMatches){
			reverse.emplace(bestRev.midName_);
		}

		VecStr inBoth;
    std::set_intersection(forwards.begin(), forwards.end(),
    		reverse.begin(), reverse.end(),
                          std::back_inserter(inBoth));
		if (1 == inBoth.size()) {
			std::vector<MidDeterminator::midPos> bestForwardBarMatchesSel;
			std::vector<MidDeterminator::midPos> bestReverseBarMatchesSel;
			for (const auto & bestFor : bestForwardBarMatches) {
				if (inBoth.front() == bestFor.midName_) {
					bestForwardBarMatchesSel.emplace_back(bestFor);
					break;
				}
			}
			for (const auto & bestRev : bestReverseBarMatches) {
				if (inBoth.front() == bestRev.midName_) {
					bestReverseBarMatchesSel.emplace_back(bestRev);
					break;
				}
			}
			bestForwardBarMatches = bestForwardBarMatchesSel;
			bestReverseBarMatches = bestReverseBarMatchesSel;
		}
	}


//	if("M04367:213:000000000-CKMCR:1:1102:15808:1331 1:N:0:TAGCTT" == seq.seqBase_.name_){
//		std::cout << "After post processing" << std::endl;
//		std::cout << "Best Forward matches: ";
//		for(const auto & f : bestForwardBarMatches){
//			std::cout << f.toJson() << std::endl << std::endl;
//		}
//		std::cout << "Best Reverse matches: ";
//		for(const auto & f : bestReverseBarMatches){
//			std::cout << f.toJson() << std::endl << std::endl;
//		}
//		exit(1);
//	}

	res.forward_ = bestForwardBarMatches;
	res.reverse_ = bestReverseBarMatches;

	return res;
}

} /* namespace njhseq */
