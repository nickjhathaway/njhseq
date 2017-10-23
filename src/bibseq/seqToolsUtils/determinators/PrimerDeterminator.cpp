/*
 * PrimerDeterminator.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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
#include "PrimerDeterminator.hpp"
#include "bibseq/helpers/seqUtil.hpp"

namespace bibseq {




size_t PrimerDeterminator::getMaxPrimerSize() const {
	size_t ret = 0;
	for (const auto & p : primers_) {
		if (len(p.second.forwardPrimerInfo_) > ret) {
			ret = len(p.second.forwardPrimerInfo_);
		}
		if (len(p.second.reversePrimerInfo_) > ret) {
			ret = len(p.second.reversePrimerInfo_);
		}
	}
	return ret;
}

PrimerDeterminator::PrimerDeterminator(const table & primers) {

//	std::cout << bib::colorBool(!bib::in(std::string("geneName"), primers.columnNames_)) << std::endl;
//	std::cout << bib::colorBool(!bib::in(std::string("targetName"), primers.columnNames_)) << std::endl;
//	std::cout << bib::colorBool((!bib::in(std::string("geneName"), primers.columnNames_) && !bib::in(std::string("targetName"), primers.columnNames_))) << std::endl;
//	std::cout << bib::colorBool(!bib::in(std::string("forwardPrimer"), primers.columnNames_)) << std::endl;
//	std::cout << bib::colorBool(!bib::in(std::string("reversePrimer"), primers.columnNames_)) << std::endl;

	if ((!bib::in(std::string("geneName"), primers.columnNames_) && !bib::in(std::string("targetName"), primers.columnNames_))
			|| !bib::in(std::string("forwardPrimer"), primers.columnNames_)
			|| !bib::in(std::string("reversePrimer"), primers.columnNames_)) {
		throw std::runtime_error {
				"Error in creating PrimerDeterminator, need to have at "
						"least the following three columns, geneName or targetName and forwardPrimer, reversePrimer, only have "
						+ bib::conToStr(primers.columnNames_, ",") };
	}

	std::string idCol = "geneName";

	if(bib::in(std::string("targetName"), primers.columnNames_)){
		idCol = "targetName";
	}

	for (const auto & row : primers.content_) {
		primers_[row[primers.getColPos(idCol)]] = {row[primers.getColPos(idCol)],
			row[primers.getColPos("forwardPrimer")], row[primers.getColPos("reversePrimer")]};
	}
}

PrimerDeterminator::primerInfo::primerInfo(const std::string & primerPairName,
		const std::string & forwardPrimer, const std::string &reversePrimer) :
		primerPairName_(primerPairName),
		forwardPrimer_(forwardPrimer),
		forwardPrimerInfo_(seqInfo { primerPairName, forwardPrimer }),
		forwardPrimerInfoRevDir_(seqInfo { primerPairName, seqUtil::reverseComplement(forwardPrimer,"DNA") }),
		reversePrimer_(reversePrimer),
		reversePrimerInfo_(seqInfo { primerPairName, seqUtil::reverseComplement(reversePrimer,"DNA") }),
		reversePrimerInfoForDir_(seqInfo { primerPairName,reversePrimer }) {

}


std::string PrimerDeterminator::determineWithReversePrimer(seqInfo & info, uint32_t withinPos,
		aligner & alignerObj, const comparison & allowable, bool forwardPrimerToLowerCase){
	for (const auto& currentPrimer : primers_) {
		// find reverse primer inforward direction or if it isn't found return unrecognized
		auto readBegin = seqInfo("readBegin",
				info.seq_.substr(0,
						withinPos + currentPrimer.second.reversePrimerInfoForDir_.seq_.size()));
		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
				currentPrimer.second.reversePrimerInfoForDir_.seq_);
		alignerObj.rearrangeObjs(readBegin,
				currentPrimer.second.reversePrimerInfoForDir_, true);
		alignerObj.profilePrimerAlignment(readBegin,
				currentPrimer.second.reversePrimerInfoForDir_);
		if (forwardPosition.first <= withinPos
				&& alignerObj.comp_.distances_.query_.coverage_
						>= allowable.distances_.query_.coverage_
				&& allowable.passErrorProfile(alignerObj.comp_)) {
			if(withinPos != 0 && forwardPosition.first != 0){
				info.setClip(forwardPosition.first, info.seq_.size() - 1);
			}
			if (forwardPrimerToLowerCase) {
				changeSubStrToLowerFromBegining(info.seq_,
						forwardPosition.second - forwardPosition.first);
			}
			info.on_ = true;
			return currentPrimer.second.reversePrimerInfoForDir_.name_;
			break;
		}
	}
	info.on_ = false;
	return "unrecognized";
}



std::string PrimerDeterminator::determineForwardPrimer(seqInfo & info,
		uint32_t withinPos, aligner & alignerObj,
		const comparison & allowable, bool forwardPrimerToLowerCase) {
//	std::cout << "Determining Forward Primers" << std::endl;
//	std::cout << __PRETTY_FUNCTION__ << std::endl;
	for (const auto& currentPrimer : primers_) {
//		std::cout << currentPrimer.first << std::endl;
//		std::cout << currentPrimer.second.forwardPrimer_ << std::endl;
//		std::cout << currentPrimer.second.reversePrimer_ << std::endl;
		// find forward primer or if it isn't found return unrecognized
		auto readBegin = seqInfo("readBegin",
				info.seq_.substr(0,
						withinPos + currentPrimer.second.forwardPrimerInfo_.seq_.size()));
		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
				currentPrimer.second.forwardPrimerInfo_.seq_);
		alignerObj.rearrangeObjs(readBegin, currentPrimer.second.forwardPrimerInfo_,
				true);
		alignerObj.profilePrimerAlignment(readBegin,
				currentPrimer.second.forwardPrimerInfo_);
//		std::cout << "forwardPosition.first: " << forwardPosition.first << std::endl;
//		std::cout << "withinPos: " << withinPos << std::endl;
//		std::cout << "alignerObj.comp_.distances_.query_.coverage_: " << alignerObj.comp_.distances_.query_.coverage_ << std::endl;
//		std::cout << "allowable.distances_.query_.coverage_: " << allowable.distances_.query_.coverage_ << std::endl;
//		std::cout << "allowable.passErrorProfile(alignerObj.comp_): " << bib::colorBool(allowable.passErrorProfile(alignerObj.comp_)) << std::endl;
		if (forwardPosition.first <= withinPos
				&& alignerObj.comp_.distances_.query_.coverage_
						>= allowable.distances_.query_.coverage_
				&& allowable.passErrorProfile(alignerObj.comp_)) {
			if(withinPos != 0 && forwardPosition.first != 0){
				info.setClip(forwardPosition.first,info.seq_.size() - 1);
			}
			if (forwardPrimerToLowerCase) {
				changeSubStrToLowerFromBegining(info.seq_,
						forwardPosition.second - forwardPosition.first);
			}
			info.on_ = true;
			return currentPrimer.second.forwardPrimerInfo_.name_;
			break;
		}
	}
	info.on_ = false;
	return "unrecognized";
}

bool PrimerDeterminator::checkForReversePrimer(seqInfo & info,
		const std::string & primerName, aligner & alignObj,
		const comparison & allowable, bool reversePrimerToLowerCase,
		uint32_t within, bool trimExtra) {
	if (!bib::in(primerName, primers_)) {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	auto trimBackSize = within + primers_[primerName].reversePrimerInfo_.seq_.size() * 2;
	if (trimBackSize < len(info)) {
		readEnd = info.getSubRead(len(info) - trimBackSize);
	} else {
		readEnd = info;
	}
	auto rPos = alignObj.findReversePrimer(readEnd.seq_,
			primers_[primerName].reversePrimerInfo_.seq_);
	if (trimBackSize < len(info)) {
		rPos.first += len(info) - trimBackSize;
		rPos.second += len(info) - trimBackSize;
	}
	alignObj.rearrangeObjs(readEnd, primers_[primerName].reversePrimerInfo_,
			true);
	alignObj.profilePrimerAlignment(readEnd,
			primers_[primerName].reversePrimerInfo_);
	bool primerGood = true;
	//std::cout << std::endl;
	//std::cout << alignObj.comp_.distances_.query_.coverage_ << std::endl;
	//std::cout << alignObj.comp_.hqMismatches_ << std::endl;
	//alignObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
	//alignObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
	//primers_[primerName].reversePrimerInfo_.outPutSeqAnsi(std::cout);
	//std::cout << std::endl;
	if (alignObj.comp_.distances_.query_.coverage_
			< allowable.distances_.query_.coverage_
			|| !allowable.passErrorProfile(alignObj.comp_)) {
		primerGood = false;
	}
	info.on_ = primerGood;

	if (primerGood) {
		if (reversePrimerToLowerCase) {
			if(trimExtra){
				while (rPos.first != 0 && info.seq_[rPos.first] == info.seq_[rPos.first - 1]) {
					--rPos.first;
				}
			}
			changeSubStrToLowerToEnd(info.seq_, rPos.first);
		}
		info.setClip(0, rPos.second);
	}
	return primerGood;
}

bool PrimerDeterminator::checkForForwardPrimerInRev(seqInfo & info, const std::string & primerName,
		aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
    uint32_t within, bool trimExtra){
	if (!bib::in(primerName, primers_)) {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	auto trimBackSize = within + primers_[primerName].forwardPrimerInfoRevDir_.seq_.size() * 2;
	if (trimBackSize < len(info)) {
		readEnd = info.getSubRead(len(info) - trimBackSize);
	} else {
		readEnd = info;
	}

	auto rPos = alignObj.findReversePrimer(readEnd.seq_,
			primers_[primerName].forwardPrimerInfoRevDir_.seq_);
	if (trimBackSize < len(info)) {
		rPos.first += len(info) - trimBackSize;
		rPos.second += len(info) - trimBackSize;
	}
	alignObj.rearrangeObjs(readEnd, primers_[primerName].forwardPrimerInfoRevDir_,
			true);
	alignObj.profilePrimerAlignment(readEnd,
			primers_[primerName].forwardPrimerInfoRevDir_);
	bool primerGood = true;
	if (alignObj.comp_.distances_.query_.coverage_
			< allowable.distances_.query_.coverage_
			|| !allowable.passErrorProfile(alignObj.comp_)) {
		primerGood = false;
	}
	info.on_ = primerGood;

	if (primerGood) {
		if (reversePrimerToLowerCase) {
			if(trimExtra){
				while (rPos.first != 0 && info.seq_[rPos.first] == info.seq_[rPos.first - 1]) {
					--rPos.first;
				}
			}
			changeSubStrToLowerToEnd(info.seq_, rPos.first);
		}
		info.setClip(0, rPos.second);
	}
	return primerGood;
}


} /* namespace bibseq */
