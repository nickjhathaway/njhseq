//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
/*
 * PrimerDeterminator.cpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */

#include "PrimerDeterminator.hpp"

namespace bibseq {

PrimerDeterminator::PrimerDeterminator(const table & primers) {
	if (!bib::in(std::string("geneName"), primers.columnNames_)
			|| !bib::in(std::string("forwardPrimer"), primers.columnNames_)
			|| !bib::in(std::string("reversePrimer"), primers.columnNames_)) {
		throw std::runtime_error {
				"Error in creating MidDeterminator, need to have at "
						"least two columns, geneName, forwardPrimer, reversePrimer, only have "
						+ bib::conToStr(primers.columnNames_, ",") };
	}
	for (const auto & row : primers.content_) {
		primers_[row[primers.getColPos("geneName")]] = {row[primers.getColPos("geneName")],
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
		aligner & alignerObj, const comparison & allowable, bool forwardPrimerToLowerCase,
		bool weighHomopolyers){
	for (const auto& currentPrimer : primers_) {
		// find reverse primer inforward direction or if it isn't found return unrecognized
		auto readBegin = seqInfo("readBegin",
				info.seq_.substr(0,
						withinPos + currentPrimer.second.reversePrimerInfoForDir_.seq_.size()));
		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
				currentPrimer.second.reversePrimerInfoForDir_.seq_);
		alignerObj.rearrangeObjs(readBegin, currentPrimer.second.reversePrimerInfoForDir_,
				true);
		alignerObj.profilePrimerAlignment(readBegin,
				currentPrimer.second.reversePrimerInfoForDir_, weighHomopolyers);
		if (forwardPosition.first <= withinPos
				&& alignerObj.comp_.distances_.queryCoverage_
						>= allowable.distances_.queryCoverage_
				&& allowable.passErrorProfile(alignerObj.comp_)) {
			if (forwardPrimerToLowerCase) {
				changeSubStrToLowerFromBegining(info.seq_,
						currentPrimer.second.reversePrimerInfoForDir_.seq_.size());
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
		const comparison & allowable, bool forwardPrimerToLowerCase,
		bool weighHomopolyers) {
	for (const auto& currentPrimer : primers_) {
		// find forward primer or if it isn't found return unrecognized
		auto readBegin = seqInfo("readBegin",
				info.seq_.substr(0,
						withinPos + currentPrimer.second.forwardPrimerInfo_.seq_.size()));
		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
				currentPrimer.second.forwardPrimerInfo_.seq_);
		alignerObj.rearrangeObjs(readBegin, currentPrimer.second.forwardPrimerInfo_,
				true);
		alignerObj.profilePrimerAlignment(readBegin,
				currentPrimer.second.forwardPrimerInfo_, weighHomopolyers);
		if (forwardPosition.first <= withinPos
				&& alignerObj.comp_.distances_.queryCoverage_
						>= allowable.distances_.queryCoverage_
				&& allowable.passErrorProfile(alignerObj.comp_)) {
			if (forwardPrimerToLowerCase) {
				changeSubStrToLowerFromBegining(info.seq_,
						currentPrimer.second.forwardPrimerInfo_.seq_.size());
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
		bool weighHomopolyers) {
	if (!bib::in(primerName, primers_)) {
		throw std::runtime_error { "checkForReversePrimer: No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	if (primers_[primerName].reversePrimer_.size() * 2 < len(info)) {
		readEnd = info.getSubRead(
				info.seq_.size() - primers_[primerName].reversePrimer_.size() * 2);
	} else {
		readEnd = info;
	}
	auto rPos = alignObj.findReversePrimer(readEnd.seq_,
			primers_[primerName].reversePrimerInfo_.seq_);
	if (primers_[primerName].reversePrimer_.size() * 2 < len(info)) {
		rPos.first += info.seq_.size()
				- primers_[primerName].reversePrimer_.size() * 2;
		rPos.second += info.seq_.size()
				- primers_[primerName].reversePrimer_.size() * 2;
	} else {
		readEnd = info;
	}
	alignObj.rearrangeObjs(readEnd, primers_[primerName].reversePrimerInfo_,
			true);
	alignObj.profilePrimerAlignment(readEnd,
			primers_[primerName].reversePrimerInfo_, weighHomopolyers);
	bool primerGood = true;
	if (alignObj.comp_.distances_.queryCoverage_
			< allowable.distances_.queryCoverage_
			|| !allowable.passErrorProfile(alignObj.comp_)) {
		primerGood = false;
	}
	info.on_ = primerGood;
	/*if(!primerGood){
		std::cout << std::endl;
		alignObj.comp_.printErrors(std::cout);
		allowable.printErrors(std::cout);
		alignObj.alignObjectA_.seqBase_.outPutSeq(std::cout);
		alignObj.alignObjectB_.seqBase_.outPutSeq(std::cout);
		primers_[primerName].reversePrimerInfo_.outPutSeq(std::cout);
		std::cout << "queryCoverage_: " << alignObj.comp_.distances_.queryCoverage_ << std::endl;
		std::cout << bib::colorBool(alignObj.comp_.distances_.queryCoverage_< allowable.distances_.queryCoverage_) << std::endl;
		std::cout << bib::colorBool(!allowable.passErrorProfile(alignObj.comp_)) << std::endl;
		std::cout << bib::colorBool(alignObj.comp_.distances_.queryCoverage_
				< allowable.distances_.queryCoverage_
				|| !allowable.passErrorProfile(alignObj.comp_)) << std::endl;
		std::cout << std::endl;
	}*/
	if (primerGood) {
		if (reversePrimerToLowerCase) {
			while (info.seq_[rPos.first] == info.seq_[rPos.first - 1]) {
				--rPos.first;
			}
			changeSubStrToLowerToEnd(info.seq_, rPos.first);
		}
		info.setClip(0, rPos.second);
	}
	return primerGood;
}

bool PrimerDeterminator::checkForForwardPrimerInRev(seqInfo & info, const std::string & primerName,
		aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
    bool weighHomopolyers){
	if (!bib::in(primerName, primers_)) {
		throw std::runtime_error { "checkForReversePrimer: No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	if (primers_[primerName].forwardPrimer_.size() * 2 < len(info)) {
		readEnd = info.getSubRead(
				info.seq_.size() - primers_[primerName].forwardPrimer_.size() * 2);
	} else {
		readEnd = info;
	}
	auto rPos = alignObj.findReversePrimer(readEnd.seq_,
			primers_[primerName].forwardPrimerInfoRevDir_.seq_);
	if (primers_[primerName].forwardPrimer_.size() * 2 < len(info)) {
		rPos.first += info.seq_.size()
				- primers_[primerName].forwardPrimer_.size() * 2;
		rPos.second += info.seq_.size()
				- primers_[primerName].forwardPrimer_.size() * 2;
	} else {
		readEnd = info;
	}
	alignObj.rearrangeObjs(readEnd, primers_[primerName].forwardPrimerInfoRevDir_,
			true);
	alignObj.profilePrimerAlignment(readEnd,
			primers_[primerName].forwardPrimerInfoRevDir_, weighHomopolyers);
	bool primerGood = true;
	if (alignObj.comp_.distances_.queryCoverage_
			< allowable.distances_.queryCoverage_
			|| !allowable.passErrorProfile(alignObj.comp_)) {
		primerGood = false;
	}
	info.on_ = primerGood;

	if (primerGood) {
		if (reversePrimerToLowerCase) {
			while (info.seq_[rPos.first] == info.seq_[rPos.first - 1]) {
				--rPos.first;
			}
			changeSubStrToLowerToEnd(info.seq_, rPos.first);
		}
		info.setClip(0, rPos.second);
	}
	return primerGood;
}


} /* namespace bibseq */
