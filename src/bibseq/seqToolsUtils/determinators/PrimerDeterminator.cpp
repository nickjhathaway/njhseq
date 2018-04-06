/*
 * PrimerDeterminator.cpp
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
		if(containsTarget(row[primers.getColPos(idCol)])){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ <<", error already contain target information for " << row[primers.getColPos(idCol)] << "\n";
			throw std::runtime_error{ss.str()};
		}
		primers_[row[primers.getColPos(idCol)]] = {
			row[primers.getColPos(idCol)],
			row[primers.getColPos("forwardPrimer")],
			row[primers.getColPos("reversePrimer")]};
	}
}

PrimerDeterminator::PrimerDeterminator(const std::unordered_map<std::string, primerInfo> & primers){
	for(const auto & primer : primers){
		primers_.emplace(primer);
	}
}

bool PrimerDeterminator::containsTarget(const std::string & targetName) const{
	return bib::in(targetName, primers_);
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


std::string PrimerDeterminator::determineWithReversePrimer(seqInfo & info, const PrimerDeterminatorPars & pars,
		aligner & alignerObj){
	/**@todo add find best primer*/
	for (const auto& currentPrimer : primers_) {
		// find reverse primer in forward direction or if it isn't found return unrecognized
		auto readBegin = seqInfo(info.name_ + "_readBegin",
				info.seq_.substr(0, pars.primerWithin_ + currentPrimer.second.reversePrimerInfoForDir_.seq_.size() + 5));
//		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
//				currentPrimer.second.reversePrimerInfoForDir_.seq_);
//		alignerObj.rearrangeObjs(readBegin,
//				currentPrimer.second.reversePrimerInfoForDir_, true);
//		alignerObj.profilePrimerAlignment(readBegin,
//				currentPrimer.second.reversePrimerInfoForDir_);

		/**@todo put in a check to make sure of semi-global alignment */
		alignerObj.alignCacheGlobal(readBegin,
				currentPrimer.second.reversePrimerInfoForDir_);
		alignerObj.rearrangeObjs(readBegin,
				currentPrimer.second.reversePrimerInfoForDir_, false);
		alignerObj.profileAlignment(readBegin,
				currentPrimer.second.reversePrimerInfoForDir_, false, true, false);
		std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-")),
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-")));
		double coverage = alignerObj.comp_.distances_.query_.coverage_;
		if (1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
						|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 1);
		} else if (2
				== alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
						|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 2);
		}
		if (forwardPosition.first <= pars.primerWithin_
				&& coverage
						>= pars.allowable_.distances_.query_.coverage_
				&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {
			if(pars.primerWithin_ != 0 && forwardPosition.first != 0){
				info.setClip(forwardPosition.first, info.seq_.size() - 1);
			}
			if (pars.primerToLowerCase_) {
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
		const PrimerDeterminatorPars & pars, aligner & alignerObj) {
	/**@todo add find best primer*/
	for (const auto& currentPrimer : primers_) {
		// find reverse primer in forward direction or if it isn't found return unrecognized
		auto readBegin = seqInfo(info.name_ + "_readBegin",
				info.seq_.substr(0, pars.primerWithin_ + currentPrimer.second.forwardPrimerInfo_.seq_.size() + 5));
//		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
//				currentPrimer.second.forwardPrimerInfo_.seq_);
//		alignerObj.rearrangeObjs(readBegin,
//				currentPrimer.second.forwardPrimerInfo_, true);
//		alignerObj.profilePrimerAlignment(readBegin,
//				currentPrimer.second.forwardPrimerInfo_);

		/**@todo put in a check to make sure of semi-global alignment */
		alignerObj.alignCacheGlobal(readBegin,
				currentPrimer.second.forwardPrimerInfo_);
		alignerObj.rearrangeObjs(readBegin,
				currentPrimer.second.forwardPrimerInfo_, false);
		alignerObj.profileAlignment(readBegin,
				currentPrimer.second.forwardPrimerInfo_, false, true, false);
		std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-")),
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-")));
		double coverage = alignerObj.comp_.distances_.query_.coverage_;
		if (1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
						|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 1);
		} else if (2
				== alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
						|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 2);
		}
		if (forwardPosition.first <= pars.primerWithin_
				&& coverage
						>= pars.allowable_.distances_.query_.coverage_
				&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {
			if(pars.primerWithin_ != 0 && forwardPosition.first != 0){
				info.setClip(forwardPosition.first, info.seq_.size() - 1);
			}
			if (pars.primerToLowerCase_) {
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
		const std::string & primerName,const PrimerDeterminatorPars & pars, aligner & alignObj) {
	if (!bib::in(primerName, primers_)) {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	auto trimBackSize = pars.primerWithin_ + primers_[primerName].reversePrimerInfo_.seq_.size() * 2;
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
			< pars.allowable_.distances_.query_.coverage_
			|| !pars.allowable_.passErrorProfile(alignObj.comp_)) {
		primerGood = false;
	}
	info.on_ = primerGood;

	if (primerGood) {
		if (pars.primerToLowerCase_) {
			if(pars.trimExtra_){
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

bool PrimerDeterminator::checkForForwardPrimerInRev(seqInfo & info, const std::string & primerName,const PrimerDeterminatorPars & pars,
		aligner & alignObj){
	if (!bib::in(primerName, primers_)) {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	auto trimBackSize = pars.primerWithin_ + primers_[primerName].forwardPrimerInfoRevDir_.seq_.size() * 2;
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
			< pars.allowable_.distances_.query_.coverage_
			|| !pars.allowable_.passErrorProfile(alignObj.comp_)) {
		primerGood = false;
	}
	info.on_ = primerGood;

	if (primerGood) {
		if (pars.primerToLowerCase_) {
			if(pars.trimExtra_){
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
