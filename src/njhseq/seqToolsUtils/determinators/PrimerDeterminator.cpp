/*
 * PrimerDeterminator.cpp
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
#include "PrimerDeterminator.hpp"
#include "njhseq/helpers/seqUtil.hpp"

namespace njhseq {

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


	if ((!njh::in(std::string("geneName"), primers.columnNames_) && !njh::in(std::string("targetName"), primers.columnNames_))
			|| !njh::in(std::string("forwardPrimer"), primers.columnNames_)
			|| !njh::in(std::string("reversePrimer"), primers.columnNames_)) {
		throw std::runtime_error {
				"Error in creating PrimerDeterminator, need to have at "
						"least the following three columns, geneName or targetName and forwardPrimer, reversePrimer, only have "
						+ njh::conToStr(primers.columnNames_, ",") };
	}

	std::string idCol = "geneName";

	if(njh::in(std::string("targetName"), primers.columnNames_)){
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
	return njh::in(targetName, primers_);
}
PrimerDeterminator::primerInfo::primerInfo():
		forwardPrimerMotif_("NNNN"),
		forwardPrimerMotifRevDir_ ("NNNN"),
		reversePrimerMotif_("NNNN"),
		reversePrimerMotifForDir_ ("NNNN"){

}

PrimerDeterminator::primerInfo::primerInfo(const std::string & primerPairName,
		const std::string & forwardPrimer, const std::string &reversePrimer) :
		primerPairName_(primerPairName),
		forwardPrimer_(forwardPrimer),
		forwardPrimerInfo_(seqInfo { primerPairName, forwardPrimer }),
		forwardPrimerInfoRevDir_(seqInfo { primerPairName, seqUtil::reverseComplement(forwardPrimer,"DNA") }),
		forwardPrimerMotif_(seqUtil::genMotifStrAccountDegenBase(forwardPrimer)),
		forwardPrimerMotifRevDir_ (seqUtil::genMotifStrAccountDegenBase(forwardPrimerInfoRevDir_.seq_)),
		reversePrimer_(reversePrimer),
		reversePrimerInfo_(seqInfo { primerPairName, seqUtil::reverseComplement(reversePrimer,"DNA") }),
		reversePrimerInfoForDir_(seqInfo { primerPairName,reversePrimer } ),
		reversePrimerMotif_(seqUtil::genMotifStrAccountDegenBase(reversePrimer)),
		reversePrimerMotifForDir_ (seqUtil::genMotifStrAccountDegenBase(reversePrimerInfoForDir_.seq_))
		{


	//create motif strings with taking into account ambiguous bases








	forwardPrimerInfoLetCounter_.increaseCountByString(forwardPrimerInfo_.seq_);
	forwardPrimerInfoLetCounter_.resetAlphabet(false);
	forwardPrimerInfoLetCounter_.setFractions();

	forwardPrimerInfoRevDirLetCounter_.increaseCountByString(forwardPrimerInfoRevDir_.seq_);
	forwardPrimerInfoRevDirLetCounter_.resetAlphabet(false);
	forwardPrimerInfoRevDirLetCounter_.setFractions();

	reversePrimerInfoLetCounter_.increaseCountByString(reversePrimerInfo_.seq_);
	reversePrimerInfoLetCounter_.resetAlphabet(false);
	reversePrimerInfoLetCounter_.setFractions();

	reversePrimerInfoForDirLetCounter_.increaseCountByString(reversePrimerInfoForDir_.seq_);
	reversePrimerInfoForDirLetCounter_.resetAlphabet(false);
	reversePrimerInfoForDirLetCounter_.setFractions();

}



std::string PrimerDeterminator::determineWithReversePrimer(seqInfo & info,
		const PrimerDeterminatorPars & pars,
		aligner & alignerObj){
	std::vector<PrimerPositionScore> determinedPrimers;
	const uint32_t maxPSize = getMaxPrimerSize();
	seqInfo readBegin(info.name_ + "_readBegin",
			info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
	charCounter letCounter;
	letCounter.increaseCountByString(readBegin.seq_);

	for (const auto& currentPrimer : primers_) {
		uint32_t basesShared = 2;
		for(const auto c : currentPrimer.second.reversePrimerInfoForDirLetCounter_.alphabet_){
			basesShared += std::min(letCounter.chars_[c], currentPrimer.second.reversePrimerInfoForDirLetCounter_.chars_[c]);
		}
		if(static_cast<double>(basesShared)/currentPrimer.second.reversePrimerInfoForDir_.seq_.size() < pars.allowable_.distances_.query_.coverage_){
			continue;
		}

		// find reverse primer in forward direction or if it isn't found return unrecognized
//		auto readBegin = seqInfo(info.name_ + "_readBegin",
//				info.seq_.substr(pars.primerStart_, pars.primerWithin_ + currentPrimer.second.reversePrimerInfoForDir_.seq_.size() + 5));

		/**@todo put in a check to make sure of semi-global alignment */

		if(pars.useMotif_){
			auto positions = currentPrimer.second.reversePrimerMotif_.findPositionsSubSets(readBegin.seq_,
					pars.allowable_.hqMismatches_ + pars.allowable_.lqMismatches_,
					0, len(readBegin),
					1, currentPrimer.second.reversePrimerMotif_.size());
			if(positions.empty()){
				continue;
			}else{
				auto motifPos = positions.front();
				alignerObj.alignObjectA_ = readBegin;
				if(0 == motifPos){
					alignerObj.alignObjectA_.seqBase_.prepend('-');
				}else{
					--motifPos;
				}
				std::string primerSeq = "";
				if(motifPos > 0){
					primerSeq = std::string(motifPos, '-');
				}
				primerSeq.append(currentPrimer.second.reversePrimerInfoForDir_.seq_);
				primerSeq.append(std::string(len(alignerObj.alignObjectA_) - currentPrimer.second.reversePrimerInfoForDir_.seq_.size() - motifPos, '-'));
				alignerObj.alignObjectB_.seqBase_ = seqInfo(currentPrimer.second.reversePrimerInfoForDir_.name_, primerSeq);
			}
		}else{
			alignerObj.alignCacheGlobal(readBegin, currentPrimer.second.reversePrimerInfoForDir_);
			alignerObj.rearrangeObjs(readBegin, currentPrimer.second.reversePrimerInfoForDir_, false);
		}

		/**@todo put in a check to make sure of semi-global alignment */
		alignerObj.profileAlignment(readBegin,
				currentPrimer.second.reversePrimerInfoForDir_, false, true, false);

		std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-")) + pars.primerStart_,
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-")) + pars.primerStart_);
		double coverage = alignerObj.comp_.distances_.query_.coverage_;
		if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
						|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.reversePrimerInfoForDir_.seq_.size() - 1);
		} else if (0 == pars.primerStart_ &&
				2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
						|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.reversePrimerInfoForDir_.seq_.size() - 2);
		}

		if (forwardPosition.first <= pars.primerWithin_
				&& coverage
						>= pars.allowable_.distances_.query_.coverage_
				&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {

			determinedPrimers.emplace_back(forwardPosition.first,
					forwardPosition.second, currentPrimer.second.reversePrimerInfoForDir_.name_,
					alignerObj.comp_);
		}
	}
	PrimerPositionScore bestPrimer;
	if (0 == determinedPrimers.size()) {
		info.on_ = false;
		return "unrecognized";
	} else if (1 == determinedPrimers.size()) {
		bestPrimer = determinedPrimers.front();
	} else {
		bestPrimer = determinedPrimers.front();
		for(const auto & pos : iter::range<uint32_t>(1, determinedPrimers.size())){
			if(determinedPrimers[pos].getNormalizedScore() > bestPrimer.getNormalizedScore()){
				bestPrimer = determinedPrimers[pos];
			}else if(determinedPrimers[pos].getNormalizedScore() == bestPrimer.getNormalizedScore()){
				if(determinedPrimers[pos].comp_.highQualityMatches_ > bestPrimer.comp_.highQualityMatches_){
					bestPrimer = determinedPrimers[pos];
				}else if(determinedPrimers[pos].comp_.highQualityMatches_ == bestPrimer.comp_.highQualityMatches_){
					if(determinedPrimers[pos].start_ < bestPrimer.start_){
						bestPrimer = determinedPrimers[pos];
					}
				}
			}
		}
	}

	if (pars.primerWithin_ != 0 && bestPrimer.start_ != 0) {
		info.setClip(bestPrimer.start_, info.seq_.size() - 1);
	}
	if (pars.primerToLowerCase_) {
		changeSubStrToLowerFromBegining(info.seq_,
				bestPrimer.end_ - bestPrimer.start_);
	}
	info.on_ = true;
	return bestPrimer.primerName_;
}



PrimerDeterminator::PrimerPositionScore PrimerDeterminator::determineBestReversePrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj){
	std::vector<PrimerPositionScore> determinedPrimers;
	const uint32_t maxPSize = getMaxPrimerSize();
	seqInfo readBegin(info.name_ + "_readBegin",
			info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
	charCounter letCounter;
	letCounter.increaseCountByString(readBegin.seq_);

	for (const auto& currentPrimer : primers_) {
		uint32_t basesShared = 2;
		for(const auto c : currentPrimer.second.reversePrimerInfoForDirLetCounter_.alphabet_){
			basesShared += std::min(letCounter.chars_[c], currentPrimer.second.reversePrimerInfoForDirLetCounter_.chars_[c]);
		}
		if(static_cast<double>(basesShared)/currentPrimer.second.reversePrimerInfoForDir_.seq_.size() < pars.allowable_.distances_.query_.coverage_){
			continue;
		}

		// find reverse primer in forward direction or if it isn't found return unrecognized
//		auto readBegin = seqInfo(info.name_ + "_readBegin",
//				info.seq_.substr(pars.primerStart_, pars.primerWithin_ + currentPrimer.second.reversePrimerInfoForDir_.seq_.size() + 5));
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
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-")) + pars.primerStart_,
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-")) + pars.primerStart_);
		double coverage = alignerObj.comp_.distances_.query_.coverage_;
		if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
						|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.reversePrimerInfoForDir_.seq_.size() - 1);
		} else if (0 == pars.primerStart_ &&
				2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
						|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.reversePrimerInfoForDir_.seq_.size() - 2);
		}

		if (forwardPosition.first <= pars.primerWithin_
				&& coverage >= pars.allowable_.distances_.query_.coverage_
				&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {
			determinedPrimers.emplace_back(forwardPosition.first,
					forwardPosition.second, currentPrimer.second.forwardPrimerInfo_.name_,
					alignerObj.comp_);
		}
	}
	PrimerPositionScore bestPrimer;
	if (0 == determinedPrimers.size()) {
	} else if (1 == determinedPrimers.size()) {
		bestPrimer = determinedPrimers.front();
	} else {
		bestPrimer = determinedPrimers.front();
		for(const auto & pos : iter::range<uint32_t>(1, determinedPrimers.size())){
			if(determinedPrimers[pos].getNormalizedScore() > bestPrimer.getNormalizedScore()){
				bestPrimer = determinedPrimers[pos];
			}else if(determinedPrimers[pos].getNormalizedScore() == bestPrimer.getNormalizedScore()){
				if(determinedPrimers[pos].comp_.highQualityMatches_ > bestPrimer.comp_.highQualityMatches_){
					bestPrimer = determinedPrimers[pos];
				}else if(determinedPrimers[pos].comp_.highQualityMatches_ == bestPrimer.comp_.highQualityMatches_){
					if(determinedPrimers[pos].start_ < bestPrimer.start_){
						bestPrimer = determinedPrimers[pos];
					}
				}
			}
		}
	}
	return bestPrimer;
}

PrimerDeterminator::PrimerPositionScore PrimerDeterminator::determineBestForwardPrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj){
	std::vector<PrimerPositionScore> determinedPrimers;
	const uint32_t maxPSize = getMaxPrimerSize();
	seqInfo readBegin(info.name_ + "_readBegin",
			info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
	charCounter letCounter;
	letCounter.increaseCountByString(readBegin.seq_);

	for (const auto& currentPrimer : primers_) {
		uint32_t basesShared = 2;
		for(const auto c : currentPrimer.second.forwardPrimerInfoLetCounter_.alphabet_){
			basesShared += std::min(letCounter.chars_[c], currentPrimer.second.forwardPrimerInfoLetCounter_.chars_[c]);
		}
		if(static_cast<double>(basesShared)/currentPrimer.second.forwardPrimerInfo_.seq_.size() < pars.allowable_.distances_.query_.coverage_){
			continue;
		}
// find reverse primer in forward direction or if it isn't found return unrecognized
//		auto readBegin = seqInfo(info.name_ + "_readBegin",
//				info.seq_.substr(pars.primerStart_, pars.primerWithin_ + currentPrimer.second.forwardPrimerInfo_.seq_.size() + 5));
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
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-")) + pars.primerStart_,
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-")) + pars.primerStart_);
		double coverage = alignerObj.comp_.distances_.query_.coverage_;
		if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
						|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 1);
		} else if (0 == pars.primerStart_ &&
				2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
						|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 2);
		}

		if (forwardPosition.first <= pars.primerWithin_
				&& coverage >= pars.allowable_.distances_.query_.coverage_
				&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {
			determinedPrimers.emplace_back(forwardPosition.first,
					forwardPosition.second, currentPrimer.second.forwardPrimerInfo_.name_,
					alignerObj.comp_);
		}
	}
	PrimerPositionScore bestPrimer;
	if (0 == determinedPrimers.size()) {
	} else if (1 == determinedPrimers.size()) {
		bestPrimer = determinedPrimers.front();
	} else {
		bestPrimer = determinedPrimers.front();
		for(const auto & pos : iter::range<uint32_t>(1, determinedPrimers.size())){
			if(determinedPrimers[pos].getNormalizedScore() > bestPrimer.getNormalizedScore()){
				bestPrimer = determinedPrimers[pos];
			}else if(determinedPrimers[pos].getNormalizedScore() == bestPrimer.getNormalizedScore()){
				if(determinedPrimers[pos].comp_.highQualityMatches_ > bestPrimer.comp_.highQualityMatches_){
					bestPrimer = determinedPrimers[pos];
				}else if(determinedPrimers[pos].comp_.highQualityMatches_ == bestPrimer.comp_.highQualityMatches_){
					if(determinedPrimers[pos].start_ < bestPrimer.start_){
						bestPrimer = determinedPrimers[pos];
					}
				}
			}
		}
	}
	return bestPrimer;
}

std::string PrimerDeterminator::determineForwardPrimer(seqInfo & info,
		const PrimerDeterminatorPars & pars, aligner & alignerObj) {
	std::vector<PrimerPositionScore> determinedPrimers;
	const uint32_t maxPSize = getMaxPrimerSize();
	seqInfo readBegin(info.name_ + "_readBegin",
			info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
	charCounter letCounter;
	letCounter.increaseCountByString(readBegin.seq_);

	for (const auto& currentPrimer : primers_) {
		uint32_t basesShared = 2;
		for(const auto c : currentPrimer.second.forwardPrimerInfoLetCounter_.alphabet_){
			basesShared += std::min(letCounter.chars_[c], currentPrimer.second.forwardPrimerInfoLetCounter_.chars_[c]);
		}
		if(static_cast<double>(basesShared)/currentPrimer.second.forwardPrimerInfo_.seq_.size() < pars.allowable_.distances_.query_.coverage_){
			continue;
		}
// find reverse primer in forward direction or if it isn't found return unrecognized
//		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
//				currentPrimer.second.forwardPrimerInfo_.seq_);
//		alignerObj.rearrangeObjs(readBegin,
//				currentPrimer.second.forwardPrimerInfo_, true);
//		alignerObj.profilePrimerAlignment(readBegin,
//				currentPrimer.second.forwardPrimerInfo_);

		/**@todo put in a check to make sure of semi-global alignment */
		if(pars.useMotif_){

			auto positions = currentPrimer.second.forwardPrimerMotif_.findPositionsSubSets(readBegin.seq_,
					pars.allowable_.hqMismatches_ + pars.allowable_.lqMismatches_,
					0, len(readBegin),
					1, currentPrimer.second.forwardPrimerMotif_.size());
			if(positions.empty()){
				continue;
			}else{
				auto motifPos = positions.front();
				alignerObj.alignObjectA_ = readBegin;
				if(0 == motifPos){
					alignerObj.alignObjectA_.seqBase_.prepend('-');
				}else{
					--motifPos;
				}
				std::string primerSeq = "";
				if(motifPos > 0){
					primerSeq = std::string(motifPos, '-');
				}
				primerSeq.append(currentPrimer.second.forwardPrimerInfo_.seq_);
				primerSeq.append(std::string(len(alignerObj.alignObjectA_) - currentPrimer.second.forwardPrimerInfo_.seq_.size() - motifPos, '-'));
				alignerObj.alignObjectB_.seqBase_ = seqInfo(currentPrimer.second.forwardPrimerInfo_.name_, primerSeq);
			}
		}else{
			alignerObj.alignCacheGlobal(readBegin, currentPrimer.second.forwardPrimerInfo_);
			alignerObj.rearrangeObjs(readBegin, currentPrimer.second.forwardPrimerInfo_, false);
		}
		alignerObj.profileAlignment(readBegin, currentPrimer.second.forwardPrimerInfo_, false, true, false);

		std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-")) + pars.primerStart_,
				alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-")) + pars.primerStart_);
		double coverage = alignerObj.comp_.distances_.query_.coverage_;
		if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
						|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 1);
		} else if (0 == pars.primerStart_ &&
				2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
				&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
						|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
			coverage =
					static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
							/ (currentPrimer.second.forwardPrimerInfo_.seq_.size() - 2);
		}
		if (forwardPosition.first <= pars.primerWithin_
				&& coverage >= pars.allowable_.distances_.query_.coverage_
				&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {
			determinedPrimers.emplace_back(forwardPosition.first,
					forwardPosition.second, currentPrimer.second.forwardPrimerInfo_.name_,
					alignerObj.comp_);
		}
	}
	PrimerPositionScore bestPrimer;
	if (0 == determinedPrimers.size()) {
		info.on_ = false;
		return "unrecognized";
	} else if (1 == determinedPrimers.size()) {
		bestPrimer = determinedPrimers.front();
	} else {
		bestPrimer = determinedPrimers.front();
		for(const auto & pos : iter::range<uint32_t>(1, determinedPrimers.size())){
			if(determinedPrimers[pos].getNormalizedScore() > bestPrimer.getNormalizedScore()){
				bestPrimer = determinedPrimers[pos];
			}else if(determinedPrimers[pos].getNormalizedScore() == bestPrimer.getNormalizedScore()){
				if(determinedPrimers[pos].comp_.highQualityMatches_ > bestPrimer.comp_.highQualityMatches_){
					bestPrimer = determinedPrimers[pos];
				}else if(determinedPrimers[pos].comp_.highQualityMatches_ == bestPrimer.comp_.highQualityMatches_){
					if(determinedPrimers[pos].start_ < bestPrimer.start_){
						bestPrimer = determinedPrimers[pos];
					}
				}
			}
		}
	}

	if (pars.primerWithin_ != 0 && bestPrimer.start_ != 0) {
		info.setClip(bestPrimer.start_, info.seq_.size() - 1);
	}
	if (pars.primerToLowerCase_) {
		changeSubStrToLowerFromBegining(info.seq_,
				bestPrimer.end_ - bestPrimer.start_);
	}
	info.on_ = true;
	return bestPrimer.primerName_;
}

bool PrimerDeterminator::checkForReversePrimer(seqInfo & info,
		const std::string & primerName,const PrimerDeterminatorPars & pars, aligner & alignObj) {
	if (!njh::in(primerName, primers_)) {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	auto trimBackSize = pars.primerWithin_ + primers_[primerName].reversePrimerInfo_.seq_.size() * 2;
	if (trimBackSize < len(info)) {
		readEnd = info.getSubRead(len(info) - trimBackSize, len(info) - pars.primerStart_);
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
	if (!njh::in(primerName, primers_)) {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	seqInfo readEnd;
	auto trimBackSize = pars.primerWithin_ + primers_[primerName].forwardPrimerInfoRevDir_.seq_.size() * 2;
	if (trimBackSize < len(info)) {
		readEnd = info.getSubRead(len(info) - trimBackSize, len(info) - pars.primerStart_);
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


} /* namespace njhseq */
