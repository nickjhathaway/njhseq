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
#include "njhseq/seqToolsUtils/seqToolsUtils.hpp"


namespace njhseq {

size_t PrimerDeterminator::getMaxPrimerSize() const {
	size_t ret = 0;
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for (const auto & p : primers_) {
		for(const auto & fwd : p.second.fwds_){
			if (len(fwd.info_) > ret) {
				ret = len(fwd.info_);
			}
		}
		for(const auto & rev : p.second.revs_){
			if (len(rev.info_) > ret) {
				ret = len(rev.info_);
			}
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return ret;
}

size_t PrimerDeterminator::getMinPrimerSize() const {
	size_t ret = std::numeric_limits<size_t>::max();
	for (const auto & p : primers_) {
		for(const auto & fwd : p.second.fwds_){
			if (len(fwd.info_) < ret) {
				ret = len(fwd.info_);
			}
		}
		for(const auto & rev : p.second.revs_){
			if (len(rev.info_) < ret) {
				ret = len(rev.info_);
			}
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
		auto fwdToks = tokenizeString(row[primers.getColPos("forwardPrimer")], ",");
		auto revToks = tokenizeString(row[primers.getColPos("reversePrimer")], ",");
		//check for duplicate primers
		for(const auto & primerSeq : fwdToks){
			for(const auto & p : primers_){
				for(const auto & fwd : p.second.fwds_){
					if (fwd.info_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for forward primer of" << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
					if (fwd.infoRC_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for reverse complement of forward primer of " << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
				for(const auto & rev : p.second.revs_){
					if (rev.info_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for reverse primer of " << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
					if (rev.infoRC_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for reverse complement of reveerse primer of " << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
			}
		}
		for(const auto & primerSeq : revToks){
			for(const auto & p : primers_){
				for(const auto & fwd : p.second.fwds_){
					if (fwd.info_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for forward primer of" << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
					if (fwd.infoRC_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for reverse complement of forward primer of " << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
				for(const auto & rev : p.second.revs_){
					if (rev.info_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for reverse primer of " << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
					if (rev.infoRC_.seq_ == primerSeq) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "adding primers for "<< row[primers.getColPos(idCol)] <<", already have " << primerSeq << "for reverse complement of reveerse primer of " << p.first << "\n";
						throw std::runtime_error{ss.str()};
					}
				}
			}
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

bool PrimerDeterminator::containsPrimerSeq(const std::string & primerSeq) const{
	for(const auto & p : primers_){
		for(const auto & fwd : p.second.fwds_){
			if(fwd.info_.seq_ == primerSeq || fwd.infoRC_.seq_ == primerSeq){
				return true;
			}
		}
		for(const auto & rev : p.second.revs_){
			if(rev.info_.seq_ == primerSeq || rev.infoRC_.seq_ == primerSeq){
				return true;
			}
		}
	}


	return false;
}



PrimerDeterminator::primerInfo::primerInfo(): primerInfo("None", "NNNN", "NNNN"){

}



PrimerDeterminator::primerInfo::PrimerSeq::PrimerSeq(const seqInfo & primer):
		primer_(primer.seq_),
		info_(seqInfo{primer.name_, primer.seq_}),
		infoRC_(seqInfo{primer.name_, seqUtil::reverseComplement(primer.seq_, "DNA")}),
		mot_(seqUtil::genMotifStrAccountDegenBase(primer.seq_)),
		motRC_ (seqUtil::genMotifStrAccountDegenBase(infoRC_.seq_))
		{

	auto degenSeqs = createDegenStrs(info_.seq_);
	for(const auto & degenSeq : degenSeqs){
		infoLetCounter_.increaseCountByString(degenSeq);
	}
	//infoLetCounter_.increaseCountByString(info_.seq_);
	infoLetCounter_.resetAlphabet(false);
	infoLetCounter_.setFractions();

	auto degenSeqsRC = createDegenStrs(infoRC_.seq_);
	for(const auto & degenSeq : degenSeqsRC){
		infoLetCounterRC_.increaseCountByString(degenSeq);
	}
	//infoLetCounterRC_.increaseCountByString(infoRC_.seq_);
	infoLetCounterRC_.resetAlphabet(false);
	infoLetCounterRC_.setFractions();
}


PrimerDeterminator::primerInfo::primerInfo(const std::string & primerPairName,
		const std::string & forwardPrimer, const std::string &reversePrimer) :
		primerPairName_(primerPairName), forwardPrimerRaw_(forwardPrimer), reversePrimerRaw_(reversePrimer)
//		forwardPrimer_(forwardPrimer),
//		forwardPrimerInfo_(seqInfo { primerPairName, forwardPrimer }),
//		forwardPrimerInfoRevDir_(seqInfo { primerPairName, seqUtil::reverseComplement(forwardPrimer,"DNA") }),
//		forwardPrimerMotif_(seqUtil::genMotifStrAccountDegenBase(forwardPrimer)),
//		forwardPrimerMotifRevDir_ (seqUtil::genMotifStrAccountDegenBase(forwardPrimerInfoRevDir_.seq_)),
//		reversePrimer_(reversePrimer),
//		reversePrimerInfo_(seqInfo { primerPairName, seqUtil::reverseComplement(reversePrimer,"DNA") }),
//		reversePrimerInfoForDir_(seqInfo { primerPairName,reversePrimer } ),
//		reversePrimerMotif_(seqUtil::genMotifStrAccountDegenBase(reversePrimer)),
//		reversePrimerMotifForDir_ (seqUtil::genMotifStrAccountDegenBase(reversePrimerInfoForDir_.seq_))
		{


	auto fwdToks = tokenizeString(forwardPrimer, ",");
	for(const auto & fwd : fwdToks){
		fwds_.emplace_back(seqInfo{primerPairName, fwd});
	}
	auto revToks = tokenizeString(reversePrimer, ",");
	for(const auto & rev : revToks){
		revs_.emplace_back(seqInfo{primerPairName, rev});
	}

}
std::string PrimerDeterminator::determineWithReversePrimer(seqInfo & info,
                                                           const PrimerDeterminatorPars & pars,
                                                           aligner & alignerObj) const{
  VecStr primers = getVectorOfMapKeys(primers_);
  return determineWithReversePrimer(info, pars, alignerObj, primers);
}


std::string PrimerDeterminator::determineWithReversePrimer(seqInfo & info,
		const PrimerDeterminatorPars & pars,
		aligner & alignerObj, const VecStr & primers) const{
	std::vector<PrimerPositionScore> determinedPrimers;
	const uint32_t maxPSize = getMaxPrimerSize();
	seqInfo readBegin(info.name_ + "_readBegin",
			info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
	charCounter letCounter;
	letCounter.increaseCountByString(readBegin.seq_);

	for (const auto& currentPrimer : primers_) {
    if(!njh::in(currentPrimer.first, primers)){
      continue;
    }
		for(const auto & rev : currentPrimer.second.revs_){


			uint32_t basesShared = 2;
			for(const auto c : rev.infoLetCounter_.alphabet_){
				basesShared += std::min(letCounter.chars_[c], rev.infoLetCounter_.chars_[c]);
			}
			if(static_cast<double>(basesShared)/rev.info_.seq_.size() < pars.allowable_.distances_.query_.coverage_){
				continue;
			}

			// find reverse primer in forward direction or if it isn't found return unrecognized
	//		auto readBegin = seqInfo(info.name_ + "_readBegin",
	//				info.seq_.substr(pars.primerStart_, pars.primerWithin_ + rev.info_.seq_.size() + 5));

			/**@todo put in a check to make sure of semi-global alignment */

			if(pars.useMotif_){
				auto positions = rev.mot_.findPositionsSubSetsBest(readBegin.seq_,
						pars.allowable_.hqMismatches_ + pars.allowable_.lqMismatches_,
						0, len(readBegin),
						1, rev.mot_.size());
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
					std::string primerSeq;
					if(motifPos > 0){
						primerSeq = std::string(motifPos, '-');
					}
					primerSeq.append(rev.info_.seq_);
					primerSeq.append(std::string(len(alignerObj.alignObjectA_) - rev.info_.seq_.size() - motifPos, '-'));
					alignerObj.alignObjectB_.seqBase_ = seqInfo(rev.info_.name_, primerSeq);
				}
			}else{
				alignerObj.alignCacheGlobal(readBegin, rev.info_);
				alignerObj.rearrangeObjs(readBegin, rev.info_, false);

			}
			/**@todo put in a check to make sure of semi-global alignment */
			alignerObj.profileAlignment(readBegin,
					rev.info_, false, true, false);

			std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
					alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-')) + pars.primerStart_,
					alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-')) + pars.primerStart_);
			double coverage = alignerObj.comp_.distances_.query_.coverage_;
			if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
					&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
							|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
				coverage =
						static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
								/ (rev.info_.seq_.size() - 1);
			} else if (0 == pars.primerStart_ &&
					2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
					&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
							|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
				coverage =
						static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
								/ (rev.info_.seq_.size() - 2);
			}

			if (forwardPosition.first <= pars.primerWithin_
					&& coverage
							>= pars.allowable_.distances_.query_.coverage_
					&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {

				determinedPrimers.emplace_back(forwardPosition.first,
						forwardPosition.second, rev.info_.name_,
						alignerObj.comp_);
			}
		}
	}
	PrimerPositionScore bestPrimer;
	if (determinedPrimers.empty()) {
		info.on_ = false;
		return "unrecognized";
	} else if (1 == determinedPrimers.size()) {
		bestPrimer = determinedPrimers.front();
	} else {
		bestPrimer = determinedPrimers.front();
		for(const auto pos : iter::range<uint32_t>(1, determinedPrimers.size())){
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




PrimerDeterminator::PrimerPositionScore PrimerDeterminator::determineBestForwardPrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj, const VecStr & primers) const {
		std::vector<PrimerPositionScore> determinedPrimers;
	const uint32_t maxPSize = getMaxPrimerSize();
	seqInfo readBegin(info.name_ + "_readBegin",
			info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
	charCounter letCounter;
	letCounter.increaseCountByString(readBegin.seq_);

	for (const auto& currentPrimer : primers_) {
		if(njh::notIn(currentPrimer.first, primers)) {
			continue;
		}
		for(const auto & fwd : currentPrimer.second.fwds_){
			uint32_t basesShared = 2;
			for(const auto c : fwd.infoLetCounter_.alphabet_){
				basesShared += static_cast<uint32_t>(std::min(letCounter.chars_[c], fwd.infoLetCounter_.chars_[c]));
			}
//			std::cout << "fwd: " << fwd.info_.seq_ << std::endl;
//			std::cout << "\t" << static_cast<double>(basesShared) << std::endl;
//			std::cout << "\t" << static_cast<double>(basesShared)/fwd.info_.seq_.size() << std::endl;
			if(static_cast<double>(basesShared)/static_cast<double>(fwd.info_.seq_.size()) < pars.allowable_.distances_.query_.coverage_){
				continue;
			}
	// find reverse primer in forward direction or if it isn't found return unrecognized
	//		auto readBegin = seqInfo(info.name_ + "_readBegin",
	//				info.seq_.substr(pars.primerStart_, pars.primerWithin_ + fwd.info_.seq_.size() + 5));
	//		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
	//				fwd.info_.seq_);
	//		alignerObj.rearrangeObjs(readBegin,
	//				fwd.info_, true);
	//		alignerObj.profilePrimerAlignment(readBegin,
	//				fwd.info_);

			/**@todo put in a check to make sure of semi-global alignment */
			alignerObj.alignCacheGlobal(readBegin,
					fwd.info_);
			alignerObj.rearrangeObjs(readBegin,
					fwd.info_, false);
			alignerObj.profileAlignment(readBegin,
					fwd.info_, false, true, false);

			std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
					alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-')) + pars.primerStart_,
					alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-')) + pars.primerStart_);
			double coverage = alignerObj.comp_.distances_.query_.coverage_;
			if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
					&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
							|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
				coverage =
						static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
								/ (fwd.info_.seq_.size() - 1);
			} else if (0 == pars.primerStart_ &&
					2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
					&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
							|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
				coverage =
						static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
								/ (fwd.info_.seq_.size() - 2);
			}

			if (forwardPosition.first <= pars.primerWithin_
					&& coverage >= pars.allowable_.distances_.query_.coverage_
					&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {
				determinedPrimers.emplace_back(forwardPosition.first,
						forwardPosition.second, fwd.info_.name_,
						alignerObj.comp_);
			}
		}
	}

	PrimerPositionScore bestPrimer;
	if (determinedPrimers.empty()) {
	} else if (1 == determinedPrimers.size()) {
		bestPrimer = determinedPrimers.front();
	} else {
		bestPrimer = determinedPrimers.front();
		for(const auto pos : iter::range<uint32_t>(1, determinedPrimers.size())){
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


PrimerDeterminator::PrimerPositionScore PrimerDeterminator::determineBestForwardPrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj) const{
	VecStr primers = getVectorOfMapKeys(primers_);
 return determineBestForwardPrimerPosFront(info, pars, alignerObj, primers);
}

std::string PrimerDeterminator::determineForwardPrimer(seqInfo & info,
                                                       const PrimerDeterminatorPars & pars, aligner & alignerObj, const VecStr & primers) const{
  std::vector<PrimerPositionScore> determinedPrimers;
  const uint32_t maxPSize = getMaxPrimerSize();
  seqInfo readBegin(info.name_ + "_readBegin",
                    info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
  charCounter letCounter;
  letCounter.increaseCountByString(readBegin.seq_);

  for (const auto& currentPrimer : primers_) {
    if(!njh::in(currentPrimer.first, primers)){
      continue;
    }
    for(const auto & fwd : currentPrimer.second.fwds_){

      uint32_t basesShared = 2;
      for(const auto c : fwd.infoLetCounter_.alphabet_){
        basesShared += std::min(letCounter.chars_[c], fwd.infoLetCounter_.chars_[c]);
      }
      ///@todo determine a faster way to limit which primers are searched, limit which ones are searched first by having a higher cut off and then progressively lower it if no match found
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "static_cast<double>(basesShared)/fwd.info_.seq_.size() < pars.allowable_.distances_.query_.coverage_: " << njh::colorBool(static_cast<double>(basesShared)/fwd.info_.seq_.size() < pars.allowable_.distances_.query_.coverage_) << std::endl;
//			std::cout << "pars.useMotif_: " << njh::colorBool(pars.useMotif_) << std::endl;
//			std::cout << "\tpars.allowable_.distances_.query_.coverage_: " << pars.allowable_.distances_.query_.coverage_ << std::endl;
//			std::cout << "\tbasesShared: " << basesShared << std::endl;
//			std::cout << "\tstatic_cast<double>(basesShared)/fwd.info_.seq_.size(): " << static_cast<double>(basesShared)/fwd.info_.seq_.size() << std::endl;

      if(static_cast<double>(basesShared)/fwd.info_.seq_.size() < pars.allowable_.distances_.query_.coverage_){
        continue;
      }
      // find reverse primer in forward direction or if it isn't found return unrecognized
      //		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
      //				fwd.info_.seq_);
      //		alignerObj.rearrangeObjs(readBegin,
      //				fwd.info_, true);
      //		alignerObj.profilePrimerAlignment(readBegin,
      //				fwd.info_);

      /**@todo put in a check to make sure of semi-global alignment */
      if(pars.useMotif_){
        auto positions = fwd.mot_.findPositionsSubSetsBest(readBegin.seq_,
                                                           pars.allowable_.hqMismatches_ + pars.allowable_.lqMismatches_,
                                                           0, len(readBegin),
                                                           1, fwd.mot_.size());
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
          std::string primerSeq;
          if(motifPos > 0){
            primerSeq = std::string(motifPos, '-');
          }
          primerSeq.append(fwd.info_.seq_);
          primerSeq.append(std::string(len(alignerObj.alignObjectA_) - fwd.info_.seq_.size() - motifPos, '-'));
          alignerObj.alignObjectB_.seqBase_ = seqInfo(fwd.info_.name_, primerSeq);
        }
      }else{
        alignerObj.alignCacheGlobal(readBegin, fwd.info_);
        alignerObj.rearrangeObjs(readBegin, fwd.info_, false);

//        alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//        alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);

      }
      alignerObj.profileAlignment(readBegin, fwd.info_, false, true, false);

      std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
          alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-')) + pars.primerStart_,
          alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-')) + pars.primerStart_);
      double coverage = alignerObj.comp_.distances_.query_.coverage_;
      if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
          && ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
              || 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
        coverage =
            static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
            / (fwd.info_.seq_.size() - 1);
      } else if (0 == pars.primerStart_ &&
                 2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
                 && ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
                     || "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
        coverage =
            static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
            / (fwd.info_.seq_.size() - 2);
      }
      if (forwardPosition.first <= pars.primerWithin_
          && coverage >= pars.allowable_.distances_.query_.coverage_
          && pars.allowable_.passErrorProfile(alignerObj.comp_)) {
        determinedPrimers.emplace_back(forwardPosition.first,
                                       forwardPosition.second, fwd.info_.name_,
                                       alignerObj.comp_);
      }
    }
  }
  PrimerPositionScore bestPrimer;
  if (determinedPrimers.empty()) {
    info.on_ = false;
    return "unrecognized";
  } else if (1 == determinedPrimers.size()) {
    bestPrimer = determinedPrimers.front();
  } else {
    bestPrimer = determinedPrimers.front();
    for(const auto pos : iter::range<uint32_t>(1, determinedPrimers.size())){
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

std::string PrimerDeterminator::determineForwardPrimer(seqInfo & info,
		const PrimerDeterminatorPars & pars, aligner & alignerObj) const{
  VecStr primers = getVectorOfMapKeys(primers_);
  return determineForwardPrimer(info, pars, alignerObj, primers);
}


PrimerDeterminator::PrimerPositionScore PrimerDeterminator::determineBestReversePrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj) const {
	VecStr primers = getVectorOfMapKeys(primers_);
	return determineBestReversePrimerPosFront(info, pars, alignerObj, primers);
}

PrimerDeterminator::PrimerPositionScore PrimerDeterminator::determineBestReversePrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj, const VecStr & primers) const {

	std::vector<PrimerPositionScore> determinedPrimers;
	const uint32_t maxPSize = getMaxPrimerSize();
	seqInfo readBegin(info.name_ + "_readBegin",
			info.seq_.substr(pars.primerStart_, pars.primerWithin_ + maxPSize + 5));
	charCounter letCounter;
	letCounter.increaseCountByString(readBegin.seq_);

	for (const auto& currentPrimer : primers_) {
		if(njh::notIn(currentPrimer.first, primers)) {
			continue;;
		}
		for(const auto & rev : currentPrimer.second.revs_){
			uint32_t basesShared = 2;
			for(const auto c : rev.infoLetCounter_.alphabet_){
				basesShared += std::min(letCounter.chars_[c], rev.infoLetCounter_.chars_[c]);
			}
			if(static_cast<double>(basesShared)/rev.info_.seq_.size() < pars.allowable_.distances_.query_.coverage_){
				continue;
			}

			// find reverse primer in forward direction or if it isn't found return unrecognized
	//		auto readBegin = seqInfo(info.name_ + "_readBegin",
	//				info.seq_.substr(pars.primerStart_, pars.primerWithin_ + rev.info_.seq_.size() + 5));
	//		auto forwardPosition = alignerObj.findReversePrimer(readBegin.seq_,
	//				rev.info_.seq_);
	//		alignerObj.rearrangeObjs(readBegin,
	//				rev.info_, true);
	//		alignerObj.profilePrimerAlignment(readBegin,
	//				rev.info_);

			/**@todo put in a check to make sure of semi-global alignment */
			alignerObj.alignCacheGlobal(readBegin,
					rev.info_);
			alignerObj.rearrangeObjs(readBegin,
					rev.info_, false);
			alignerObj.profileAlignment(readBegin,
					rev.info_, false, true, false);

			std::pair<uint32_t, uint32_t> forwardPosition = std::make_pair(
					alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-')) + pars.primerStart_,
					alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-')) + pars.primerStart_);
			double coverage = alignerObj.comp_.distances_.query_.coverage_;
			if (0 == pars.primerStart_ && 1 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
					&& ('T' == alignerObj.alignObjectB_.seqBase_.seq_.front()
							|| 'A' == alignerObj.alignObjectB_.seqBase_.seq_.front())) {
				coverage =
						static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
								/ (rev.info_.seq_.size() - 1);
			} else if (0 == pars.primerStart_ &&
					2 == alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-')
					&& ("TT" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2)
							|| "AA" == alignerObj.alignObjectB_.seqBase_.seq_.substr(0,2))) {
				coverage =
						static_cast<double>(alignerObj.comp_.distances_.query_.covered_)
								/ (rev.info_.seq_.size() - 2);
			}

			if (forwardPosition.first <= pars.primerWithin_
					&& coverage >= pars.allowable_.distances_.query_.coverage_
					&& pars.allowable_.passErrorProfile(alignerObj.comp_)) {
				determinedPrimers.emplace_back(forwardPosition.first,
						forwardPosition.second, currentPrimer.second.primerPairName_,
						alignerObj.comp_);
			}

		}
	}
	PrimerPositionScore bestPrimer;
	if (determinedPrimers.empty()) {
	} else if (1 == determinedPrimers.size()) {
		bestPrimer = determinedPrimers.front();
	} else {
		bestPrimer = determinedPrimers.front();
		for(const auto pos : iter::range<uint32_t>(1, determinedPrimers.size())){
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


bool PrimerDeterminator::checkForReversePrimer(seqInfo & info,
		const std::string & primerName,const PrimerDeterminatorPars & pars, aligner & alignObj) const{
	if (!njh::in(primerName, primers_)) {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	std::vector<PrimerPositionScore> determinedPrimers;
	for(const auto & rev : primers_.at(primerName).revs_){
		seqInfo readEnd;
		auto trimBackSize = pars.primerWithin_ + rev.infoRC_.seq_.size() * 2;
		if (trimBackSize < len(info)) {
			readEnd = info.getSubRead(len(info) - trimBackSize, len(info) - pars.primerStart_);
		} else {
			readEnd = info;
		}
		auto rPos = alignObj.findReversePrimer(readEnd.seq_,
				rev.infoRC_.seq_);
		if (trimBackSize < len(info)) {
			rPos.first += len(info) - trimBackSize;
			rPos.second += len(info) - trimBackSize;
		}
		alignObj.rearrangeObjs(readEnd, rev.infoRC_,
				true);
		alignObj.profilePrimerAlignment(readEnd,
				rev.infoRC_);

		if (alignObj.comp_.distances_.query_.coverage_
				< pars.allowable_.distances_.query_.coverage_
				|| !pars.allowable_.passErrorProfile(alignObj.comp_)) {
		} else {
			determinedPrimers.emplace_back(rPos.first, rPos.second, rev.info_.name_,
					alignObj.comp_);
		}

//		info.on_ = primerGood;
//
//		if (primerGood) {
//			if (pars.primerToLowerCase_) {
//				if(pars.trimExtra_){
//					while (rPos.first != 0 && info.seq_[rPos.first] == info.seq_[rPos.first - 1]) {
//						--rPos.first;
//					}
//				}
//				changeSubStrToLowerToEnd(info.seq_, rPos.first);
//			}
//			info.setClip(0, rPos.second);
//		}
//		if(primerGood){
//			break;
//		}
	}


	PrimerPositionScore bestPrimer;
	bool primerGood = true;

	if (0 == determinedPrimers.size()) {
		primerGood = false;
	} else if (1 == determinedPrimers.size()) {
		primerGood = true;
		bestPrimer = determinedPrimers.front();
	} else {
		primerGood =true;
		bestPrimer = determinedPrimers.front();
		for(const auto pos : iter::range<uint32_t>(1, determinedPrimers.size())){
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

	info.on_ = primerGood;

	if (primerGood) {
		if (pars.primerToLowerCase_) {
			if(pars.trimExtra_){
				while (bestPrimer.start_!= 0 && info.seq_[bestPrimer.start_] == info.seq_[bestPrimer.start_ - 1]) {
					--bestPrimer.start_;
				}
			}
			changeSubStrToLowerToEnd(info.seq_, bestPrimer.start_);
		}
		info.setClip(0, bestPrimer.end_);
	}
	return primerGood;
}



bool PrimerDeterminator::checkForForwardPrimerInRev(seqInfo & info, const std::string & primerName,const PrimerDeterminatorPars & pars,
		aligner & alignObj) const{
	if (!njh::in(primerName, primers_))  {
		throw std::runtime_error { std::string(__PRETTY_FUNCTION__) + ": No primer info for: "
				+ primerName };
	}
	std::vector<PrimerPositionScore> determinedPrimers;
	for(const auto & fwd : primers_.at(primerName).fwds_){
		seqInfo readEnd;
		auto trimBackSize = pars.primerWithin_ + fwd.infoRC_.seq_.size() * 2;
		if (trimBackSize < len(info)) {
			readEnd = info.getSubRead(len(info) - trimBackSize, len(info) - pars.primerStart_);
		} else {
			readEnd = info;
		}

		auto rPos = alignObj.findReversePrimer(readEnd.seq_,
				fwd.infoRC_.seq_);
		if (trimBackSize < len(info)) {
			rPos.first += len(info) - trimBackSize;
			rPos.second += len(info) - trimBackSize;
		}
		alignObj.rearrangeObjs(readEnd, fwd.infoRC_,
				true);
		alignObj.profilePrimerAlignment(readEnd,
				fwd.infoRC_);

		if (alignObj.comp_.distances_.query_.coverage_
				< pars.allowable_.distances_.query_.coverage_
				|| !pars.allowable_.passErrorProfile(alignObj.comp_)) {
		}else{
			determinedPrimers.emplace_back(rPos.first, rPos.second, fwd.info_.name_,
					alignObj.comp_);
		}
//		info.on_ = primerGood;
//
//		if (primerGood) {
//			if (pars.primerToLowerCase_) {
//				if(pars.trimExtra_){
//					while (rPos.first != 0 && info.seq_[rPos.first] == info.seq_[rPos.first - 1]) {
//						--rPos.first;
//					}
//				}
//				changeSubStrToLowerToEnd(info.seq_, rPos.first);
//			}
//			info.setClip(0, rPos.second);
//		}
//		if(primerGood){
//			break;
//		}
	}

	PrimerPositionScore bestPrimer;
	bool primerGood = true;

	if (determinedPrimers.empty()) {
		primerGood = false;
	} else if (1 == determinedPrimers.size()) {
		primerGood = true;
		bestPrimer = determinedPrimers.front();
	} else {
		primerGood = true;
		bestPrimer = determinedPrimers.front();
		for(const auto pos : iter::range<uint32_t>(1, determinedPrimers.size())){
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

	info.on_ = primerGood;
	if (primerGood) {
		if (pars.primerToLowerCase_) {
			if(pars.trimExtra_){
				while (bestPrimer.start_!= 0 && info.seq_[bestPrimer.start_] == info.seq_[bestPrimer.start_ - 1]) {
					--bestPrimer.start_;
				}
			}
			changeSubStrToLowerToEnd(info.seq_, bestPrimer.start_);
		}
		info.setClip(0, bestPrimer.end_);
	}
	return primerGood;
}


} /* namespace njhseq */
