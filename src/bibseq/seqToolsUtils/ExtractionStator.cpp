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
//

/*
 * ExtractionStator.cpp
 *
 *  Created on: Sep 29, 2015
 *      Author: nick
 */


#include "ExtractionStator.hpp"


namespace bibseq {

ExtractionStator::ExtractionStator() {

}
;
ExtractionStator::ExtractionStator(uint32_t totalReadCount,
		uint32_t readsUnrecBarcode, uint32_t smallFrags) :
		totalReadCount_(totalReadCount), readsUnrecBarcode_(readsUnrecBarcode), smallFrags_(
				smallFrags) {

}


void ExtractionStator::increaseFailedForward(const std::string & midName, const std::string & seqName){
	bool rComp = containsSubString(seqName, "_Comp");
	++failedForward_[midName][rComp];
}

void ExtractionStator::increaseCounts(const std::string & midName, const std::string & seqName,
		extractCase eCase) {
	bool rComp = containsSubString(seqName, "_Comp");
	switch (eCase) {
	case extractCase::GOOD:
		++counts_[midName][rComp].good_;
		break;
	case extractCase::BADREVERSE:
		++counts_[midName][rComp].bad_;
		++counts_[midName][rComp].badReverse_;
		break;
	case extractCase::CONTAINSNS:
		++counts_[midName][rComp].bad_;
		++counts_[midName][rComp].containsNs_;
		break;
	case extractCase::MINLENBAD:
		++counts_[midName][rComp].bad_;
		++counts_[midName][rComp].minLenBad_;
		break;
	case extractCase::MAXLENBAD:
		++counts_[midName][rComp].bad_;
		++counts_[midName][rComp].maxLenBad_;
		break;
	case extractCase::QUALITYFAILED:
		++counts_[midName][rComp].bad_;
		++counts_[midName][rComp].qualityFailed_;
		break;
	case extractCase::CONTAMINATION:
		++counts_[midName][rComp].contamination_;
		break;
	default:
		std::stringstream ss;
		ss << bib::bashCT::boldBlack(__PRETTY_FUNCTION__)
			 << bib::bashCT::boldRed(": shouldn't be happending..., unknown case: ")
			 << std::endl;
		throw std::runtime_error{ss.str()};
		break;
	}
}
void ExtractionStator::outStatsPerName(std::ostream & out, const std::string & delim){

	for(auto & mid : counts_){
		uint32_t totalReads = mid.second[true].getTotal() + mid.second[false].getTotal();
		uint32_t totalGoodReads = mid.second[true].good_ + mid.second[false].good_;
		uint32_t totalBadReads = mid.second[true].bad_ + mid.second[false].bad_;
		uint32_t totalContam = mid.second[true].contamination_ + mid.second[false].contamination_;
		uint32_t totalBadRev = mid.second[true].badReverse_ + mid.second[false].badReverse_;
		uint32_t totalConN = mid.second[true].containsNs_ + mid.second[false].containsNs_;
		uint32_t totalMinLen = mid.second[true].minLenBad_ + mid.second[false].minLenBad_;
		uint32_t totalMaxLen = mid.second[true].maxLenBad_ + mid.second[false].maxLenBad_;
		uint32_t totalQualFail = mid.second[true].qualityFailed_ + mid.second[false].qualityFailed_;
		out << vectorToString(toVecStr(mid.first, totalReads,
				getPercentageString(totalGoodReads, totalReads),
				getPercentageString(mid.second[false].good_, totalGoodReads),
				getPercentageString(mid.second[true].good_, totalGoodReads),
				getPercentageString(totalBadReads, totalReads),
				getPercentageString(totalBadRev, totalBadReads),
				getPercentageString(totalConN, totalBadReads),
				getPercentageString(totalMinLen, totalBadReads),
				getPercentageString(totalMaxLen, totalBadReads),
				getPercentageString(totalQualFail, totalBadReads),
				getPercentageString(totalContam, totalReads)), delim) << "\n";

	}
}
void ExtractionStator::outFailedForwardStats(std::ostream & out, const std::string & delim){
	for(auto & ff : failedForward_){
		uint32_t total = 0;
		for(auto & mid : counts_){
			if(endsWith(mid.first, ff.first)){
				total += mid.second[true].getTotal() + mid.second[false].getTotal();
			}
		}
		total += ff.second[true] + ff.second[false];
		out << ff.first
				<< delim << getPercentageString(ff.second[true] + ff.second[false], total)
				<< delim << getPercentageString(ff.second[true] ,ff.second[true] + ff.second[false])
				<< delim << getPercentageString(ff.second[false],ff.second[true] + ff.second[false]) << std::endl;
	}
}

void ExtractionStator::outTotalStats(std::ostream & out, const std::string & delim) {
	uint32_t totalBadReads = 0;
	uint32_t totalGoodReads = 0;
	uint32_t totalContam = 0;
	uint32_t totalFailedForward = 0;
	for(auto & ff : failedForward_){
		totalFailedForward += ff.second[true] + ff.second[false];
	}
	for (auto & mid : counts_) {
		totalGoodReads += mid.second[true].good_ + mid.second[false].good_;
		totalBadReads += mid.second[true].bad_ + mid.second[false].bad_;
		totalContam += mid.second[true].contamination_
				+ mid.second[false].contamination_;
	}
	out << vectorToString(toVecStr(totalReadCount_,
			getPercentageString(readsUnrecBarcode_, totalReadCount_),
			getPercentageString(smallFrags_, totalReadCount_),
			getPercentageString(totalFailedForward, totalReadCount_),
			getPercentageString(totalBadReads, totalReadCount_),
			getPercentageString(totalGoodReads, totalReadCount_),
			getPercentageString(totalContam, totalReadCount_)), delim) << "\n";
}

}  // namespace bibseq
