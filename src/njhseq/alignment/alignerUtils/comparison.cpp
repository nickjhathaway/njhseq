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
#include "njhseq/alignment/alignerUtils/comparison.hpp"

namespace njhseq {
void DistanceMet::reset() {
	identities_ = 0;
	identity_ = 0;
	covered_ = 0;
	coverage_ = 0;
}

Json::Value DistanceMet::toJson()const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["identities_"] = njh::json::toJson(identities_);
	ret["identity_"] = njh::json::toJson(identity_);
	ret["coverage_"] = njh::json::toJson(coverage_);
	ret["covered_"] = njh::json::toJson(covered_);
	return ret;
}




void DistanceComp::reset() {
	basesInAln_ = 0;
	percentMismatch_ = 0;
	percentMatch_ = 0;
	percentGaps_ = 0;
	overLappingEvents_ = 0;
	eventBasedIdentity_ = 0;
	overLappingEventsHq_ = 0;
	eventBasedIdentityHq_ = 0;


	ref_.reset();
	query_.reset();
	mismatches_.clear();
	lowKmerMismatches_.clear();
	alignmentGaps_.clear();

}

uint32_t DistanceComp::getNumOfEvents(bool countLowKmer)const{
	if(countLowKmer){
		return alignmentGaps_.size() + mismatches_.size() + lowKmerMismatches_.size();
	}
	return alignmentGaps_.size() + mismatches_.size();
}


Json::Value DistanceComp::toJson()const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["basesInAln_"] = njh::json::toJson(basesInAln_);
	ret["percentMismatch_"] = njh::json::toJson(percentMismatch_);
	ret["percentMatch_"] = njh::json::toJson(percentMatch_);
	ret["percentGaps_"] = njh::json::toJson(percentGaps_);
	ret["overLappingEvents_"] = njh::json::toJson(overLappingEvents_);
	ret["eventBasedIdentity_"] = njh::json::toJson(eventBasedIdentity_);
	ret["overLappingEventsHq_"] = njh::json::toJson(overLappingEventsHq_);
	ret["eventBasedIdentityHq_"] = njh::json::toJson(eventBasedIdentityHq_);
	ret["ref_"] = njh::json::toJson(ref_);
	ret["query_"] = njh::json::toJson(query_);
	ret["mismatches_"] = njh::json::toJson(mismatches_);
	ret["lowKmerMismatches_"] = njh::json::toJson(lowKmerMismatches_);
	ret["alignmentGaps_"] = njh::json::toJson(alignmentGaps_);
	return ret;
}

void comparison::setEventBaseIdentity() {
	distances_.overLappingEvents_ = highQualityMatches_ + lowQualityMatches_
			+ hqMismatches_ + lqMismatches_ + lowKmerMismatches_
			+ distances_.alignmentGaps_.size();
	if (distances_.overLappingEvents_ == 0) {
		distances_.eventBasedIdentity_ = 0;
	} else {
		distances_.eventBasedIdentity_ = (highQualityMatches_ + lowQualityMatches_)
				/ static_cast<double>(distances_.overLappingEvents_);
	}

}

void comparison::setEventBaseIdentityHq() {
	//distances_.overLappingEventsHq_ = highQualityMatches_ + hqMismatches_
	//		+ distances_.alignmentGaps_.size();
	//distances_.overLappingEventsHq_ = highQualityMatches_ + hqMismatches_
	//		+ lowKmerMismatches_ + distances_.alignmentGaps_.size();

	//high quality events being just high quality mismatches and high quality mismatches and indels
	//if weighing for indel in homopolymer this will be taken into account and so will low freq k-mer mismatches;
	double indelEvents = oneBaseIndel_ + twoBaseIndel_ + largeBaseIndel_;
	distances_.overLappingEventsHq_ = highQualityMatches_ + hqMismatches_ + indelEvents;
	//distances_.overLappingEventsHq_ = highQualityMatches_ + hqMismatches_
	//		+ distances_.alignmentGaps_.size();
	if (distances_.overLappingEventsHq_ == 0) {
		distances_.eventBasedIdentityHq_ = 0;
	} else {
		distances_.eventBasedIdentityHq_ = (highQualityMatches_)
				/ static_cast<double>(distances_.overLappingEventsHq_);
	}
}

void comparison::recalcMismatchQuality(const QualScorePars & pars){
  hqMismatches_ = 0;
  lqMismatches_ = 0;
  for(const auto & mis : distances_.mismatches_){
  	if(mis.second.highQuality(pars)){
  		++hqMismatches_;
  	}else{
  		++lqMismatches_;
  	}
  }
}

void comparison::resetCounts() {
  oneBaseIndel_ = 0;
  twoBaseIndel_ = 0;
  largeBaseIndel_ = 0;
  hqMismatches_ = 0;
  lqMismatches_ = 0;
  lowKmerMismatches_ = 0;

  highQualityMatches_ = 0;
  lowQualityMatches_ = 0;

  refName_ = "";
  queryName_ = "";

  distances_.reset();
}


bool comparison::passIdAndErrorThreshold(const comparison& generatedError) const {
	return (oneBaseIndel_ >= generatedError.oneBaseIndel_ &&
	          twoBaseIndel_ >= generatedError.twoBaseIndel_ &&
	          largeBaseIndel_ >= generatedError.largeBaseIndel_ &&
	          hqMismatches_ >= generatedError.hqMismatches_ &&
	          lqMismatches_ >= generatedError.lqMismatches_ &&
	          lowKmerMismatches_ >= generatedError.lowKmerMismatches_ &&
						generatedError.distances_.eventBasedIdentity_ >= distances_.eventBasedIdentity_);
}

bool comparison::passErrorProfile(const comparison& generatedError) const {
  return (oneBaseIndel_ >= generatedError.oneBaseIndel_ &&
          twoBaseIndel_ >= generatedError.twoBaseIndel_ &&
          largeBaseIndel_ >= generatedError.largeBaseIndel_ &&
          hqMismatches_ >= generatedError.hqMismatches_ &&
          lqMismatches_ >= generatedError.lqMismatches_ &&
          lowKmerMismatches_ >= generatedError.lowKmerMismatches_);
}


bool comparison::passIdThreshold(const comparison& generatedError) const {
	//rounding to precision possible with normal read length for short amplicons;
	/**@todo determine this with input read length*/
  return roundDecPlaces(generatedError.distances_.eventBasedIdentity_, 3)  >= distances_.eventBasedIdentity_;
}

bool comparison::passIdThresholdHq(const comparison& generatedError) const {
	//rounding to precision possible with normal read length for short amplicons;
	/**@todo determine this with input read length*/
  return roundDecPlaces(generatedError.distances_.eventBasedIdentityHq_,3) >= distances_.eventBasedIdentityHq_;
}


Json::Value comparison::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["oneBaseIndel_"] = njh::json::toJson(oneBaseIndel_);
	ret["twoBaseIndel_"] = njh::json::toJson(twoBaseIndel_);
	ret["largeBaseIndel_"] = njh::json::toJson(largeBaseIndel_);
	ret["hqMismatches_"] = njh::json::toJson(hqMismatches_);
	ret["lqMismatches_"] = njh::json::toJson(lqMismatches_);
	ret["lowKmerMismatches_"] = njh::json::toJson(lowKmerMismatches_);
	ret["highQualityMatches_"] = njh::json::toJson(highQualityMatches_);
	ret["lowQualityMatches_"] = njh::json::toJson(lowQualityMatches_);
	ret["distances_"] = njh::json::toJson(distances_);
	ret["refName_"] = njh::json::toJson(refName_);
	ret["queryName_"] = njh::json::toJson(queryName_);
	return ret;
}

std::ostream & operator <<(std::ostream & out, const comparison & comp) {
	out << njh::json::writeAsOneLine(comp.toJson()) << std::endl;
	return out;
}

}  // namespace njhseq
