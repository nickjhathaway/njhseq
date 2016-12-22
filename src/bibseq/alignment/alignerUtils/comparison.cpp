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
#include "bibseq/alignment/alignerUtils/comparison.hpp"

namespace bibseq {
void DistanceMet::reset() {
	identities_ = 0;
	identity_ = 0;
	covered_ = 0;
	coverage_ = 0;
}

Json::Value DistanceMet::toJson()const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::DistanceMet");
	ret["identities_"] = bib::json::toJson(identities_);
	ret["identity_"] = bib::json::toJson(identity_);
	ret["coverage_"] = bib::json::toJson(coverage_);
	ret["covered_"] = bib::json::toJson(covered_);
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
	ret["class"] = bib::json::toJson("bibseq::DistanceComp");
	ret["basesInAln_"] = bib::json::toJson(basesInAln_);
	ret["percentMismatch_"] = bib::json::toJson(percentMismatch_);
	ret["percentMatch_"] = bib::json::toJson(percentMatch_);
	ret["percentGaps_"] = bib::json::toJson(percentGaps_);
	ret["overLappingEvents_"] = bib::json::toJson(overLappingEvents_);
	ret["eventBasedIdentity_"] = bib::json::toJson(eventBasedIdentity_);
	ret["overLappingEventsHq_"] = bib::json::toJson(overLappingEventsHq_);
	ret["eventBasedIdentityHq_"] = bib::json::toJson(eventBasedIdentityHq_);
	ret["ref_"] = bib::json::toJson(ref_);
	ret["query_"] = bib::json::toJson(query_);
	ret["mismatches_"] = bib::json::toJson(mismatches_);
	ret["lowKmerMismatches_"] = bib::json::toJson(lowKmerMismatches_);
	ret["alignmentGaps_"] = bib::json::toJson(alignmentGaps_);
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
	distances_.overLappingEventsHq_ = highQualityMatches_ + hqMismatches_
			+ indelEvents;
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
  return generatedError.distances_.eventBasedIdentity_ >= distances_.eventBasedIdentity_;
}

bool comparison::passIdThresholdHq(const comparison& generatedError) const {
  return generatedError.distances_.eventBasedIdentityHq_ >= distances_.eventBasedIdentityHq_;
}


Json::Value comparison::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::comparison");
	ret["oneBaseIndel_"] = bib::json::toJson(oneBaseIndel_);
	ret["twoBaseIndel_"] = bib::json::toJson(twoBaseIndel_);
	ret["largeBaseIndel_"] = bib::json::toJson(largeBaseIndel_);
	ret["hqMismatches_"] = bib::json::toJson(hqMismatches_);
	ret["lqMismatches_"] = bib::json::toJson(lqMismatches_);
	ret["lowKmerMismatches_"] = bib::json::toJson(lowKmerMismatches_);
	ret["highQualityMatches_"] = bib::json::toJson(highQualityMatches_);
	ret["lowQualityMatches_"] = bib::json::toJson(lowQualityMatches_);
	ret["distances_"] = bib::json::toJson(distances_);
	ret["refName_"] = bib::json::toJson(refName_);
	ret["queryName_"] = bib::json::toJson(queryName_);
	return ret;
}

std::ostream & operator <<(std::ostream & out, const comparison & comp) {
	Json::FastWriter jWriter;
	out << jWriter.write(comp.toJson()) << std::endl;
	return out;
}

}  // namespace bibseq
