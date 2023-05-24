
/*
 * Primer3Results.cpp
 *
 *  Created on: Aug 16, 2017
 *      Author: nick
 */
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


#include "Primer3Results.hpp"
#include "njhseq/IO/InputStream.hpp"
#include "njhseq/helpers/seqUtil.hpp"

namespace njhseq {



void Primer3Runner::Primer3Options::setPrimerSizeOpts(seqSetUp & setUp){
	setUp.setOption(PRIMER_MAX_SIZE, "--PRIMER_MAX_SIZE", "PRIMER max SIZE");
	setUp.setOption(PRIMER_MIN_SIZE, "--PRIMER_MIN_SIZE", "PRIMER min SIZE");
	setUp.setOption(PRIMER_OPT_SIZE, "--PRIMER_OPT_SIZE", "PRIMER optimal SIZE");

	if(PRIMER_MIN_SIZE > PRIMER_MAX_SIZE){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("PRIMER_MIN_SIZE, ", PRIMER_MIN_SIZE, ", can't be larger than PRIMER_MAX_SIZE", PRIMER_MAX_SIZE));
	}

	if(PRIMER_OPT_SIZE < PRIMER_MIN_SIZE){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("PRIMER_OPT_SIZE, ", PRIMER_OPT_SIZE, ", can't be smaller than PRIMER_MIN_SIZE", PRIMER_MIN_SIZE));
	}
	if(PRIMER_OPT_SIZE > PRIMER_MAX_SIZE){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("PRIMER_OPT_SIZE, ", PRIMER_OPT_SIZE, ", can't be larger than PRIMER_MAX_SIZE", PRIMER_MAX_SIZE));
	}
}

void Primer3Runner::Primer3Options::setInsertSizeOptions(seqSetUp & setUp){
	setUp.setOption(minSize, "--minAmpSize", "min amplicon Size");
	setUp.setOption(maxSize, "--maxAmpSize", "max amplicon Size");
	if(minSize > maxSize){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("minSize, ", minSize, ", can't be larger than maxSize", maxSize));
	}
}

void Primer3Runner::Primer3Options::setPrimaryOptions(seqSetUp & setUp){
	setUp.setOption(primer3ConfigPath, "--primer3ConfigPath", "primer3 Config Path for the thermodynamic data", true);
	primer3ConfigPath = njh::appendAsNeededRet(primer3ConfigPath.string(), "/");

	setUp.setOption(PRIMER_OPT_GC_PERCENT, "--PRIMER_OPT_GC_PERCENT", "optimal primer GC content");
	setUp.setOption(PRIMER_INTERNAL_MIN_GC, "--PRIMER_INTERNAL_MIN_GC", "minimal internal GC content");
}

void Primer3Runner::Primer3Options::setReturnOptions(seqSetUp & setUp){
	setUp.setOption(PRIMER_NUM_RETURN, "--PRIMER_NUM_RETURN", "PRIMER NUM RETURN");
}


void Primer3Runner::Primer3Options::loadAdditionalOptions(const bfs::path & jsonOpts){
	auto values =	njh::json::parseFile(jsonOpts.string());
	additionalOpts = njh::json::JsonToOMap<std::string,std::string>(values);
}


Json::Value Primer3Runner::Primer3Options::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["PRIMER_MAX_SIZE"] = njh::json::toJson(PRIMER_MAX_SIZE);
	ret["PRIMER_MIN_SIZE"] = njh::json::toJson(PRIMER_MIN_SIZE);
	ret["PRIMER_OPT_SIZE"] = njh::json::toJson(PRIMER_OPT_SIZE);
	ret["PRIMER_OPT_GC_PERCENT"] = njh::json::toJson(PRIMER_OPT_GC_PERCENT);
	ret["PRIMER_INTERNAL_MIN_GC"] = njh::json::toJson(PRIMER_INTERNAL_MIN_GC);
	ret["primer3ConfigPath"] = njh::json::toJson(primer3ConfigPath);
	ret["minSize"] = njh::json::toJson(minSize);
	ret["maxSize"] = njh::json::toJson(maxSize);
	ret["PRIMER_NUM_RETURN"] = njh::json::toJson(PRIMER_NUM_RETURN);
	ret["task"] = njh::json::toJson(task);

	ret["additionalOpts"] = njh::json::toJson(additionalOpts);
	return ret;
}


Json::Value Primer3Runner::region::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["start_"] =njh::json::toJson(start_);
	ret["size_"] =njh::json::toJson(size_);
	return ret;
}


GenomicRegion Primer3Runner::Primer::genRegion(const std::string & templateId) const {
	return GenomicRegion(
			name_,
			templateId,
			forwardOrientationPos_.start_,
			forwardOrientationPos_.start_ + forwardOrientationPos_.size_,
			right_);
}


Json::Value Primer3Runner::Primer::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["name_"] = njh::json::toJson(name_);
	ret["penalty_"] = njh::json::toJson(penalty_);
	ret["seq_"] = njh::json::toJson(seq_);
	ret["problems_"] = njh::json::toJson(problems_);
	ret["tm_"] = njh::json::toJson(tm_);
	ret["gc_percent_"] = njh::json::toJson(gc_percent_);
	ret["self_any_th_"] = njh::json::toJson(self_any_th_);
	ret["self_end_th_"] = njh::json::toJson(self_end_th_);
	ret["hairpin_th_"] = njh::json::toJson(hairpin_th_);
	ret["end_stability_"] = njh::json::toJson(end_stability_);
	ret["right_"] = njh::json::toJson(right_);
	ret["originalPos_"] = njh::json::toJson(originalPos_);
	ret["forwardOrientationPos_"] = njh::json::toJson(forwardOrientationPos_);
	return ret;
}



uint32_t Primer3Runner::PrimerPair::getStart() const{
	return left_.forwardOrientationPos_.start_;
}
uint32_t Primer3Runner::PrimerPair::getEnd() const{
	return right_.forwardOrientationPos_.start_ + right_.forwardOrientationPos_.size_;
}

bool Primer3Runner::PrimerPair::overlaps(const PrimerPair & other, uint32_t minOverlap)const{
	return getOverlapLen(other) >= minOverlap;
}

uint32_t Primer3Runner::PrimerPair::getOverlapLen(const PrimerPair & other)const{
	uint32_t overlapLength = 0;


	if( (other.getStart() >= getStart() && other.getStart() < getEnd()) ||
			(other.getEnd() > getStart() && other.getEnd() <= getEnd()) ||
			(getStart() >= other.getStart() &&  getStart() <  other.getEnd()) ||
			(getEnd() > other.getStart() && getEnd() <= other.getEnd())  ) {

		auto overlapStart = std::max(other.getStart(), getStart());
		auto overlapEnd = std::min(other.getEnd(), getEnd());
		overlapLength =  overlapEnd - overlapStart;
	}
	return overlapLength;
}



Json::Value Primer3Runner::PrimerPair::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["name_"] = njh::json::toJson(name_);
	ret["left_"] = njh::json::toJson(left_);
	ret["right_"] = njh::json::toJson(right_);
	return ret;
}


Json::Value Primer3Runner::PrimerPairGeneric::toJson() const{
	Json::Value ret = PrimerPair::toJson();
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["compl_any_th_"] = njh::json::toJson(compl_any_th_);
	ret["compl_end_th_"] = njh::json::toJson(compl_end_th_);
	ret["penalty_"] = njh::json::toJson(penalty_);
	ret["product_size_"] = njh::json::toJson(product_size_);
	return ret;
}


Json::Value Primer3Runner::Primer3Results::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["sequence_id_"] =njh::json::toJson(sequence_id_);
	ret["sequence_template_"] =njh::json::toJson(sequence_template_);
	ret["sequence_target_"] =njh::json::toJson(sequence_target_);
	ret["sequence_excluded_region_"] =njh::json::toJson(sequence_excluded_region_);
	ret["primer_left_num_returned_"] =njh::json::toJson(primer_left_num_returned_);
	ret["primer_right_num_returned_"] =njh::json::toJson(primer_right_num_returned_);
	ret["primer_internal_num_returned_"] =njh::json::toJson(primer_internal_num_returned_);
	ret["primer_pair_num_returned_"] =njh::json::toJson(primer_pair_num_returned_);

	ret["warnings_"] =njh::json::toJson(warnings_);
	ret["misc_"] =njh::json::toJson(misc_);

	return ret;
}

Json::Value Primer3Runner::Primer3ResultsGeneric::toJson() const{
	Json::Value ret = Primer3Results::toJson();

	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["primers_"] =njh::json::toJson(primerPairs_);

	return ret;
}

Json::Value Primer3Runner::Primer3ResultsList::toJson() const{
	Json::Value ret = Primer3Results::toJson();

	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["leftPrimers_"] =njh::json::toJson(leftPrimers_);
	ret["rightPrimers_"] =njh::json::toJson(rightPrimers_);

	return ret;
}



std::unordered_map<std::string, GenomicRegion> Primer3Runner::Primer3ResultsGeneric::genPrimersRegions() const {
	std::unordered_map<std::string, GenomicRegion> ret;
	for (const auto &p : primerPairs_) {
		ret[p->right_.name_] = p->right_.genRegion(sequence_id_);
		ret[p->left_.name_] = p->left_.genRegion(sequence_id_);
	}
	return ret;
}

std::unordered_map<std::string, GenomicRegion> Primer3Runner::Primer3ResultsList::genPrimersRegions() const {
	std::unordered_map<std::string, GenomicRegion> ret;
	for (const auto &p : leftPrimers_) {
		ret[p->name_] = p->genRegion(sequence_id_);
	}
	for (const auto &p : rightPrimers_) {
		ret[p->name_] = p->genRegion(sequence_id_);
	}
	return ret;
}



std::shared_ptr<Primer3Runner::PrimerPairGeneric> Primer3Runner::Primer3ResultsGeneric::getLowestPenaltyPair() const{
	std::shared_ptr<Primer3Runner::PrimerPairGeneric> ret;
	double lowestPenality = std::numeric_limits<double>::max();
	for(const auto & pair : primerPairs_){
		if(pair->penalty_ < lowestPenality){
			lowestPenality = pair->penalty_;
			ret = pair;
		}
	}
	return ret;
}
std::vector<std::shared_ptr<Primer3Runner::PrimerPairGeneric>> Primer3Runner::Primer3ResultsGeneric::getLowestPenaltyNonOverlappingPairs() const{
	std::vector<std::shared_ptr<Primer3Runner::PrimerPairGeneric>> ret;

	std::vector<uint32_t> primerPositons(primerPairs_.size());
	njh::iota(primerPositons, 0U);

	njh::sort(primerPositons,[this](const uint32_t & pair1Pos, const uint32_t & pair2Pos){
		return primerPairs_[pair1Pos]->penalty_ <= primerPairs_[pair2Pos]->penalty_;
	});

	for(const auto & pairPos : primerPositons){
		bool pass = true;
		for(const auto & already : ret){
			if(primerPairs_[pairPos]->overlaps(*already)){
				pass = false;
				break;
			}
		}

		if(pass){
			ret.emplace_back(primerPairs_[pairPos]);
		}
	}

	return ret;
}



std::vector<std::shared_ptr<Primer3Runner::Primer3ResultsList>> Primer3Runner::Primer3ResultsList::parsePrimer3OutputResults(const bfs::path & input, bool ignoreUnexpected){
	std::vector<std::shared_ptr<Primer3Runner::Primer3ResultsList>> results;


	InputStream in(input);
	std::string line;
	std::shared_ptr<Primer3Runner::Primer3ResultsList> currentResult = std::make_shared<Primer3Runner::Primer3ResultsList>();
	std::shared_ptr<Primer> currentLeft = std::make_shared<Primer>();
	std::shared_ptr<Primer> currentRight= std::make_shared<Primer>();

	currentResult->leftPrimers_.emplace_back(currentLeft);
	currentResult->rightPrimers_.emplace_back(currentRight);

	uint32_t currentLeftPrimerPairID = 0;
	uint32_t currentRightPrimerPairID = 0;



	std::regex primerInfoPat { "^(PRIMER_(LEFT|RIGHT|PAIR)_([0-9]+))_?(.*)=.*" };
	while (njh::files::crossPlatGetline(in, line)) {
		if ("=" == line) {
			results.emplace_back(currentResult);
			currentResult = std::make_shared<Primer3ResultsList>();
			currentLeft = std::make_shared<Primer>();
			currentRight = std::make_shared<Primer>();
			currentResult->leftPrimers_.emplace_back(currentLeft);
			currentResult->rightPrimers_.emplace_back(currentRight);
			currentLeftPrimerPairID = 0;
			currentRightPrimerPairID = 0;
		} else {
			auto toks = tokenizeString(line, "=");
			if (1 == toks.size() && '-' != line.back()) {
				//empty value
				continue;
			}
			if (2 != toks.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__
						<< ", error line should be key=value, didn't find two values separated by a =, found "
						<< toks.size() << " objects instead for line " << line << "\n";
				throw std::runtime_error { ss.str() };
			}
			std::smatch match;
			if (std::regex_match(line, match, primerInfoPat)) {
				//matched primer
				std::string name = match[1];
				std::string leftOrRightOrPair = match[2];
				std::string pairName = match[3];

				uint32_t primerNameNum = njh::StrToNumConverter::stoToNum<uint32_t>(pairName);
				bool right = leftOrRightOrPair == "RIGHT";

				if (right) {
					if (primerNameNum != currentRightPrimerPairID) {
						currentRight = std::make_shared<Primer>();
						currentRightPrimerPairID = primerNameNum;
						currentResult->rightPrimers_.emplace_back(currentRight);
					}
				} else {
					if (primerNameNum != currentLeftPrimerPairID) {
						currentLeft = std::make_shared<Primer>();
						currentLeftPrimerPairID = primerNameNum;
						currentResult->leftPrimers_.emplace_back(currentLeft);
					}
				}


				if (leftOrRightOrPair != "LEFT" && leftOrRightOrPair != "RIGHT" && leftOrRightOrPair != "PAIR") {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", unexpected designation for line: " << line
							<< ", was expecting RIGHT or LEFT not: " << match[2] << "\n";
					throw std::runtime_error { ss.str() };
				}

				if ("" == match[4]) {
					//position
					auto positionToks = tokenizeString(toks[1], ",");
					if (2 != positionToks.size()) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__
								<< ", error position values should be two values separated by a \",\", found "
								<< toks.size() << " objects instead for line " << line
								<< " for tok: " << toks[1] << "\n";
						throw std::runtime_error { ss.str() };
					}
					auto start = njh::StrToNumConverter::stoToNum<uint32_t>(positionToks[0]);
					auto size = njh::StrToNumConverter::stoToNum<uint32_t>(positionToks[1]);
					if (right) {
						currentRight->originalPos_.start_ = start;
						currentRight->originalPos_.size_ = size;
						currentRight->forwardOrientationPos_.start_ = start+ 1 - size;
						currentRight->forwardOrientationPos_.size_ = size;
						currentRight->name_ = name;
					}else{
						currentLeft->originalPos_.start_ = start;
						currentLeft->originalPos_.size_ = size;
						currentLeft->forwardOrientationPos_.start_ = start;
						currentLeft->forwardOrientationPos_.size_ = size;
						currentLeft->name_ = name;
					}
				} else if ("PENALTY" == match[4]) {
					if (leftOrRightOrPair == "RIGHT") {
						currentRight->penalty_= njh::StrToNumConverter::stoToNum<double>(toks[1]);
					} else {
						currentLeft->penalty_= njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("SEQUENCE" == match[4]) {
					if(right){
						currentRight->seq_ = toks[1];
					}else{
						currentLeft->seq_ = toks[1];
					}
				}  else if ("PROBLEMS" == match[4]) {
					if(right){
						currentRight->problems_ = njh::tokenizeString(toks[1], ";");
					}else{
						currentLeft->problems_ = njh::tokenizeString(toks[1], ";");
					}
				} else if ("TM" == match[4]) {
					if(right){
						currentRight->tm_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentLeft->tm_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("GC_PERCENT" == match[4]) {
					if(right){
						currentRight->gc_percent_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentLeft->gc_percent_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("SELF_ANY_TH" == match[4]) {
					if(right){
						currentRight->self_any_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentLeft->self_any_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}

				} else if ("SELF_END_TH" == match[4]) {
					if(right){
						currentRight->self_end_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentLeft->self_end_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("HAIRPIN_TH" == match[4]) {
					if(right){
						currentRight->hairpin_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentLeft->hairpin_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("END_STABILITY" == match[4]) {
					if(right){
						currentRight->end_stability_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentLeft->end_stability_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else {
					if(ignoreUnexpected){
						currentResult->misc_.emplace(toks[0], toks[1]);
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", unexpected sub case for : " << line
								<< ", case: " << match[4] << "\n";
						throw std::runtime_error { ss.str() };
					}
				}
			} else {
				/*
				 * 	std::string ;
	std::string ;
	std::string ;
				 */
				if ("SEQUENCE_ID" == toks[0]) {
					currentResult->sequence_id_ = toks[1];
				} else if ("SEQUENCE_TEMPLATE" == toks[0]) {
					currentResult->sequence_template_ = toks[1];
				} else if ("PRIMER_WARNING" == toks[0]) {
					currentResult->warnings_ = tokenizeString(toks[1], ";"); //=pick_sequencing_primers
				} else if ("PRIMER_TASK" == toks[0]) {
					currentResult->primer_task_ = toks[1];
				} else if ("PRIMER_LEFT_EXPLAIN" == toks[0]) {
					currentResult->primer_left_explain_ = toks[1];
				} else if ("PRIMER_RIGHT_EXPLAIN" == toks[0]) {
					currentResult->primer_right_explain_ = toks[1];
				} else if ("PRIMER_PAIR_EXPLAIN" == toks[0]) {
					currentResult->primer_pair_explain_ = toks[1];
				} else if ("SEQUENCE_TARGET" == toks[0]) {
					auto secondToks = tokenizeString(toks[1], " ");
					for (const auto & tok : secondToks) {
						auto regToks = tokenizeString(tok, ",");
						if (2 != regToks.size()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__
									<< ", error position values should be two values separated by a \",\" , found "
									<< toks.size() << " objects instead for line " << line
									<< " for tok: " << toks[1] << " for sub toks: " << tok
									<< "\n";
							throw std::runtime_error { ss.str() };
						}
						currentResult->sequence_target_.emplace_back(
								Primer3Runner::region { njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[0]), njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[1]) });
					}
				} else if ("SEQUENCE_EXCLUDED_REGION" == toks[0]) {
					auto secondToks = tokenizeString(toks[1], " ");
					for (const auto & tok : secondToks) {
						auto regToks = tokenizeString(tok, ",");
						if (2 != regToks.size()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__
									<< ", error position values should be two values separated by a \",\" , found "
									<< toks.size() << " objects instead for line " << line
									<< " for tok: " << toks[1] << " for sub toks: " << tok
									<< "\n";
							throw std::runtime_error { ss.str() };
						}
						currentResult->sequence_excluded_region_.emplace_back(
								Primer3Runner::region { njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[0]), njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[1]) });
					}
				} else if ("PRIMER_LEFT_NUM_RETURNED" == toks[0]) {
					currentResult->primer_left_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_RIGHT_NUM_RETURNED" == toks[0]) {
					currentResult->primer_right_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_INTERNAL_NUM_RETURNED" == toks[0]) {
					currentResult->primer_internal_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_PAIR_NUM_RETURNED" == toks[0]) {
					currentResult->primer_pair_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else {
					if (ignoreUnexpected) {
						currentResult->misc_.emplace(toks[0], toks[1]);
					} else {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", unhandled case for line: " << line << ", case: " << toks[0] << "\n";
						throw std::runtime_error { ss.str() };
					}
				}
			}
		}
	}
	return results;
}





std::vector<std::shared_ptr<Primer3Runner::Primer3ResultsGeneric>> Primer3Runner::Primer3ResultsGeneric::parsePrimer3OutputResults(
		const bfs::path & input, bool ignoreUnexpected) {
	std::vector<std::shared_ptr<Primer3Runner::Primer3ResultsGeneric>> results;

	InputStream in(input);
	std::string line = "";
	std::shared_ptr<Primer3ResultsGeneric> currentResult = std::make_shared<Primer3Runner::Primer3ResultsGeneric>();
	std::shared_ptr<PrimerPairGeneric> currentPair = std::make_shared<PrimerPairGeneric>();
	currentResult->primerPairs_.emplace_back(currentPair);
	uint32_t currentPrimerPairID = 0;



	std::regex primerInfoPat { "^(PRIMER_(LEFT|RIGHT|PAIR)_([0-9]+))_?(.*)=.*" };
	while (njh::files::crossPlatGetline(in, line)) {
		if ("=" == line) {
			results.emplace_back(currentResult);
			currentResult = std::make_shared<Primer3ResultsGeneric>();
			currentPair = std::make_shared<PrimerPairGeneric>();
			currentResult->primerPairs_.emplace_back(currentPair);
			currentPrimerPairID = 0;
		} else {
			auto toks = tokenizeString(line, "=");
			if (1 == toks.size() && '-' != line.back()) {
				//empty value
				continue;
			}
			if (2 != toks.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__
						<< ", error line should be key=value, didn't find two values separated by a =, found "
						<< toks.size() << " objects instead for line " << line << "\n";
				throw std::runtime_error { ss.str() };
			}
			std::smatch match;
			if (std::regex_match(line, match, primerInfoPat)) {
				//matched primer
				std::string name = match[1];
				std::string leftOrRightOrPair = match[2];
				std::string pairName = match[3];

				uint32_t pairNameNum = njh::StrToNumConverter::stoToNum<uint32_t>(pairName);
				if(pairNameNum != currentPrimerPairID){
					currentPair = std::make_shared<PrimerPairGeneric>();
					currentPrimerPairID = pairNameNum;
					currentResult->primerPairs_.emplace_back(currentPair);
				}

				if (leftOrRightOrPair != "LEFT" && leftOrRightOrPair != "RIGHT" && leftOrRightOrPair != "PAIR") {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", unexpected designation for line: " << line
							<< ", was expecting RIGHT or LEFT not: " << match[2] << "\n";
					throw std::runtime_error { ss.str() };
				}
				bool right = leftOrRightOrPair == "RIGHT";

				if ("" == match[4]) {
					//position
					auto positionToks = tokenizeString(toks[1], ",");
					if (2 != positionToks.size()) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__
								<< ", error position values should be two values separated by a \",\", found "
								<< toks.size() << " objects instead for line " << line
								<< " for tok: " << toks[1] << "\n";
						throw std::runtime_error { ss.str() };
					}
					auto start = njh::StrToNumConverter::stoToNum<uint32_t>(positionToks[0]);
					auto size = njh::StrToNumConverter::stoToNum<uint32_t>(positionToks[1]);
					if (right) {
						currentPair->right_.originalPos_.start_ = start;
						currentPair->right_.originalPos_.size_ = size;
						currentPair->right_.forwardOrientationPos_.start_ = start+ 1 - size;
						currentPair->right_.forwardOrientationPos_.size_ = size;
						currentPair->right_.name_ = name;
					}else{
						currentPair->left_.originalPos_.start_ = start;
						currentPair->left_.originalPos_.size_ = size;
						currentPair->left_.forwardOrientationPos_.start_ = start;
						currentPair->left_.forwardOrientationPos_.size_ = size;
						currentPair->left_.name_ = name;
					}
				} else if ("PENALTY" == match[4]) {
					if (leftOrRightOrPair == "PAIR") {
						currentPair->penalty_= njh::StrToNumConverter::stoToNum<double>(toks[1]);
					} else if (leftOrRightOrPair == "RIGHT") {
						currentPair->right_.penalty_= njh::StrToNumConverter::stoToNum<double>(toks[1]);
					} else {
						currentPair->left_.penalty_= njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("SEQUENCE" == match[4]) {
					if(right){
						currentPair->right_.seq_ = toks[1];
					}else{
						currentPair->left_.seq_ = toks[1];
					}
				}  else if ("PROBLEMS" == match[4]) {
					if(right){
						currentPair->right_.problems_ = njh::tokenizeString(toks[1], ";");
					}else{
						currentPair->left_.problems_ = njh::tokenizeString(toks[1], ";");
					}
				} else if ("TM" == match[4]) {
					if(right){
						currentPair->right_.tm_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentPair->left_.tm_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("GC_PERCENT" == match[4]) {
					if(right){
						currentPair->right_.gc_percent_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentPair->left_.gc_percent_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("SELF_ANY_TH" == match[4]) {
					if(right){
						currentPair->right_.self_any_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentPair->left_.self_any_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}

				} else if ("SELF_END_TH" == match[4]) {
					if(right){
						currentPair->right_.self_end_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentPair->left_.self_end_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("HAIRPIN_TH" == match[4]) {
					if(right){
						currentPair->right_.hairpin_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentPair->left_.hairpin_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("END_STABILITY" == match[4]) {
					if(right){
						currentPair->right_.end_stability_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentPair->left_.end_stability_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
				} else if ("COMPL_ANY_TH" == match[4]) {
					currentPair->compl_any_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("COMPL_END_TH" == match[4]) {
					currentPair->compl_end_th_ =
							njh::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("PRODUCT_SIZE" == match[4]) {
					currentPair->product_size_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
					currentPair->name_ = name;
				} else {
					if(ignoreUnexpected){
						currentResult->misc_.emplace(toks[0], toks[1]);
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", unexpected sub case for : " << line
								<< ", case: " << match[4] << "\n";
						throw std::runtime_error { ss.str() };
					}
				}
			} else {
				/*
				 * 	std::string ;
	std::string ;
	std::string ;
				 */
				if ("SEQUENCE_ID" == toks[0]) {
					currentResult->sequence_id_ = toks[1];
				} else if ("SEQUENCE_TEMPLATE" == toks[0]) {
					currentResult->sequence_template_ = toks[1];
				} else if ("PRIMER_WARNING" == toks[0]) {
					currentResult->warnings_ = tokenizeString(toks[1], ";"); //=pick_sequencing_primers
				} else if ("PRIMER_TASK" == toks[0]) {
					currentResult->primer_task_ = toks[1];
				} else if ("PRIMER_LEFT_EXPLAIN" == toks[0]) {
					currentResult->primer_left_explain_ = toks[1];
				} else if ("PRIMER_RIGHT_EXPLAIN" == toks[0]) {
					currentResult->primer_right_explain_ = toks[1];
				} else if ("PRIMER_PAIR_EXPLAIN" == toks[0]) {
					currentResult->primer_pair_explain_ = toks[1];
				} else if ("SEQUENCE_TARGET" == toks[0]) {
					auto secondToks = tokenizeString(toks[1], " ");
					for (const auto & tok : secondToks) {
						auto regToks = tokenizeString(tok, ",");
						if (2 != regToks.size()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__
									<< ", error position values should be two values separated by a \",\" , found "
									<< toks.size() << " objects instead for line " << line
									<< " for tok: " << toks[1] << " for sub toks: " << tok
									<< "\n";
							throw std::runtime_error { ss.str() };
						}
						currentResult->sequence_target_.emplace_back(
								Primer3Runner::region { njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[0]), njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[1]) });
					}
				} else if ("SEQUENCE_EXCLUDED_REGION" == toks[0]) {
					auto secondToks = tokenizeString(toks[1], " ");
					for (const auto & tok : secondToks) {
						auto regToks = tokenizeString(tok, ",");
						if (2 != regToks.size()) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__
									<< ", error position values should be two values separated by a \",\" , found "
									<< toks.size() << " objects instead for line " << line
									<< " for tok: " << toks[1] << " for sub toks: " << tok
									<< "\n";
							throw std::runtime_error { ss.str() };
						}
						currentResult->sequence_excluded_region_.emplace_back(
								Primer3Runner::region { njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[0]), njh::StrToNumConverter::stoToNum<
										uint32_t>(regToks[1]) });
					}
				} else if ("PRIMER_LEFT_NUM_RETURNED" == toks[0]) {
					currentResult->primer_left_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_RIGHT_NUM_RETURNED" == toks[0]) {
					currentResult->primer_right_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_INTERNAL_NUM_RETURNED" == toks[0]) {
					currentResult->primer_internal_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_PAIR_NUM_RETURNED" == toks[0]) {
					currentResult->primer_pair_num_returned_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else {
					if (ignoreUnexpected) {
						currentResult->misc_.emplace(toks[0], toks[1]);
					} else {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", unhandled case for line: " << line << ", case: " << toks[0] << "\n";
						throw std::runtime_error { ss.str() };
					}
				}
			}
		}
	}

	for( auto & res : results){
		if(res->primer_task_ == "pick_sequencing_primers"){
			for(auto & primerPair : res->primerPairs_){
				primerPair->penalty_ = primerPair->left_.penalty_ + primerPair->right_.penalty_;
				primerPair->name_ = njh::replaceString(primerPair->left_.name_, "LEFT", "PAIR");
				uint32_t targetStart = primerPair->left_.forwardOrientationPos_.start_;
				uint32_t targetEnd = primerPair->right_.forwardOrientationPos_.start_ + primerPair->right_.forwardOrientationPos_.size_;
				primerPair->product_size_ = 1 + targetEnd - targetStart;
				auto misMapFind = res->misc_.find("PRIMER_FIRST_BASE_INDEX");
				if(misMapFind != res->misc_.end() && "0" == misMapFind->second){
					primerPair->product_size_ = targetEnd - targetStart;
				}
			}
		}
	}
	return results;
}
}  // namespace njhseq


