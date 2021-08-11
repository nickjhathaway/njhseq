
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

Json::Value Primer3Results::region::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["start_"] =njh::json::toJson(start_);
	ret["size_"] =njh::json::toJson(size_);
	return ret;
}


GenomicRegion Primer3Results::Primer::genRegion(const std::string & templateId) const {
	return GenomicRegion(name_, templateId, forwardOrientationPos_.start_,
			forwardOrientationPos_.start_ + forwardOrientationPos_.size_, right_);
}


Json::Value Primer3Results::PrimerPair::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::typeStr(*this));
	ret["name_"] = njh::json::toJson(name_);
	ret["left_"] = njh::json::toJson(left_);
	ret["right_"] = njh::json::toJson(right_);
	ret["compl_any_th_"] = njh::json::toJson(compl_any_th_);
	ret["compl_end_th_"] = njh::json::toJson(compl_end_th_);
	ret["penalty_"] = njh::json::toJson(penalty_);
	ret["product_size_"] = njh::json::toJson(product_size_);
	return ret;
}


Json::Value Primer3Results::Primer::toJson() const {
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

Json::Value Primer3Results::toJson() const{
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

	ret["primers_"] =njh::json::toJson(primerPairs_);
	return ret;
}

std::unordered_map<std::string, GenomicRegion> Primer3Results::genPrimersRegions() const {
	std::unordered_map<std::string, GenomicRegion> ret;
	for (const auto &p : primerPairs_) {
		ret[p->right_.name_] = p->right_.genRegion(sequence_id_);
		ret[p->left_.name_] = p->left_.genRegion(sequence_id_);
	}
	return ret;
}

std::vector<std::shared_ptr<Primer3Results>> Primer3Results::parsePrimer3OutputResults(
		const bfs::path & input, bool ignoreUnexpected) {
	std::vector<std::shared_ptr<Primer3Results>> results;

	InputStream in(input);
	std::string line = "";
	std::unordered_multimap<std::string, std::string> misc;
	std::shared_ptr<Primer3Results> currentResult = std::make_shared<Primer3Results>();
	std::shared_ptr<PrimerPair> currentPair = std::make_shared<PrimerPair>();
	currentResult->primerPairs_.emplace_back(currentPair);
	uint32_t currentPrimerPairID = 0;



	std::regex primerInfoPat { "^(PRIMER_(LEFT|RIGHT|PAIR)_([0-9]+))_?(.*)=.*" };
	while (njh::files::crossPlatGetline(in, line)) {
		if ("=" == line) {
			results.emplace_back(currentResult);
			currentResult = std::make_shared<Primer3Results>();
			currentPair = std::make_shared<PrimerPair>();
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
					currentPair = std::make_shared<PrimerPair>();
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
					if(right){
						currentPair->right_.self_any_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}else{
						currentPair->left_.self_any_th_ = njh::StrToNumConverter::stoToNum<double>(toks[1]);
					}
					currentPair->compl_any_th_ =
							njh::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("COMPL_END_TH" == match[4]) {
					currentPair->compl_end_th_ =
							njh::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("PRODUCT_SIZE" == match[4]) {
					currentPair->product_size_ =
							njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
					currentPair->name_ = name;
				} else {
					if(ignoreUnexpected){
						misc.emplace(toks[0], toks[1]);
					}else{
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", unexpected sub case for : " << line
								<< ", case: " << match[4] << "\n";
						throw std::runtime_error { ss.str() };
					}
				}
			} else {
				if ("SEQUENCE_ID" == toks[0]) {
					currentResult->sequence_id_ = toks[1];
				} else if ("SEQUENCE_TEMPLATE" == toks[0]) {
					currentResult->sequence_template_ = toks[1];
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
								Primer3Results::region { njh::StrToNumConverter::stoToNum<
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
								Primer3Results::region { njh::StrToNumConverter::stoToNum<
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
						misc.emplace(toks[0], toks[1]);
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
}  // namespace njhseq


