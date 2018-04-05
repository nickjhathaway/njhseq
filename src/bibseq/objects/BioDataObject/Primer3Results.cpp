
/*
 * Primer3Results.cpp
 *
 *  Created on: Aug 16, 2017
 *      Author: nick
 */
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


#include "Primer3Results.hpp"
#include "bibseq/IO/InputStream.hpp"
#include "bibseq/helpers/seqUtil.hpp"

namespace bibseq {

Json::Value Primer3Results::region::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::typeStr(*this));
	ret["start_"] =bib::json::toJson(start_);
	ret["size_"] =bib::json::toJson(size_);
	return ret;
}


GenomicRegion Primer3Results::Primer::genRegion(const std::string & templateId) const {
	return GenomicRegion(name_, templateId, forwardOrientationPos_.start_,
			forwardOrientationPos_.start_ + forwardOrientationPos_.size_, right_);
}


Json::Value Primer3Results::Primer::toJson() const {
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::typeStr(*this));
	ret["name_"] = bib::json::toJson(name_);
	ret["penality_"] = bib::json::toJson(penality_);
	ret["originalSeq_"] = bib::json::toJson(originalSeq_);
	ret["fowardOrientationSeq_"] = bib::json::toJson(fowardOrientationSeq_);
	ret["tm_"] = bib::json::toJson(tm_);
	ret["gc_percent_"] = bib::json::toJson(gc_percent_);
	ret["self_any_th_"] = bib::json::toJson(self_any_th_);
	ret["self_end_th_"] = bib::json::toJson(self_end_th_);
	ret["hairpin_th_"] = bib::json::toJson(hairpin_th_);
	ret["end_stability_"] = bib::json::toJson(end_stability_);
	ret["right_"] = bib::json::toJson(right_);
	ret["originalPos_"] = bib::json::toJson(originalPos_);
	ret["forwardOrientationPos_"] = bib::json::toJson(forwardOrientationPos_);
	return ret;
}

Json::Value Primer3Results::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::typeStr(*this));
	ret["sequence_id_"] =bib::json::toJson(sequence_id_);
	ret["sequence_template_"] =bib::json::toJson(sequence_template_);
	ret["sequence_target_"] =bib::json::toJson(sequence_target_);
	ret["sequence_excluded_region_"] =bib::json::toJson(sequence_excluded_region_);
	ret["primer_left_num_returned_"] =bib::json::toJson(primer_left_num_returned_);
	ret["primer_right_num_returned_"] =bib::json::toJson(primer_right_num_returned_);
	ret["primer_internal_num_returned_"] =bib::json::toJson(primer_internal_num_returned_);
	ret["primer_pair_num_returned_"] =bib::json::toJson(primer_pair_num_returned_);

	ret["primers_"] =bib::json::toJson(primers_);
	return ret;
}

std::unordered_map<std::string, GenomicRegion> Primer3Results::genPrimersRegions() const{
	std::unordered_map<std::string, GenomicRegion> ret;
	for(const auto & p : primers_){
		ret[p.first] = p.second->genRegion(sequence_id_);
	}
	return ret;
}

std::vector<std::shared_ptr<Primer3Results>> Primer3Results::parsePrimer3OutputResults(
		const bfs::path & input) {
	std::vector<std::shared_ptr<Primer3Results>> results;

	InputStream in(input);
	std::string line = "";
	std::shared_ptr<Primer3Results> currentResult = std::make_shared<
			Primer3Results>();
	std::regex primerInfoPat { "^(PRIMER_(LEFT|RIGHT)_[0-9]+)_?(.*)=.*" };
	while (bib::files::crossPlatGetline(in, line)) {
		if ("=" == line) {
			results.emplace_back(currentResult);
			currentResult = std::make_shared<Primer3Results>();
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
				std::string name = match[1];
				std::string leftOrRight = match[2];
				if (leftOrRight != "LEFT" && leftOrRight != "RIGHT") {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", unhandled case for line: " << line
							<< ", was expecting RIGHT or LEFT not: " << match[2] << "\n";
					throw std::runtime_error { ss.str() };
				}
				bool right = leftOrRight == "RIGHT";
				if (!bib::in(name, currentResult->primers_)) {
					currentResult->primers_[name] = std::make_shared<Primer>();
					currentResult->primers_[name]->right_ = right;
					currentResult->primers_[name]->name_ = name;
				}
				if ("" == match[3]) {
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
					auto start = bib::StrToNumConverter::stoToNum<uint32_t>(
							positionToks[0]);
					auto size = bib::StrToNumConverter::stoToNum<uint32_t>(
							positionToks[1]);
					currentResult->primers_[name]->originalPos_.start_ = start;
					currentResult->primers_[name]->originalPos_.size_ = size;
					currentResult->primers_[name]->forwardOrientationPos_ =
							currentResult->primers_[name]->originalPos_;
					if (right) {
						currentResult->primers_[name]->forwardOrientationPos_.start_ = start
								+ 1 - size;
					}
				} else if ("PENALTY" == match[3]) {
					currentResult->primers_[name]->penality_ =
							bib::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("SEQUENCE" == match[3]) {
					currentResult->primers_[name]->originalSeq_ = toks[1];
					if (right) {
						currentResult->primers_[name]->fowardOrientationSeq_ =
								seqUtil::reverseComplement(
										currentResult->primers_[name]->originalSeq_, "DNA");
					} else {
						currentResult->primers_[name]->fowardOrientationSeq_ =
								currentResult->primers_[name]->originalSeq_;
					}
				} else if ("TM" == match[3]) {
					currentResult->primers_[name]->tm_ = bib::StrToNumConverter::stoToNum<
							double>(toks[1]);
				} else if ("GC_PERCENT" == match[3]) {
					currentResult->primers_[name]->gc_percent_ =
							bib::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("SELF_ANY_TH" == match[3]) {
					currentResult->primers_[name]->self_any_th_ =
							bib::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("SELF_END_TH" == match[3]) {
					currentResult->primers_[name]->self_end_th_ =
							bib::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("HAIRPIN_TH" == match[3]) {
					currentResult->primers_[name]->hairpin_th_ =
							bib::StrToNumConverter::stoToNum<double>(toks[1]);
				} else if ("END_STABILITY" == match[3]) {
					currentResult->primers_[name]->end_stability_ =
							bib::StrToNumConverter::stoToNum<double>(toks[1]);
				} else {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", unhandled case for line: " << line
							<< ", case: " << match[3] << "\n";
					throw std::runtime_error { ss.str() };
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
								Primer3Results::region { bib::StrToNumConverter::stoToNum<
										uint32_t>(regToks[0]), bib::StrToNumConverter::stoToNum<
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
								Primer3Results::region { bib::StrToNumConverter::stoToNum<
										uint32_t>(regToks[0]), bib::StrToNumConverter::stoToNum<
										uint32_t>(regToks[1]) });
					}
				} else if ("PRIMER_LEFT_NUM_RETURNED" == toks[0]) {
					currentResult->primer_left_num_returned_ =
							bib::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_RIGHT_NUM_RETURNED" == toks[0]) {
					currentResult->primer_right_num_returned_ =
							bib::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_INTERNAL_NUM_RETURNED" == toks[0]) {
					currentResult->primer_internal_num_returned_ =
							bib::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else if ("PRIMER_PAIR_NUM_RETURNED" == toks[0]) {
					currentResult->primer_pair_num_returned_ =
							bib::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
				} else {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", unhandled case for line: " << line
							<< ", case: " << match[3] << "\n";
					throw std::runtime_error { ss.str() };
				}
			}
		}
	}
	return results;
}
}  // namespace bibseq


