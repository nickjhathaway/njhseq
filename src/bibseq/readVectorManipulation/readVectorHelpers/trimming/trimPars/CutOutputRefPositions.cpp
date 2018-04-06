/*
 * CutOutputRefPositions.cpp
 *
 *  Created on: Dec 13, 2017
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
#include "CutOutputRefPositions.hpp"

namespace bibseq {

Json::Value CutOutputRefPositions::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["refStart_"] = bib::json::toJson(refStart_);
	ret["refStartLength_"] = bib::json::toJson(refStartLength_);
	ret["refStop_"] = bib::json::toJson(refStop_);
	ret["refStopLength_"] = bib::json::toJson(refStopLength_);
	ret["refStartStr_"] = bib::json::toJson(bib::pasteAsStr(refStart_, ":", refStartLength_));
	ret["refStopStr_"] = bib::json::toJson(bib::pasteAsStr(refStop_, ":", refStopLength_));
	ret["checkComp_"] = bib::json::toJson(checkComp_);
	ret["comp_"] = bib::json::toJson(comp_);
	return ret;
}

std::vector<CutOutputRefPositions> CutOutputRefPositions::readInPositions(
		const bfs::path & positionsJsonFnp) {
	std::vector<CutOutputRefPositions> ret;

	auto refPositionsJson = bib::json::parseFile(positionsJsonFnp.string());
	if (!refPositionsJson.isArray()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, expected " << positionsJsonFnp
				<< " to contain a json array" << "\n";
		throw std::runtime_error { ss.str() };
	}
	for (const auto & val : refPositionsJson) {
		CutOutputRefPositions currentPos;
		if (val.isMember("name_")) {
			currentPos.name_ = val["name_"].asString();
		}
		if (val.isMember("refStart_")) {
			currentPos.refStart_ = val["refStart_"].asUInt();
		}
		if (val.isMember("refStartLength_")) {
			currentPos.refStartLength_ = val["refStartLength_"].asUInt();
		}
		if (val.isMember("refStop_")) {
			currentPos.refStop_ = val["refStop_"].asUInt();
		}
		if (val.isMember("refStopLength_")) {
			currentPos.refStopLength_ = val["refStopLength_"].asUInt();
		}
		if (val.isMember("checkComp_")) {
			currentPos.checkComp_ = val["checkComp_"].asBool();
		}
		if (val.isMember("comp_")) {
			if (val["comp_"].isMember("oneBaseIndel_")) {
				currentPos.comp_.oneBaseIndel_ =
						val["comp_"]["oneBaseIndel_"].asDouble();
			}
			if (val["comp_"].isMember("twoBaseIndel_")) {
				currentPos.comp_.twoBaseIndel_ =
						val["comp_"]["twoBaseIndel_"].asDouble();
			}
			if (val["comp_"].isMember("largeBaseIndel_")) {
				currentPos.comp_.largeBaseIndel_ =
						val["comp_"]["largeBaseIndel_"].asDouble();
			}
			if (val["comp_"].isMember("hqMismatches_")) {
				currentPos.comp_.hqMismatches_ =
						val["comp_"]["hqMismatches_"].asDouble();
			}
		}
		ret.emplace_back(currentPos);
	}
	for(auto & pos : ret){
		 setPositions(pos, __PRETTY_FUNCTION__);
	}
	return ret;
}


void CutOutputRefPositions::checkPositionsThrow(const std::string & funcName) const{
	if (refStopLength_ == 0) {
		std::stringstream ss;
		ss << funcName << ", error refStopLength_ can't be zero " << '\n';
		throw std::runtime_error { ss.str() };
	}
	if (refStartLength_ == 0) {
		std::stringstream ss;
		ss << funcName << ", error refStartLength_ can't be zero " << '\n';
		throw std::runtime_error { ss.str() };
	}
	if (refStopLength_ - 1 > refStop_) {
		std::stringstream ss;
		ss << funcName << ", error refStopLength_ is greater than stop "
				<< refStopLength_ - 1 << ", stop: " << refStop_
				<< '\n';
		throw std::runtime_error { ss.str() };
	}
	if (refStart_ > refStop_) {
		std::stringstream ss;
		ss << funcName << ", error start is greater than stop "
				<< refStart_ << ", stop: " << refStop_
				<< '\n';
		throw std::runtime_error { ss.str() };
	}

	if (refStart_ + refStartLength_ - 1
			> refStop_) {
		std::stringstream ss;
		ss << funcName
				<< ", error refStart_ + refStartLength_ - 1, "
				<< refStart_ + refStartLength_ - 1
				<< ", is greater than refStop_ + refStopLength_ - 1, "
				<< refStop_ << '\n';
		throw std::runtime_error { ss.str() };
	}
}

void CutOutputRefPositions::checkPositionsThrow(const std::string & funcName,
		const size_t refLength,
		const std::string & refName) const {
	checkPositionsThrow(funcName);
	std::stringstream errorStream;
	bool failPositionChecks = false;
	if (refStart_ >= refLength) {
		errorStream << "Ref Start, " << refStart_
				<< ", is greater than or equal to ref length for " << refName << ", length: "
				<< refLength << "\n";
		failPositionChecks = true;
	}

	if (refStop_ >= refLength) {
		errorStream << "Ref Stop, " << refStop_
				<< ", is greater than or equal to ref length for " << refName << ", length: "
				<< refLength << "\n";
		failPositionChecks = true;
	}

	if(failPositionChecks){
		throw std::runtime_error{errorStream.str()};
	}
}

void CutOutputRefPositions::setPositions(CutOutputRefPositions & refPositions,
		std::string funcName) {
	if("" != refPositions.refStopStr_){
		if (std::string::npos != refPositions.refStopStr_.find(":")) {
			auto toks = tokenizeString(refPositions.refStopStr_, ":");
			if (2 != toks.size()) {
				std::stringstream ss;
				ss << funcName
						<< ", error, should be two numbers separated by a colon, not "
						<< refPositions.refStopStr_ << ", for --refStop" << "\n";
			}
			refPositions.refStop_ = bib::StrToNumConverter::stoToNum<uint32_t>(
					toks[0]);
			refPositions.refStopLength_ = bib::StrToNumConverter::stoToNum<uint32_t>(
					toks[1]);
		} else {
			refPositions.refStop_ = bib::StrToNumConverter::stoToNum<uint32_t>(
					refPositions.refStopStr_);
		}
	}
	if("" != refPositions.refStartStr_){
		if (std::string::npos != refPositions.refStartStr_.find(":")) {
			auto toks = tokenizeString(refPositions.refStartStr_, ":");
			if (2 != toks.size()) {
				std::stringstream ss;
				ss << funcName
						<< ", error, should be two numbers separated by a colon, not "
						<< refPositions.refStartStr_ << ", for --refStart" << "\n";
			}
			refPositions.refStart_ = bib::StrToNumConverter::stoToNum<uint32_t>(
					toks[0]);
			refPositions.refStartLength_ = bib::StrToNumConverter::stoToNum<uint32_t>(
					toks[1]);
		} else {
			refPositions.refStart_ = bib::StrToNumConverter::stoToNum<uint32_t>(
					refPositions.refStartStr_);
		}
	}
	refPositions.checkPositionsThrow(__PRETTY_FUNCTION__);
}

std::string CutOutputRefPositions::getId() const{
	return bib::pasteAsStr(refStart_, ":", refStartLength_, "-", refStop_, ":", refStopLength_);
}

}  // namespace bibseq



