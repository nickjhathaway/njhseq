/*
 * BLASTHitTabular.cpp
 *
 *  Created on: Apr 20, 2018
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
#include "BLASTHitTabular.hpp"

namespace njhseq {

BLASTHitTab::BLASTHitTab():
		queryName_(""),
		subjectName_(""),
		perId_(std::numeric_limits<double>::max()),
		alignLen_(std::numeric_limits<uint32_t>::max()),
		mismatches_(std::numeric_limits<uint32_t>::max()),
		gapOpens_(std::numeric_limits<uint32_t>::max()),
		qStart_(std::numeric_limits<uint32_t>::max()),
		qEnd_(std::numeric_limits<uint32_t>::max()),
		sStart_(std::numeric_limits<uint32_t>::max()),
		sEnd_(std::numeric_limits<uint32_t>::max()),
		evalue_(std::numeric_limits<uint32_t>::max()),
		bitScore_(std::numeric_limits<uint32_t>::max())
		{

}

BLASTHitTab::BLASTHitTab(const std::string & line){

	auto toks = tokenizeString(line, "\t");
	if(12 != toks.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in processing line: " << "\n";
		ss << line << "\n";
		throw std::runtime_error{ss.str()};
	}
	queryName_ = toks[0];
	subjectName_ = toks[1];
	perId_ = njh::StrToNumConverter::stoToNum<double>(toks[2]);
	alignLen_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[3]);
	mismatches_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[4]);
	gapOpens_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[5]);
	qStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[6]);
	qEnd_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[7]);
	sStart_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[8]);
	sEnd_ = njh::StrToNumConverter::stoToNum<uint32_t>(toks[9]);
	evalue_ = njh::StrToNumConverter::stoToNum<double>(toks[10]);
	bitScore_ = njh::StrToNumConverter::stoToNum<double>(toks[11]);
}



Json::Value BLASTHitTab::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["queryName_"] = njh::json::toJson(queryName_);
	ret["subjectName_"] = njh::json::toJson(subjectName_);
	ret["perId_"] = njh::json::toJson(perId_);
	ret["alignLen_"] = njh::json::toJson(alignLen_);
	ret["mismatches_"] = njh::json::toJson(mismatches_);
	ret["gapOpens_"] = njh::json::toJson(gapOpens_);
	ret["qStart_"] = njh::json::toJson(qStart_);
	ret["qEnd_"] = njh::json::toJson(qEnd_);
	ret["sStart_"] = njh::json::toJson(sStart_);
	ret["sEnd_"] = njh::json::toJson(sEnd_);
	ret["evalue_"] = njh::json::toJson(evalue_);
	ret["bitScore_"] = njh::json::toJson(bitScore_);
	return ret;
}


}  // namespace njhseq

