/*
 * BLASTHitTabular.cpp
 *
 *  Created on: Apr 20, 2018
 *      Author: nick
 */
//
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
#include "BLASTHitTabular.hpp"

namespace bibseq {

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
	perId_ = bib::StrToNumConverter::stoToNum<double>(toks[2]);
	alignLen_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[3]);
	mismatches_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[4]);
	gapOpens_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[5]);
	qStart_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[6]);
	qEnd_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[7]);
	sStart_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[8]);
	sEnd_ = bib::StrToNumConverter::stoToNum<uint32_t>(toks[9]);
	evalue_ = bib::StrToNumConverter::stoToNum<double>(toks[10]);
	bitScore_ = bib::StrToNumConverter::stoToNum<double>(toks[11]);
}



Json::Value BLASTHitTab::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["queryName_"] = bib::json::toJson(queryName_);
	ret["subjectName_"] = bib::json::toJson(subjectName_);
	ret["perId_"] = bib::json::toJson(perId_);
	ret["alignLen_"] = bib::json::toJson(alignLen_);
	ret["mismatches_"] = bib::json::toJson(mismatches_);
	ret["gapOpens_"] = bib::json::toJson(gapOpens_);
	ret["qStart_"] = bib::json::toJson(qStart_);
	ret["qEnd_"] = bib::json::toJson(qEnd_);
	ret["sStart_"] = bib::json::toJson(sStart_);
	ret["sEnd_"] = bib::json::toJson(sEnd_);
	ret["evalue_"] = bib::json::toJson(evalue_);
	ret["bitScore_"] = bib::json::toJson(bitScore_);
	return ret;
}


}  // namespace bibseq

