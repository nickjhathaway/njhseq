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
/*
 * BedRecord.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "bibseq/objects/BioDataObject/BedRecordCore.hpp"

namespace bibseq {
BedRecordCore::BedRecordCore(const std::string & line) {
	auto toks = tokenizeString(line, "\t");
	if (toks.size() < 6) {
		std::stringstream ss;
		ss << "Error in parsing line: " << line << "\n";
		ss << "Error in read bed file, need to have at least first 6 fields not: "
				<< toks.size() << std::endl;
		ss << "1)chrom,2)chromStart,3)chromEnd,4)name,5)score,6)strand"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	chrom_ = toks[0];
	chromStart_ = bib::lexical_cast<uint32_t>(toks[1]);
	chromEnd_ = bib::lexical_cast<uint32_t>(toks[2]);
	name_ = toks[3];
	score_ = bib::lexical_cast<double>(toks[4]);
	if (toks[5] == "-") {
		strand_ = '-';
	} else if (toks[5] == "+") {
		strand_ = '+';
	} else {
		std::stringstream ss;
		ss << "Error in parsing strand field: " << toks[5] << "\n";
		ss << "should be either + or - not " << toks[5] << "\n";
		throw std::runtime_error { ss.str() };
	}
}

BedRecordCore::BedRecordCore(std::string chrom, uint32_t chromStart,
		uint32_t chromEnd, std::string name, double score, char strand) :
		chrom_(chrom), chromStart_(chromStart), chromEnd_(chromEnd), name_(name), score_(
				score), strand_(strand) {
	if (strand_ != '+' && strand_ != '-') {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error strand must be + or - not " << strand_
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
}

BedRecordCore::BedRecordCore() :
		chrom_(""), chromStart_(0), chromEnd_(0), name_(""), score_(0), strand_('+') {

}

bool BedRecordCore::reverseStrand() const {
	return '-' == strand_;
}

std::string BedRecordCore::toDelimStr() const {
	return vectorToString(
			toVecStr(chrom_, chromStart_, chromEnd_, name_, score_, strand_), "\t");
}

Json::Value BedRecordCore::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::BedRecordCore");
	ret["chrom_"] = bib::json::toJson(chrom_);
	ret["chromStart_"] = bib::json::toJson(chromStart_);
	ret["chromEnd_"] = bib::json::toJson(chromEnd_);
	ret["name_"] = bib::json::toJson(name_);
	ret["score_"] = bib::json::toJson(score_);
	ret["strand_"] = bib::json::toJson(strand_);
	return ret;
}

} /* namespace bibseq */
