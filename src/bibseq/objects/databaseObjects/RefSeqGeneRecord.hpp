#pragma once
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
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"

namespace bibseq {
class RefSeqGeneRecord{
public:

	RefSeqGeneRecord(const std::string & line);
	//members
	//see https://genome.ucsc.edu/cgi-bin/hgTables for description
	uint32_t bin_;
	std::string name_;
	std::string chrom_;
	char strand_; //should be - or +
	uint32_t txStart_;
	uint32_t txEnd_;
	uint32_t cdsStart_;
	uint32_t cdsEnd_;
	uint32_t exonCount_;
	std::vector<uint32_t> exonStarts_;
	std::vector<uint32_t> exonEnds_;
	int32_t score_;
	std::string name2_; //alternate name
	enum class completeness : char{
		none='n',
		unk ='u',
		incmpl = 'i',
		cmpl = 'c'
	};
	completeness cdsStartStat_;
	completeness cdsEndStat_;
	static completeness parseComplStr(const std::string & str){
		completeness ret;
		if(str == "none"){
			ret = completeness::none;
		} else if(str == "unk"){
			ret = completeness::unk;
		} else if(str == "incmpl"){
			ret = completeness::incmpl;
		} else if(str == "cmpl"){
			ret = completeness::cmpl;
		} else {
			std::stringstream ss;
			ss << "Error in processing completeness str: "  << str << "\n";
			ss << "be none, unk, incmp, or cmpl not " << str << "\n";
			throw std::runtime_error{bib::bashCT::boldRed(ss.str())};
		}
		return ret;
	}

	static std::string complToStr(const completeness & comInfo){
		std::string ret = "";
		switch (comInfo) {
			case completeness::cmpl:
				ret = "cmpl";
				break;
			case completeness::incmpl:
				ret = "incmpl";
				break;
			case completeness::none:
				ret = "none";
				break;
			case completeness::unk:
				ret = "unk";
				break;
			default:
				break;
		}
		return ret;
	}
	std::vector<int32_t> exonFrames_;


	std::string toStrLine();
};



std::unordered_map<std::string, std::shared_ptr<RefSeqGeneRecord>> getRefSeqRecs(
		const std::string & refSeqFilename, const VecStr & names,
		const std::string & aliasDictJsonFile = "") ;

} /* namespace bibseq */



