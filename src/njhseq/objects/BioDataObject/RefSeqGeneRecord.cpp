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
/*
 * RefSeqGeneRecord.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "RefSeqGeneRecord.hpp"
#include "njhseq/objects/BioDataObject/BioDataFileIO.hpp"

namespace njhseq {

RefSeqGeneRecord::RefSeqGeneRecord(const std::string & line){
	auto toks = tokenizeString(line, "\t");
	if(toks.size() != 16){
		std::stringstream ss;
		ss << "Error in processing line: "  << line << "\n";
		ss << "should have 16 columns not " << toks.size() << "\n";
		throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
	}
	if(toks[0] == "none"){
		bin_ = 0;
	}else{
		bin_ = estd::stou(toks[0]);
	}
	name_ = toks[1];
	chrom_ = toks[2];
	if(toks[3] == "-"){
		strand_ = '-';
	}else if (toks[3] == "+"){
		strand_ = '+';
	}else{
		std::stringstream ss;
		ss << "Error in processing strand column: "  << toks[3] << "\n";
		ss << "be either + or - not " << toks[3] << "\n";
		throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
	}
	txStart_ = estd::stou(toks[4]);
	txEnd_ = estd::stou(toks[5]);
	if(toks[6] == "none"){
		cdsStart_ = std::numeric_limits<uint32_t>::max();
	}else{
		cdsStart_ = estd::stou(toks[6]);
	}
	if(toks[7] == "none"){
		cdsEnd_ = std::numeric_limits<uint32_t>::max();
	}else{
		cdsEnd_ = estd::stou(toks[7]);
	}
	exonCount_ = estd::stou(toks[8]);
	exonStarts_  = vecStrToVecNum<uint32_t>(tokenizeString(toks[9], ","));
	exonEnds_ = vecStrToVecNum<uint32_t>(tokenizeString(toks[10], ","));
	if(toks[11] == "none"){
		score_ = 0;
	}else{
		score_ = njh::lexical_cast<int32_t>(toks[11]);
	}
	name2_ = toks[12];
	cdsStartStat_ = parseComplStr(toks[13]);
	cdsEndStat_ = parseComplStr(toks[14]);
	if(toks[15] == "none"){
		exonFrames_ = std::vector<int32_t>(exonCount_, -1);
	}else{
		exonFrames_  = vecStrToVecNum<int32_t>(tokenizeString(toks[15], ","));
	}
}


std::string RefSeqGeneRecord::toStrLine() {
	return vectorToString(
			toVecStr(bin_, name_, chrom_, strand_, txStart_, txEnd_, cdsStart_,
					cdsEnd_, exonCount_, vectorToString(exonStarts_, ",") + ",",
					vectorToString(exonEnds_, ",") + ",", score_, name2_,
					complToStr(cdsStartStat_), complToStr(cdsEndStat_),
					vectorToString(exonFrames_, ",") + ","), "\t");
}


std::unordered_map<std::string, std::shared_ptr<RefSeqGeneRecord>> getRefSeqRecs(
		const bfs::path & refSeqFilename, const VecStr & names,
		const bfs::path & aliasDictJsonFile) {
	std::unordered_map<std::string, std::shared_ptr<RefSeqGeneRecord>> ret;



	VecStr codedExtractNames = names;
	std::unordered_map<std::string, std::string> coder;
	//get name alias if any, should be a json dictionary, seq alias to name in refSeqGene file
	if(aliasDictJsonFile != ""){
		codedExtractNames.clear();
		Json::Value root = njh::json::parseFile(aliasDictJsonFile.string());
		for (const auto & n : names) {
			coder[root[n].asString()] = n;
			codedExtractNames.emplace_back(root[n].asString());
		}
	}
	BioDataFileIO <RefSeqGeneRecord> reader((IoOptions(InOptions(refSeqFilename))));
	reader.openIn();
	auto currentRecord = reader.readNextRecord();
	while (nullptr != currentRecord) {
		if(names.empty()){
			ret.emplace(currentRecord->name_, std::move(currentRecord));
		}else{
			if (njh::in(currentRecord->name_, codedExtractNames)) {
				//currentRecord->name2_ = coder[currentRecord->name_];
				ret.emplace(currentRecord->name_, std::move(currentRecord));
			}
		}
		currentRecord = reader.readNextRecord();
	}
	//add alias name if there are any
	if(aliasDictJsonFile != ""){
		for(const auto & code : coder){
			auto search = ret.find(code.first);
			if(search != ret.end()){
				search->second->name2_ = code.second;
			}
		}
	}
	return ret;
}


RefSeqGeneRecord::completeness RefSeqGeneRecord::parseComplStr(const std::string & str){
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
		throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
	}
	return ret;
}

std::string RefSeqGeneRecord::complToStr(const completeness & comInfo){
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

} /* namespace njhseq */
