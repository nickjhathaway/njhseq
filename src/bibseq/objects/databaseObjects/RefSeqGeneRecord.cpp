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
 * RefSeqGeneRecord.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "RefSeqGeneRecord.hpp"
#include "bibseq/objects/databaseObjects/DataFileReader.hpp"
namespace bibseq {

RefSeqGeneRecord::RefSeqGeneRecord(const std::string & line){
	auto toks = tokenizeString(line, "\t");
	if(toks.size() != 16){
		std::stringstream ss;
		ss << "Error in processing line: "  << line << "\n";
		ss << "should have 16 columns not " << toks.size() << "\n";
		throw std::runtime_error{bib::bashCT::boldRed(ss.str())};
	}
	if(toks[0] == "none"){
		bin_ = 0;
	}else{
		bin_ = bib::lexical_cast<uint32_t>(toks[0]);
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
		throw std::runtime_error{bib::bashCT::boldRed(ss.str())};
	}
	txStart_ = bib::lexical_cast<uint32_t>(toks[4]);
	txEnd_ = bib::lexical_cast<uint32_t>(toks[5]);
	if(toks[6] == "none"){
		cdsStart_ = std::numeric_limits<uint32_t>::max();
	}else{
		cdsStart_ = bib::lexical_cast<uint32_t>(toks[6]);
	}
	if(toks[7] == "none"){
		cdsEnd_ = std::numeric_limits<uint32_t>::max();
	}else{
		cdsEnd_ = bib::lexical_cast<uint32_t>(toks[7]);
	}
	exonCount_ = bib::lexical_cast<uint32_t>(toks[8]);
	exonStarts_  = vecStrToVecNum<uint32_t>(tokenizeString(toks[9], ","));
	exonEnds_ = vecStrToVecNum<uint32_t>(tokenizeString(toks[10], ","));
	if(toks[11] == "none"){
		score_ = 0;
	}else{
		score_ = bib::lexical_cast<int32_t>(toks[11]);
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
		const std::string & refSeqFilename, const VecStr & names,
		const std::string & aliasDictJsonFile) {
	std::unordered_map<std::string, std::shared_ptr<RefSeqGeneRecord>> ret;
	std::ifstream inputFile(refSeqFilename);
	if (!inputFile) {
		std::stringstream ss;
		ss << bib::bashCT::boldRed("Error in opening " + refSeqFilename)
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}

	std::string line = "";
	VecStr codedExtractNames = names;
	std::unordered_map<std::string, std::string> coder;
	//get name alias if any, should be a json dictionary, seq alias to name in refSeqGene file
	if(aliasDictJsonFile != ""){
		codedExtractNames.clear();
		Json::Value root;
		Json::Reader jReader;
		std::ifstream aliasFile(aliasDictJsonFile);
		jReader.parse(aliasFile, root);
		for (const auto & n : names) {
			coder[root[n].asString()] = n;
			codedExtractNames.emplace_back(root[n].asString());
		}
	}

	while (std::getline(inputFile, line)) {
		//skip comment lines
		if(beginsWith(line, "#") || line == ""){
			continue;
		}
		try {
			auto currentRecord = std::make_shared<RefSeqGeneRecord>(line);
			if(names.empty()){
				ret.emplace(currentRecord->name_, std::move(currentRecord));
			}else{
				if (bib::in(currentRecord->name_, codedExtractNames)) {
					//currentRecord->name2_ = coder[currentRecord->name_];
					ret.emplace(currentRecord->name_, std::move(currentRecord));
				}
			}
		} catch (std::exception & e) {
			std::stringstream ss;
			ss << e.what() << "\n";
			ss << "on line: " << "\n" << line << "\n";
			throw std::runtime_error { ss.str() };
		}
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

} /* namespace bibseq */
