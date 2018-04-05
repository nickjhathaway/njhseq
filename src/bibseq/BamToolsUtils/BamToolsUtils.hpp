#pragma once
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
/*
 * BamToolsUtils.hpp
 *
 *  Created on: May 28, 2015
 *      Author: nick
 */

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>


#include "bibseq/common.h"
#include "bibseq/objects/BioDataObject/GenomicRegion.hpp"
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "bibseq/alignment/alnCache.h"
#include "bibseq/objects/dataContainers/tables/table.hpp"

namespace bibseq {


std::vector<GenomicRegion> genGenRegionsFromRefData(const BamTools::RefVector & rData);

seqInfo bamAlnToSeqInfo(const BamTools::BamAlignment & aln, bool keepPlusStrandOrientation = false);

void bamWriteFasta(const BamTools::BamAlignment & aln, std::ostream & out);
void bamWriteFastq(const BamTools::BamAlignment & aln, std::ostream & out);

alnInfoLocal bamAlnToAlnInfoLocal(const std::vector<BamTools::CigarOp> & cigarData);

std::unordered_map<size_t,alnInfoLocal> bamAlnToAlnInfoLocal(const BamTools::BamAlignment & bAln);


namespace bamAlnFlagCheck {

inline bool IsDuplicate(uint32_t flag) {
	return ((flag & 0x0400) != 0);
}
inline bool IsFailedQC(uint32_t flag) {
	return ((flag & 0x0200) != 0);
}
inline bool IsFirstMate(uint32_t flag) {
	return ((flag & 0x0040) != 0);
}
inline bool IsMapped(uint32_t flag) {
	return ((flag & 0x0004) == 0);
}
inline bool IsMateMapped(uint32_t flag) {
	return ((flag & 0x0008) == 0);
}
inline bool IsMateReverseStrand(uint32_t flag) {
	return ((flag & 0x0020) != 0);
}
inline bool IsPaired(uint32_t flag) {
	return ((flag & 0x0001) != 0);
}
inline bool IsPrimaryAlignment(uint32_t flag) {
	return ((flag & 0x0100) == 0);
}
inline bool IsProperPair(uint32_t flag) {
	return ((flag & 0x0002) != 0);
}
inline bool IsReverseStrand(uint32_t flag) {
	return ((flag & 0x0010) != 0);
}
inline bool IsSecondMate(uint32_t flag) {
	return ((flag & 0x0080) != 0);
}
}

bool IsBamSorted(const std::string & filename, bool verbose = false);


void checkBamOpenThrow(BamTools::BamReader & bReader, const bfs::path & bamFnp);

void loadBamIndexThrow(BamTools::BamReader & bReader);

table refDataVecToTab(const std::vector<BamTools::RefData> & refInfos);
std::vector<BamTools::RefData> tabToRefDataVec(const table & refDataTab);

void logAlnInfo(std::ostream & out, BamTools::RefVector & refInfo,
		const BamTools::BamAlignment & aln, std::string indent = "");



void setBamFileRegionThrow(BamTools::BamReader & bReader, const GenomicRegion & region);

void checkBamFilesForIndexesAndAbilityToOpen(const std::vector<bfs::path> & bamFnps);


} /* namespace bibseq */



namespace bib {
namespace json {
namespace JsonConversion {
template<>
inline Json::Value toJson(const BamTools::CigarOp & cOps){
	Json::Value ret;
	ret["class"] = bib::getTypeName(cOps);
	ret["Length"] = bib::json::toJson(cOps.Length);
	ret["Type"] = bib::json::toJson(cOps.Type);
	return ret;
}


/* Tag Data info
 *
 * A = printable character
 * i = Signed 32-bit integer
 * f = Single-precision floating number
 * Z = Printable string, including space
 * H = byte array in the Hex format
 * B[cCsSiIf] = numeric array with the second letter being the numeric type (see below)
 *
 * numeric types
 *
 * c = int8_t
 * C = uint8_t
 * s = int16_t
 * S = uint16_t
 * i = int32_t
 * I = uint32_t
 * f = float
 *
 */

template<>
inline Json::Value toJson(const BamTools::BamAlignment & bAln){
	Json::Value ret;
	ret["class"] = bib::getTypeName(bAln);
	ret["AlignedBases"] = bib::json::toJson(bAln.AlignedBases);
	ret["AlignmentFlag"] = bib::json::toJson(bAln.AlignmentFlag);
	ret["Bin"] = bib::json::toJson(bAln.Bin);
	ret["CigarData"] = bib::json::toJson(bAln.CigarData);
	ret["Filename"] = bib::json::toJson(bAln.Filename);
	ret["InsertSize"] = bib::json::toJson(bAln.InsertSize);
	ret["Length"] = bib::json::toJson(bAln.Length);
	ret["MapQuality"] = bib::json::toJson(bAln.MapQuality);
	ret["MatePosition"] = bib::json::toJson(bAln.MatePosition);
	ret["MateRefID"] = bib::json::toJson(bAln.MateRefID);
	ret["Name"] = bib::json::toJson(bAln.Name);
	ret["Position"] = bib::json::toJson(bAln.Position);
	ret["Qualities"] = bib::json::toJson(bAln.Qualities);
	ret["QueryBases"] = bib::json::toJson(bAln.QueryBases);
	ret["RefID"] = bib::json::toJson(bAln.RefID);
	ret["TagDataRaw"] = bib::json::toJson(bAln.TagData);
	ret["TagNames"] = bib::json::toJson(bAln.GetTagNames());
	std::unordered_map<std::string, Json::Value> tagData;
	for(const auto & tag : bAln.GetTagNames()){
		//tagData[tag] = bAln.GetTag()
		//std::cout << tag << std::endl;
		char type = ' ';
		bAln.GetTagType(tag, type);
		std::string val;
		if('c' == type || 'C' == type){
			int32_t numVal = 0;
			bAln.GetTag(tag, numVal);
			val = estd::to_string(numVal);
		}else{
			bAln.GetTag(tag, val);
		}
		Json::Value jVal;
		jVal["type"] = bib::json::toJson(type);
		jVal["val"] = bib::json::toJson(val);
		tagData[tag] = jVal;
		//std::cout << "\t" << type << std::endl;
		//std::cout << "\t" << val << std::endl;
	}
	ret["TagData"] = bib::json::toJson(tagData);
	ret["endPosition"] = bib::json::toJson(bAln.GetEndPosition());
	ret["IsDuplicate"] = bib::json::toJson(bAln.IsDuplicate());
	ret["IsFailedQC"] = bib::json::toJson(bAln.IsFailedQC());
	ret["IsFirstMate"] = bib::json::toJson(bAln.IsFirstMate());
	ret["IsMapped"] = bib::json::toJson(bAln.IsMapped());
	ret["IsMateMapped"] = bib::json::toJson(bAln.IsMateMapped());
	ret["IsMateReverseStrand"] = bib::json::toJson(bAln.IsMateReverseStrand());
	ret["IsPaired"] = bib::json::toJson(bAln.IsPaired());
	ret["IsPrimaryAlignment"] = bib::json::toJson(bAln.IsPrimaryAlignment());
	ret["IsProperPair"] = bib::json::toJson(bAln.IsProperPair());
	ret["IsReverseStrand"] = bib::json::toJson(bAln.IsReverseStrand());
	ret["IsSecondMate"] = bib::json::toJson(bAln.IsSecondMate());
	return ret;
}

}  // namespace JsonConversion
}  // namespace json
}  // namespace bib




