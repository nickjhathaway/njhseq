#pragma once
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
 * BamToolsUtils.hpp
 *
 *  Created on: May 28, 2015
 *      Author: nick
 */

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>


#include "njhseq/common.h"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "njhseq/alignment/alnCache.h"
#include "njhseq/objects/dataContainers/tables/table.hpp"

namespace njhseq {


uint32_t getSoftClipAmount(const BamTools::BamAlignment & bAln);



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
void checkBamOpenThrow(BamTools::BamReader & bReader, const bfs::path & bamFnp, const std::string & funcName);

void loadBamIndexThrow(BamTools::BamReader & bReader);
void loadBamIndexThrow(BamTools::BamReader & bReader, const std::string & funcName);

table refDataVecToTab(const std::vector<BamTools::RefData> & refInfos);
std::vector<BamTools::RefData> tabToRefDataVec(const table & refDataTab);

void logAlnInfo(std::ostream & out, BamTools::RefVector & refInfo,
		const BamTools::BamAlignment & aln, std::string indent = "");



void setBamFileRegionThrow(BamTools::BamReader & bReader, const GenomicRegion & region);

void checkBamFilesForIndexesAndAbilityToOpen(const std::vector<bfs::path> & bamFnps);


} /* namespace njhseq */



namespace njh {
namespace json {
namespace JsonConversion {
template<>
inline Json::Value toJson(const BamTools::CigarOp & cOps){
	Json::Value ret;
	ret["class"] = njh::getTypeName(cOps);
	ret["Length"] = njh::json::toJson(cOps.Length);
	ret["Type"] = njh::json::toJson(cOps.Type);
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
	ret["class"] = njh::getTypeName(bAln);
	ret["AlignedBases"] = njh::json::toJson(bAln.AlignedBases);
	ret["AlignmentFlag"] = njh::json::toJson(bAln.AlignmentFlag);
	ret["Bin"] = njh::json::toJson(bAln.Bin);
	ret["CigarData"] = njh::json::toJson(bAln.CigarData);
	ret["Filename"] = njh::json::toJson(bAln.Filename);
	ret["InsertSize"] = njh::json::toJson(bAln.InsertSize);
	ret["Length"] = njh::json::toJson(bAln.Length);
	ret["MapQuality"] = njh::json::toJson(bAln.MapQuality);
	ret["MatePosition"] = njh::json::toJson(bAln.MatePosition);
	ret["MateRefID"] = njh::json::toJson(bAln.MateRefID);
	ret["Name"] = njh::json::toJson(bAln.Name);
	ret["Position"] = njh::json::toJson(bAln.Position);
	ret["Qualities"] = njh::json::toJson(bAln.Qualities);
	ret["QueryBases"] = njh::json::toJson(bAln.QueryBases);
	ret["RefID"] = njh::json::toJson(bAln.RefID);
	ret["TagDataRaw"] = njh::json::toJson(bAln.TagData);
	ret["TagNames"] = njh::json::toJson(bAln.GetTagNames());
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
		jVal["type"] = njh::json::toJson(type);
		jVal["val"] = njh::json::toJson(val);
		tagData[tag] = jVal;
		//std::cout << "\t" << type << std::endl;
		//std::cout << "\t" << val << std::endl;
	}
	ret["TagData"] = njh::json::toJson(tagData);
	ret["endPosition"] = njh::json::toJson(bAln.GetEndPosition());
	ret["IsDuplicate"] = njh::json::toJson(bAln.IsDuplicate());
	ret["IsFailedQC"] = njh::json::toJson(bAln.IsFailedQC());
	ret["IsFirstMate"] = njh::json::toJson(bAln.IsFirstMate());
	ret["IsMapped"] = njh::json::toJson(bAln.IsMapped());
	ret["IsMateMapped"] = njh::json::toJson(bAln.IsMateMapped());
	ret["IsMateReverseStrand"] = njh::json::toJson(bAln.IsMateReverseStrand());
	ret["IsPaired"] = njh::json::toJson(bAln.IsPaired());
	ret["IsPrimaryAlignment"] = njh::json::toJson(bAln.IsPrimaryAlignment());
	ret["IsProperPair"] = njh::json::toJson(bAln.IsProperPair());
	ret["IsReverseStrand"] = njh::json::toJson(bAln.IsReverseStrand());
	ret["IsSecondMate"] = njh::json::toJson(bAln.IsSecondMate());
	return ret;
}

}  // namespace JsonConversion
}  // namespace json
}  // namespace njh




