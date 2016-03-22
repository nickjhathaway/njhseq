#pragma once
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
 * BamToolsUtils.hpp
 *
 *  Created on: May 28, 2015
 *      Author: nick
 */

#include <api/BamAlignment.h>

#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
#include "bibseq/alignment/alnCache/alnInfoHolder.hpp"

namespace bibseq {

seqInfo bamAlnToSeqInfo(const BamTools::BamAlignment & aln, bool complement);

void bamWriteFasta(const BamTools::BamAlignment & aln, std::ostream & out,
		bool complement);
void bamWriteFastq(const BamTools::BamAlignment & aln, std::ostream & out,
		bool complement);

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

bool IsBamSorted(const std::string & filename);

} /* namespace bibseq */

