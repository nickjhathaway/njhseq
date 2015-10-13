#pragma once
/*
 * BamToolsUtils.hpp
 *
 *  Created on: May 28, 2015
 *      Author: nick
 */

#include <api/BamAlignment.h>

#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"

namespace bibseq {

seqInfo bamAlnToSeqInfo(const BamTools::BamAlignment & aln);

void bamWriteFasta(const BamTools::BamAlignment & aln,
		std::ostream & out);
void bamWriteFastq(const BamTools::BamAlignment & aln,
		std::ostream & out);


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

