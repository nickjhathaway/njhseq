/*
 * longestTranslation.cpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "longestTranslation.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp"
#include "njhseq/objects/Meta/MetaDataInName.hpp"

namespace njhseq {

LongestTranslation LongestTranslation::translate(const seqInfo & seq, const Pars & pars) {
	LongestTranslation ret;
	ret.seqLen_ = len(seq);
	std::shared_ptr<readVecTrimmer::BreakUpRes> longest;
	njh::PatPosFinder pFinder(R"(\*+)");
	for (const auto comp : {false, true}) {
		for (uint32_t start = 0; start < 3; ++start) {
			auto prot = seq.translateRet(comp, comp, start, pars.forceMStart_);
			auto possibleSubSeqs = readVecTrimmer::breakUpSeqOnPat(prot, pFinder);
			for (const auto & possible : possibleSubSeqs) {
				if (nullptr == longest || len(possible.seqBase_) > len(longest->seqBase_)) {
					longest = std::make_shared<readVecTrimmer::BreakUpRes>(possible);
					ret.frame_ = start;
					ret.frameProteinLen_ = len(prot);
					ret.isRevComp_ = comp;
				}
			}
		}
	}
	if (nullptr == longest) {
		ret.protein_ = seqInfo(seq.name_, "");
		return ret;
	}
	ret.protein_ = longest->seqBase_;
	ret.aaStart_ = longest->start_;
	ret.aaStop_ = longest->end_;
	ret.seqStart_ = longest->start_ * 3 + ret.frame_;
	ret.seqStop_ = longest->end_ * 3 + ret.frame_;
	//not ending at the end of the frame's translation means it ended on a stop
	if (!pars.removeTrailingStop_ && longest->end_ != ret.frameProteinLen_) {
		ret.protein_.append("*");
		++ret.aaStop_;
		ret.seqStop_ += 3;
	}
	return ret;
}

void LongestTranslation::markName(bool overWriteMeta) {
	MetaDataInName meta;
	bool nameHasMeta = false;
	if (MetaDataInName::nameHasMetaData(protein_.name_)) {
		meta = MetaDataInName(protein_.name_);
		nameHasMeta = true;
	}
	meta.addMeta("frame", frame_, overWriteMeta);
	meta.addMeta("seqStart", seqStart_, overWriteMeta);
	meta.addMeta("seqStop", seqStop_, overWriteMeta);
	meta.addMeta("seqLen", seqLen_, overWriteMeta);
	meta.addMeta("aaStart", aaStart_, overWriteMeta);
	meta.addMeta("aaStop", aaStop_, overWriteMeta);
	meta.addMeta("aaLen", frameProteinLen_, overWriteMeta);
	meta.addMeta("isRevComp", isRevComp_);
	if (nameHasMeta) {
		meta.resetMetaInName(protein_.name_);
	} else {
		protein_.name_.append(" " + meta.createMetaName());
	}
}

}  // namespace njhseq
