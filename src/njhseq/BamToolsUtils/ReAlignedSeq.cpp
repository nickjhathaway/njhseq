/*
 * ReAlignedSeq.cpp
 *
 *  Created on: Aug 16, 2019
 *      Author: nicholashathaway
 */

#include "ReAlignedSeq.hpp"

#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"
#include "njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp"
#include "njhseq/readVectorManipulation/readVectorOperations/massGetters.hpp"


namespace njhseq {



ReAlignedSeq ReAlignedSeq::genRealignment(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData,
		aligner & alignerObj,
		const std::unordered_map<std::string, uint32_t> & chromLengths,
		TwoBit::TwoBitFile & tReader,
	const genRealignmentPars & pars){
	GenomicRegion gRegion(bAln, refData);
	gRegion.reverseSrand_ = false;
	uint32_t softClipLeft = 0;
	uint32_t softClipRight = 0;
	if(pars.adjustForSoftClipping){
		if('S' == bAln.CigarData.front().Type ){
			softClipLeft += bAln.CigarData.front().Length;
		}
		if('S' == bAln.CigarData.back().Type ){
			softClipRight += bAln.CigarData.front().Length;
		}
	}
	BedUtility::extendLeftRight(gRegion, pars.extendAmount + softClipLeft, pars.extendAmount + softClipRight, chromLengths.at(gRegion.chrom_));
	auto qSeq = bamAlnToSeqInfo(bAln, true);
	auto rSeq = gRegion.extractSeq(tReader);
	uint64_t maxLen = alignerObj.parts_.maxSize_;
	readVec::getMaxLength(qSeq, maxLen);
	readVec::getMaxLength(rSeq, maxLen);
	alignerObj.parts_.setMaxSize(maxLen);
	alignerObj.alignCacheGlobal(rSeq, qSeq);
	uint32_t queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-");
	uint32_t queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-");
	uint32_t queryAlnEnd = queryAlnLastBase + 1;
	uint32_t realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
	uint32_t realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
	uint32_t realRefEnd = realRefLastBase + 1;
	seqInfo referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
	seqInfo queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
	alignerObj.alignObjectA_.seqBase_ = referenceAln;
	alignerObj.alignObjectB_.seqBase_ = queryAln;
	seqInfo refSeq = referenceAln;
	refSeq.removeGaps();
	gRegion.start_ = gRegion.start_ + realRefStart;
	gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;
	alignerObj.profileAlignment(rSeq, qSeq, false, false, false);
	ReAlignedSeq ret;
	ret.bAln_ = bAln;
	ret.gRegion_ = gRegion;
	ret.refSeq_ = refSeq;
	ret.querySeq_ = qSeq;
	ret.alnRefSeq_ = referenceAln;
	ret.alnQuerySeq_ = queryAln;
	ret.comp_ = alignerObj.comp_;
	return ret;
}


}  // namespace njhseq
