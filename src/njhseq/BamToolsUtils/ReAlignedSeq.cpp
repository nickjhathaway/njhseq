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
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"
#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"

namespace njhseq {



ReAlignedSeq ReAlignedSeq::genRealignment(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData,
		aligner & alignerObj,
		const std::unordered_map<std::string, uint32_t> & chromLengths,
		TwoBit::TwoBitFile & tReader,
	const genRealignmentPars & pars){
	GenomicRegion gRegion(bAln, refData);
//
//	auto endPositon = bAln.GetEndPosition();
//	std::cout << "endPositon:" << endPositon << std::endl;
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

//	std::cout << bAln.AlignedBases << std::endl;
	uint32_t insertAmount = 0;
	for(const auto & cig : bAln.CigarData){
//		std::cout << cig.Type << cig.Length;
		if(cig.Type == 'I'){
			insertAmount += cig.Length;
		}
	}
//	std::cout << std::endl;
	uint32_t extend = std::max(pars.extendAmount, insertAmount);
//	std::cout << "extend: " << extend << std::endl;
//	std::cout << "pars.extendAmount: " << pars.extendAmount << std::endl;
	auto alnLocal = bamAlnToAlnInfoLocal(bAln);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	gRegion.reverseSrand_ = false;
	uint32_t softClipLeft = 0;
	uint32_t softClipRight = 0;
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(pars.adjustForSoftClipping){
		if('S' == bAln.CigarData.front().Type ){
			softClipLeft += bAln.CigarData.front().Length;
		}
		if('S' == bAln.CigarData.back().Type ){
			softClipRight += bAln.CigarData.front().Length;
		}
	}
//	std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
	BedUtility::extendLeftRight(gRegion, extend + softClipLeft, extend + softClipRight, chromLengths.at(gRegion.chrom_));
//	std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
	auto rSeq = gRegion.extractSeq(tReader);
	rSeq.name_ = gRegion.createUidFromCoordsStrand();
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	auto qSeq = bamAlnToSeqInfo(bAln, true);
//	rSeq.outPutSeqAnsi(std::cout);
//	qSeq.outPutSeqAnsi(std::cout);

//	auto rSeqCopy = rSeq;
//	auto qSeqCopy = qSeq;
//	alignCalc::rearrangeLocal(rSeqCopy.seq_,
//			qSeqCopy.seq_, '-', alnLocal.begin()->second);
//	alignCalc::rearrangeLocal(rSeqCopy.qual_,
//			qSeqCopy.qual_, 0, alnLocal.begin()->second);
//	rSeqCopy.name_ = gRegion.createUidFromCoordsStrand();
//	rSeqCopy.outPutSeqAnsi(std::cout);
//	qSeqCopy.outPutSeqAnsi(std::cout);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	uint64_t maxLen = alignerObj.parts_.maxSize_;
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << "maxLen: " << maxLen << std::endl;
	readVec::getMaxLength(qSeq, maxLen);
	readVec::getMaxLength(rSeq, maxLen);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << "maxLen: " << maxLen << std::endl;
	alignerObj.parts_.setMaxSize(maxLen);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	rSeq.outPutSeqAnsi(std::cout);
//	qSeq.outPutSeqAnsi(std::cout);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	alignerObj.alignCacheGlobal(rSeq, qSeq);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	uint32_t queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-");
	uint32_t queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-");
	uint32_t queryAlnEnd = queryAlnLastBase + 1;
	uint32_t realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
	uint32_t realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
	uint32_t realRefEnd = realRefLastBase + 1;
////std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//	alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
	seqInfo referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
	seqInfo queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
//	alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//	alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	alignerObj.alignObjectA_.seqBase_ = referenceAln;
	alignerObj.alignObjectB_.seqBase_ = queryAln;
	seqInfo refSeq = referenceAln;
	refSeq.removeGaps();
//	std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
	gRegion.start_ = gRegion.start_ + realRefStart;
	gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;
//	std::cout << gRegion.genBedRecordCore().toDelimStr()  << std::endl;
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	alignerObj.profileAlignment(rSeq, qSeq, false, false, false);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front() || '-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
		uint32_t extraExtendFront = 25;
		uint32_t extraExtendEnd = 25;

		if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front()){
			extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
		}
		if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
			extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.size() - alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
		}
//		std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
		BedUtility::extendLeftRight(gRegion, extraExtendFront, extraExtendEnd, chromLengths.at(gRegion.chrom_));
//		std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
		auto rSeq = gRegion.extractSeq(tReader);
		rSeq.name_ = gRegion.createUidFromCoordsStrand();

		auto qSeq = bamAlnToSeqInfo(bAln, true);
//		rSeq.outPutSeqAnsi(std::cout);
//		qSeq.outPutSeqAnsi(std::cout);

//		auto rSeqCopy = rSeq;
//		auto qSeqCopy = qSeq;
//		alignCalc::rearrangeLocal(rSeqCopy.seq_,
//				qSeqCopy.seq_, '-', alnLocal.begin()->second);
//		alignCalc::rearrangeLocal(rSeqCopy.qual_,
//				qSeqCopy.qual_, 0, alnLocal.begin()->second);
//		rSeqCopy.name_ = gRegion.createUidFromCoordsStrand();
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		rSeqCopy.outPutSeqAnsi(std::cout);
//		qSeqCopy.outPutSeqAnsi(std::cout);

		readVec::getMaxLength(qSeq, maxLen);
		readVec::getMaxLength(rSeq, maxLen);
		alignerObj.parts_.setMaxSize(maxLen);
		alignerObj.alignCacheGlobal(rSeq, qSeq);
		queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-");
		queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-");
		queryAlnEnd = queryAlnLastBase + 1;
		realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
		realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
		realRefEnd = realRefLastBase + 1;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//		alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
		queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
//		alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//		alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

		alignerObj.alignObjectA_.seqBase_ = referenceAln;
		alignerObj.alignObjectB_.seqBase_ = queryAln;
		refSeq = referenceAln;
		refSeq.removeGaps();
//		std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
		gRegion.start_ = gRegion.start_ + realRefStart;
		gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;
//		std::cout << gRegion.genBedRecordCore().toDelimStr()  << std::endl;

		alignerObj.profileAlignment(rSeq, qSeq, false, false, false);
	}


	ReAlignedSeq ret;
	ret.bAln_ = bAln;
	ret.gRegion_ = gRegion;
	ret.refSeq_ = refSeq;
	ret.querySeq_ = qSeq;
	ret.alnRefSeq_ = referenceAln;
	ret.alnQuerySeq_ = queryAln;
	ret.comp_ = alignerObj.comp_;
//	if("Pf3D7_04_v3-127475-129575.487[HapPopUIDCount=1]" == bAln.Name){
//		exit(1);
//	}
	return ret;
}


}  // namespace njhseq
