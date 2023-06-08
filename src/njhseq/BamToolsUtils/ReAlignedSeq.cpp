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
	uint32_t insertAmount = 0;
	for(const auto & cig : bAln.CigarData){
		if(cig.Type == 'I'){
			insertAmount += cig.Length;
		}
	}

	uint32_t extend = std::max(pars.extendAmount, insertAmount);
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
	BedUtility::extendLeftRight(gRegion, extend + softClipLeft, extend + softClipRight, chromLengths.at(gRegion.chrom_));
	auto rSeq = gRegion.extractSeq(tReader);
	rSeq.name_ = gRegion.createUidFromCoordsStrand();
	auto qSeq = bamAlnToSeqInfo(bAln, true);
	uint64_t maxLen = alignerObj.parts_.maxSize_ - 1;
	readVec::getMaxLength(qSeq, maxLen);
	readVec::getMaxLength(rSeq, maxLen);
	alignerObj.parts_.setMaxSize(maxLen);
	alignerObj.alignCacheGlobal(rSeq, qSeq);
	uint32_t queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
	uint32_t queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
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
	if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front() || '-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
		uint32_t extraExtendFront = 25;
		uint32_t extraExtendEnd = 25;

		if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front()){
			extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
		}
		if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
			extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.size() - alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
		}
		BedUtility::extendLeftRight(gRegion, extraExtendFront, extraExtendEnd, chromLengths.at(gRegion.chrom_));
    rSeq = gRegion.extractSeq(tReader);
		rSeq.name_ = gRegion.createUidFromCoordsStrand();

		readVec::getMaxLength(qSeq, maxLen);
		readVec::getMaxLength(rSeq, maxLen);
		alignerObj.parts_.setMaxSize(maxLen);
		alignerObj.alignCacheGlobal(rSeq, qSeq);
		queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
		queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
		queryAlnEnd = queryAlnLastBase + 1;
		realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
		realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
		realRefEnd = realRefLastBase + 1;
		referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
		queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
		alignerObj.alignObjectA_.seqBase_ = referenceAln;
		alignerObj.alignObjectB_.seqBase_ = queryAln;
		refSeq = referenceAln;
		refSeq.removeGaps();
		gRegion.start_ = gRegion.start_ + realRefStart;
		gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;

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
	return ret;
}


ReAlignedSeq ReAlignedSeq::genRealignment(const BLASTHitTab &blastHit,
                                          const std::string &originalQuery,
                                          aligner &alignerObj,
                                          const std::unordered_map<std::string, uint32_t> &chromLengths,
                                          TwoBit::TwoBitFile &tReader,
                                          const genRealignmentPars &pars) {


  GenomicRegion gRegion = blastHit.genSubjectBed6();
	auto originalRegion = gRegion;
  gRegion.meta_.meta_.clear();

//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
//	{
//		auto rSeq = gRegion.extractSeq(tReader);
//		rSeq.name_ = gRegion.createUidFromCoordsStrand();
//		rSeq.outPutSeqAnsi(std::cout);
//	}
  uint32_t extend = pars.extendAmount;
  gRegion.reverseSrand_ = false;
  uint32_t softClipLeft = 0;
  uint32_t softClipRight = 0;
  if(pars.adjustForSoftClipping){
		if(originalRegion.reverseSrand_){
			softClipLeft += originalQuery.size() - blastHit.qEnd_;
			softClipRight += blastHit.qStart_ - 1;
		}else{
			softClipLeft += blastHit.qStart_ - 1;
			softClipRight += originalQuery.size() - blastHit.qEnd_;
		}
  }
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;
  BedUtility::extendLeftRight(gRegion, extend + softClipLeft, extend + softClipRight, chromLengths.at(gRegion.chrom_));
//	std::cout << gRegion.genBedRecordCore().toDelimStr() << std::endl;

  auto rSeq = gRegion.extractSeq(tReader);
  rSeq.name_ = gRegion.createUidFromCoordsStrand();
//	rSeq.outPutSeqAnsi(std::cout);
  auto qSeq = seqInfo(blastHit.queryName_, originalQuery);
  if(blastHit.reverseStrand()){
    qSeq.reverseComplementRead(false, true);
  }
//	qSeq.outPutSeqAnsi(std::cout);
//	std::cout << std::endl;

	uint64_t maxLen = alignerObj.parts_.maxSize_ - 1;
  readVec::getMaxLength(qSeq, maxLen);
  readVec::getMaxLength(rSeq, maxLen);
  alignerObj.parts_.setMaxSize(maxLen);
  alignerObj.alignCacheGlobal(rSeq, qSeq);
  uint32_t queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
  uint32_t queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
  uint32_t queryAlnEnd = queryAlnLastBase + 1;

	uint32_t refAlnStart = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
	uint32_t refAlnLastBase = alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
  uint32_t realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
  uint32_t realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
	//have to do these checks if the reference is shorter than the query
	if(refAlnStart > queryAlnStart){
		realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, refAlnStart);
	}
	if(refAlnLastBase < queryAlnLastBase){
		realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, refAlnLastBase);
	}
  uint32_t realRefEnd = realRefLastBase + 1;


//  alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//  alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);


  seqInfo referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
  seqInfo queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
  alignerObj.alignObjectA_.seqBase_ = referenceAln;
  alignerObj.alignObjectB_.seqBase_ = queryAln;
  seqInfo refSeq = referenceAln;
  refSeq.removeGaps();
  gRegion.start_ = gRegion.start_ + realRefStart;
  gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;
  alignerObj.profileAlignment(rSeq, qSeq, false, false, false);
	if (pars.adjustForSoftClipping &&
			('-' == alignerObj.alignObjectA_.seqBase_.seq_.front() || '-' == alignerObj.alignObjectA_.seqBase_.seq_.back())) {
		uint32_t extraExtendFront = 25;
		uint32_t extraExtendEnd = 25;

    if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front()){
      extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
    }
    if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
      extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.size() - alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
    }
    BedUtility::extendLeftRight(gRegion, extraExtendFront, extraExtendEnd, chromLengths.at(gRegion.chrom_));
    rSeq = gRegion.extractSeq(tReader);
    rSeq.name_ = gRegion.createUidFromCoordsStrand();

    readVec::getMaxLength(qSeq, maxLen);
    readVec::getMaxLength(rSeq, maxLen);
    alignerObj.parts_.setMaxSize(maxLen);
    alignerObj.alignCacheGlobal(rSeq, qSeq);
    queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
    queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
    queryAlnEnd = queryAlnLastBase + 1;
    realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
    realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
    realRefEnd = realRefLastBase + 1;
    referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
    queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
    alignerObj.alignObjectA_.seqBase_ = referenceAln;
    alignerObj.alignObjectB_.seqBase_ = queryAln;
    refSeq = referenceAln;
    refSeq.removeGaps();
    gRegion.start_ = gRegion.start_ + realRefStart;
    gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;

    alignerObj.profileAlignment(rSeq, qSeq, false, false, false);
  }

  ReAlignedSeq ret;
  //ret.bAln_ = bAln; /**@todo need to generate a bAln_ object */

  ret.gRegion_ = gRegion;
  if(blastHit.reverseStrand()){
    ret.gRegion_.reverseSrand_ = true;
  }
  ret.refSeq_ = refSeq;
  ret.querySeq_ = qSeq;
  ret.alnRefSeq_ = referenceAln;
  ret.alnQuerySeq_ = queryAln;
  ret.comp_ = alignerObj.comp_;
//	std::cout << originalRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//	std::cout << ret.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
  return ret;
}

ReAlignedSeq ReAlignedSeq::genRealignment(const BLASTHitTab &blastHit,
																					const std::string &originalQuery,
																					aligner &alignerObj,
																					const seqInfo & seq,
																					const genRealignmentPars &pars) {


	GenomicRegion gRegion = blastHit.genSubjectBed6();
	auto originalRegion = gRegion;
	gRegion.meta_.meta_.clear();

	uint32_t extend = pars.extendAmount;
	gRegion.reverseSrand_ = false;
	uint32_t softClipLeft = 0;
	uint32_t softClipRight = 0;
	if (pars.adjustForSoftClipping) {
		if (originalRegion.reverseSrand_) {
			softClipLeft += originalQuery.size() - blastHit.qEnd_;
			softClipRight += blastHit.qStart_ - 1;
		} else {
			softClipLeft += blastHit.qStart_ - 1;
			softClipRight += originalQuery.size() - blastHit.qEnd_;
		}
	}
	BedUtility::extendLeftRight(gRegion, extend + softClipLeft, extend + softClipRight, seq.seq_.length());
	auto rSeq = seq.getSubRead(gRegion.start_, gRegion.getLen());
	rSeq.name_ = gRegion.createUidFromCoordsStrand();
	auto qSeq = seqInfo(blastHit.queryName_, originalQuery);
	if(blastHit.reverseStrand()){
		qSeq.reverseComplementRead(false, true);
	}
	uint64_t maxLen = alignerObj.parts_.maxSize_ - 1;
	readVec::getMaxLength(qSeq, maxLen);
	readVec::getMaxLength(rSeq, maxLen);
	alignerObj.parts_.setMaxSize(maxLen);
	alignerObj.alignCacheGlobal(rSeq, qSeq);
	uint32_t queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
	uint32_t queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
	uint32_t queryAlnEnd = queryAlnLastBase + 1;
	uint32_t realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
	uint32_t realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
	uint32_t realRefEnd = realRefLastBase + 1;

	//
//  alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//  alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
	//

	seqInfo referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
	seqInfo queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
	alignerObj.alignObjectA_.seqBase_ = referenceAln;
	alignerObj.alignObjectB_.seqBase_ = queryAln;
	seqInfo refSeq = referenceAln;
	refSeq.removeGaps();
	gRegion.start_ = gRegion.start_ + realRefStart;
	gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;
	alignerObj.profileAlignment(rSeq, qSeq, false, false, false);
	if(pars.adjustForSoftClipping && ('-' == alignerObj.alignObjectA_.seqBase_.seq_.front() || '-' == alignerObj.alignObjectA_.seqBase_.seq_.back())){
		uint32_t extraExtendFront = 25;
		uint32_t extraExtendEnd = 25;

		if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front()){
			extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
		}
		if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
			extraExtendFront = alignerObj.alignObjectA_.seqBase_.seq_.size() - alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
		}
		BedUtility::extendLeftRight(gRegion, extraExtendFront, extraExtendEnd, len(seq));
		rSeq = seq.getSubRead(gRegion.start_, gRegion.getLen());
		rSeq.name_ = gRegion.createUidFromCoordsStrand();

		readVec::getMaxLength(qSeq, maxLen);
		readVec::getMaxLength(rSeq, maxLen);
		alignerObj.parts_.setMaxSize(maxLen);
		alignerObj.alignCacheGlobal(rSeq, qSeq);
		queryAlnStart = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
		queryAlnLastBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
		queryAlnEnd = queryAlnLastBase + 1;
		realRefStart  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnStart);
		realRefLastBase  = getRealPosForAlnPos(alignerObj.alignObjectA_.seqBase_.seq_, queryAlnLastBase);
		realRefEnd = realRefLastBase + 1;
		referenceAln = alignerObj.alignObjectA_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
		queryAln = alignerObj.alignObjectB_.seqBase_.getSubRead(queryAlnStart, queryAlnEnd - queryAlnStart);
		alignerObj.alignObjectA_.seqBase_ = referenceAln;
		alignerObj.alignObjectB_.seqBase_ = queryAln;
		refSeq = referenceAln;
		refSeq.removeGaps();
		gRegion.start_ = gRegion.start_ + realRefStart;
		gRegion.end_ = gRegion.start_ + realRefEnd - realRefStart;

		alignerObj.profileAlignment(rSeq, qSeq, false, false, false);
	}

	ReAlignedSeq ret;
	//ret.bAln_ = bAln; /**@todo need to generate a bAln_ object */
	ret.gRegion_ = gRegion;
	if(blastHit.reverseStrand()){
		ret.gRegion_.reverseSrand_ = true;
	}
	ret.refSeq_ = refSeq;
	ret.querySeq_ = qSeq;
	ret.alnRefSeq_ = referenceAln;
	ret.alnQuerySeq_ = queryAln;
	ret.comp_ = alignerObj.comp_;
	return ret;
}

std::vector<std::shared_ptr<ReAlignedSeq>> ReAlignedSeq::getUniqueLocationResults(
    std::vector<std::shared_ptr<ReAlignedSeq>> & alnResults) {
  if (alnResults.size() <= 1) {
    return alnResults;
  }
  njh::sort(alnResults,
            [](const std::shared_ptr<ReAlignedSeq> &results1,
               const std::shared_ptr<ReAlignedSeq> &results2) {
              if (results1->gRegion_.createUidFromCoords() == results2->gRegion_.createUidFromCoords()) {
                return results1->comp_.distances_.eventBasedIdentity_ > results2->comp_.distances_.eventBasedIdentity_;
              } else {
                return results1->gRegion_.createUidFromCoords() < results2->gRegion_.createUidFromCoords();
              }
            });
//  std::cout << "alnResults.size(): " << alnResults.size() << std::endl;
  std::vector<std::shared_ptr<ReAlignedSeq>>  ret;
  ret.emplace_back(alnResults.front());
  for(const auto pos : iter::range<uint32_t>(1, alnResults.size())){
//    std::cout << "ret.back()->gRegion_.createUidFromCoords(): " << ret.back()->gRegion_.createUidFromCoords() << std::endl;
//    std::cout << "alnResults[pos]->gRegion_.createUidFromCoords(): " << alnResults[pos]->gRegion_.createUidFromCoords() << std::endl;

    if(ret.back()->gRegion_.createUidFromCoords() != alnResults[pos]->gRegion_.createUidFromCoords()){
      ret.emplace_back(alnResults[pos]);
    }
  }
//  std::cout << "ret.size(): " << ret.size() << std::endl;
//  std::cout << std::endl;
  return ret;
}


}  // namespace njhseq
