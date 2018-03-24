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
#include "bibseq/alignment/aligner/alignCalc.hpp"
#include "aligner.hpp"
#include "bibseq/helpers/seqUtil.hpp"

namespace bibseq {


void aligner::resetAlnCache(){
	alnHolder_ = alnInfoMasterHolder(parts_.gapScores_, parts_.scoring_);
}

aligner::aligner() :
		parts_(alnParts()), alnHolder_(
				alnInfoMasterHolder(parts_.gapScores_, parts_.scoring_)) {
	countEndGaps_ = false;
	weighHomopolymers_ = false;
	setDefaultQualities();
}

aligner::aligner(uint64_t maxSize, const gapScoringParameters& gapPars,
		const substituteMatrix& scoreMatrix) :
		parts_(maxSize, gapPars, scoreMatrix), alnHolder_(
				alnInfoMasterHolder(gapPars, scoreMatrix)) {
	countEndGaps_ = false;
	weighHomopolymers_ = false;
	setDefaultQualities();
}

aligner::aligner(uint64_t maxSize, const gapScoringParameters& gapPars,
		const substituteMatrix& scoreMatrix, bool countEndGaps) :
		parts_(maxSize, gapPars, scoreMatrix), alnHolder_(
				alnInfoMasterHolder(gapPars, scoreMatrix)), countEndGaps_(countEndGaps) {
	weighHomopolymers_ = false;
	setDefaultQualities();
}

aligner::aligner(uint64_t maxSize, const gapScoringParameters & gapPars,
		const substituteMatrix& subMatrix, const KmerMaps& kmaps,
		QualScorePars qScorePars, bool countEndGaps, bool weighHomopolymers) :
		parts_(maxSize, gapPars, subMatrix), alnHolder_(gapPars, subMatrix), kMaps_(
				kmaps), qScorePars_(qScorePars), countEndGaps_(countEndGaps), weighHomopolymers_(
				weighHomopolymers) {
}

void aligner::setGapScoring(const gapScoringParameters & gapPars){
	parts_.gapScores_ = gapPars;
	if(alnHolder_.globalHolder_.find(gapPars.uniqueIdentifer_) == alnHolder_.globalHolder_.end()){
		alnHolder_.addHolder(gapPars, parts_.scoring_);
	}
}

//void aligner::setMatchScoring(const substituteMatrix& subMatrix){
//	parts_.scoring_ = subMatrix;
//	alnHolder_.addHolder(parts_.gapScores_, subMatrix);
//}



void aligner::alignScoreLocal(const std::string& firstSeq,
		const std::string& secondSeq) {
	alignCalc::runSmithSave(firstSeq, secondSeq, parts_);
	++numberOfAlingmentsDone_;
}

void aligner::alignScoreCacheLocal(const std::string& firstSeq,
		const std::string& secondSeq) {
	if (alnHolder_.localHolder_[parts_.gapScores_.uniqueIdentifer_].getAlnInfo(
			firstSeq, secondSeq, parts_.lHolder_)) {
		parts_.score_ = parts_.lHolder_.score_;
		comp_.alnScore_ = parts_.score_;
	} else {
		alignCalc::runSmithSave(firstSeq, secondSeq, parts_);
		alnHolder_.localHolder_[parts_.gapScores_.uniqueIdentifer_].addAlnInfo(
				firstSeq, secondSeq, parts_.lHolder_);
		++numberOfAlingmentsDone_;
	}
}
void aligner::alignScoreGlobal(const std::string& firstSeq,
		const std::string& secondSeq) {
	alignCalc::runNeedleSave(firstSeq, secondSeq, parts_);
	++numberOfAlingmentsDone_;
}

void aligner::alignScoreGlobalNoInternalGaps(const std::string& firstSeq,
		const std::string& secondSeq) {
	alignCalc::runNeedleOnlyEndGapsSave(firstSeq, secondSeq, parts_);
	++numberOfAlingmentsDone_;
}


void aligner::alignScoreCacheGlobal(const std::string& firstSeq,
		const std::string& secondSeq) {
	if (alnHolder_.globalHolder_[parts_.gapScores_.uniqueIdentifer_].getAlnInfo(
			firstSeq, secondSeq, parts_.gHolder_)) {
		parts_.score_ = parts_.gHolder_.score_;
		comp_.alnScore_ = parts_.score_;
	} else {
		alignCalc::runNeedleSave(firstSeq, secondSeq, parts_);
		alnHolder_.globalHolder_[parts_.gapScores_.uniqueIdentifer_].addAlnInfo(
				firstSeq, secondSeq, parts_.gHolder_);
		++numberOfAlingmentsDone_;
	}
}

void aligner::alignScore(const std::string& firstSeq, const std::string& secondSeq,
		bool local) {
	if (local) {
		alignCalc::runSmithSave(firstSeq, secondSeq, parts_);
	} else {
		alignCalc::runNeedleSave(firstSeq, secondSeq, parts_);
	}
	++numberOfAlingmentsDone_;
}



void aligner::alignScoreCache(const std::string& firstSeq,
		const std::string& secondSeq, bool local) {
	if (local) {
		if (alnHolder_.localHolder_[parts_.gapScores_.uniqueIdentifer_].getAlnInfo(
				firstSeq, secondSeq, parts_.lHolder_)) {
			parts_.score_ = parts_.lHolder_.score_;
			comp_.alnScore_ = parts_.score_;
		} else {
			alignCalc::runSmithSave(firstSeq, secondSeq, parts_);
			alnHolder_.localHolder_[parts_.gapScores_.uniqueIdentifer_].addAlnInfo(
					firstSeq, secondSeq, parts_.lHolder_);
			++numberOfAlingmentsDone_;
		}
	} else {
		if (alnHolder_.globalHolder_[parts_.gapScores_.uniqueIdentifer_].getAlnInfo(
				firstSeq, secondSeq, parts_.gHolder_)) {
			parts_.score_ = parts_.gHolder_.score_;
			comp_.alnScore_ = parts_.score_;
		} else {
			alignCalc::runNeedleSave(firstSeq, secondSeq, parts_);
			alnHolder_.globalHolder_[parts_.gapScores_.uniqueIdentifer_].addAlnInfo(
					firstSeq, secondSeq, parts_.gHolder_);
			++numberOfAlingmentsDone_;
		}
	}
}

void aligner::alignCacheLocal(const seqInfo & ref, const seqInfo & read){
	alignScoreCacheLocal(ref.seq_, read.seq_);
	rearrangeObjsLocal(ref, read);
}

void aligner::alignCacheGlobal(const seqInfo & ref, const seqInfo & read){
	alignScoreCacheGlobal(ref.seq_, read.seq_);
	rearrangeObjsGlobal(ref, read);
}

void aligner::alignCache(const seqInfo & ref, const seqInfo & read, bool local){
	alignScoreCache(ref.seq_, read.seq_, local);
	rearrangeObjs(ref, read, local);
}

void aligner::alignReg(const seqInfo & ref, const seqInfo & read, bool local){
	alignScore(ref.seq_, read.seq_, local);
	rearrangeObjs(ref, read, local);
}



void aligner::alignRegGlobalNoInternalGaps(const seqInfo & ref, const seqInfo & read){
	alignScoreGlobalNoInternalGaps(ref.seq_, read.seq_);
	rearrangeObjsGlobal(ref, read);
}

void aligner::alignRegGlobal(const seqInfo & ref, const seqInfo & read){
	alignScoreGlobal(ref.seq_, read.seq_);
	rearrangeObjsGlobal(ref, read);
}

void aligner::alignRegLocal(const seqInfo & ref, const seqInfo & read){
	alignScoreLocal(ref.seq_, read.seq_);
	rearrangeObjsLocal(ref, read);
}





std::pair<uint32_t, uint32_t> aligner::findReversePrimer(const std::string& read,
                                        				const std::string& primer){
	alignScoreCache(read, primer, true);
	return {parts_.lHolder_.localAStart_,
		parts_.lHolder_.localAStart_ + parts_.lHolder_.localASize_ - 1};
}
std::pair<uint32_t, uint32_t> aligner::findReversePrimer(const baseReadObject& read,
                                        				const baseReadObject& primer){
	return findReversePrimer(read.seqBase_.seq_, primer.seqBase_.seq_);
}

void aligner::rearrangeSeq(const std::string& firstRead,
		const std::string& secondRead, bool local) {
	alignObjectA_.seqBase_ = seqInfo("A", firstRead);
	alignObjectB_.seqBase_ = seqInfo("B", secondRead);
	if (local) {
		alignCalc::rearrangeLocal(alignObjectA_.seqBase_.seq_,
				alignObjectB_.seqBase_.seq_, '-', parts_.lHolder_);
		alignCalc::rearrangeLocal(alignObjectA_.seqBase_.qual_,
				alignObjectB_.seqBase_.qual_, 0, parts_.lHolder_);
	} else {
		alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.seq_,
				alignObjectB_.seqBase_.seq_, '-', parts_.gHolder_);
		alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.qual_,
				alignObjectB_.seqBase_.qual_, 0, parts_.gHolder_);
	}
}

void aligner::rearrangeObjs(const seqInfo& firstRead, const seqInfo& secondRead,
		bool local) {
	alignObjectA_.seqBase_ = firstRead;
	alignObjectB_.seqBase_ = secondRead;
	if (local) {
		alignCalc::rearrangeLocal(alignObjectA_.seqBase_.seq_,
				alignObjectB_.seqBase_.seq_, '-', parts_.lHolder_);
		alignCalc::rearrangeLocal(alignObjectA_.seqBase_.qual_,
				alignObjectB_.seqBase_.qual_, 0, parts_.lHolder_);
	} else {
		alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.seq_,
				alignObjectB_.seqBase_.seq_, '-', parts_.gHolder_);
		alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.qual_,
				alignObjectB_.seqBase_.qual_, 0, parts_.gHolder_);
	}
}

void aligner::rearrangeObjsLocal(const seqInfo& firstRead, const seqInfo& secondRead){
	alignObjectA_.seqBase_ = firstRead;
	alignObjectB_.seqBase_ = secondRead;
	alignCalc::rearrangeLocal(alignObjectA_.seqBase_.seq_,
			alignObjectB_.seqBase_.seq_, '-', parts_.lHolder_);
	alignCalc::rearrangeLocal(alignObjectA_.seqBase_.qual_,
			alignObjectB_.seqBase_.qual_, 0, parts_.lHolder_);
}

void aligner::rearrangeObjsGlobal(const seqInfo& firstRead, const seqInfo& secondRead){
	alignObjectA_.seqBase_ = firstRead;
	alignObjectB_.seqBase_ = secondRead;
	alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.seq_,
			alignObjectB_.seqBase_.seq_, '-', parts_.gHolder_);
	alignCalc::rearrangeGlobal(alignObjectA_.seqBase_.qual_,
			alignObjectB_.seqBase_.qual_, 0, parts_.gHolder_);
}


void aligner::resetCounts() {
  comp_.resetCounts();
}

void aligner::resetAlignmentInfo() {
  resetCounts();
}

void aligner::setQual(QualScorePars pars) {
  qScorePars_ = pars;
}


const comparison & aligner::profilePrimerAlignment(const seqInfo& objectA,
                            const seqInfo& objectB){
  resetAlignmentInfo();
	uint32_t firstOffset = 0;
	uint32_t secondOffset = 0;
	uint32_t gappedBasesInA = 0;
	uint32_t gappedBasesInB = 0;

  for (uint32_t i = 0; i < len(alignObjectA_); ++i) {
  	//gap in alignObjectA, normally reference sequence
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectB_.seqBase_.seq_.substr(i, 1),
					alignObjectB_.seqBase_.qual_[i], true);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
        ++i;
        ++firstOffset;
      }
			bool endGap = (newGap.startPos_ + newGap.size_ >= len(alignObjectA_)) || 0 == newGap.startPos_;
			if (!endGap) {
				gappedBasesInA += newGap.size_;
				handleGapCountingInA(newGap);
				comp_.distances_.alignmentGaps_.insert(
						std::make_pair(newGap.startPos_, newGap));
			} else if (countEndGaps_) {
				handleGapCountingInA(newGap);
				comp_.distances_.alignmentGaps_.insert(
						std::make_pair(newGap.startPos_, newGap));
			}
      continue;
    }
    //gap in alignObjectB, normally query sequence
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectA_.seqBase_.seq_.substr(i, 1),
					alignObjectA_.seqBase_.qual_[i], false);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
        ++i;
        ++secondOffset;
      }
      bool endGap = (newGap.startPos_ + newGap.size_ >= len(alignObjectB_)) || 0 == newGap.startPos_;
			if (!endGap) {
				gappedBasesInB += newGap.size_;
				handleGapCountingInB(newGap);
				comp_.distances_.alignmentGaps_.insert(
										std::make_pair(newGap.startPos_, newGap));
			}else if (countEndGaps_) {
				handleGapCountingInB(newGap);
				comp_.distances_.alignmentGaps_.insert(
										std::make_pair(newGap.startPos_, newGap));
			}
      continue;
    }
    if (0 > parts_.scoring_.mat_[alignObjectA_.seqBase_.seq_[i]]
                         [alignObjectB_.seqBase_.seq_[i]]) {
      ++comp_.hqMismatches_;
      comp_.distances_.mismatches_.insert(std::make_pair(
                  i,
                  mismatch(
                      alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                      objectA.getLeadQual(i - firstOffset, qScorePars_.qualThresWindow_),
                      objectA.getTrailQual(i - firstOffset,qScorePars_.qualThresWindow_), i - firstOffset,
                      alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                      objectB.getLeadQual(i - secondOffset,qScorePars_.qualThresWindow_),
                      objectB.getTrailQual(i - secondOffset,qScorePars_.qualThresWindow_), i - secondOffset,
                      0,0)));
    } else {
      ++comp_.highQualityMatches_;
    }
  }

	uint32_t gappedLeftEndA = 0;
	uint32_t gappedRightEndA = 0;
	if (alignObjectA_.seqBase_.seq_.front() == '-') {
		gappedLeftEndA += countBeginChar(alignObjectA_.seqBase_.seq_);
	}
	if (alignObjectA_.seqBase_.seq_.back() == '-') {
		gappedRightEndA += countEndChar(alignObjectA_.seqBase_.seq_);
	}
	uint32_t gappedEndA = gappedRightEndA + gappedLeftEndA;

	uint32_t gappedLeftEndB = 0;
	uint32_t gappedRightEndB = 0;
	if (alignObjectB_.seqBase_.seq_.front() == '-') {
		gappedLeftEndB += countBeginChar(alignObjectB_.seqBase_.seq_);
	}
	if (alignObjectB_.seqBase_.seq_.back() == '-') {
		gappedRightEndB += countEndChar(alignObjectB_.seqBase_.seq_);
	}
	uint32_t gappedEndB = gappedLeftEndB + gappedRightEndB;

	uint32_t gappedEnd = gappedEndA + gappedEndB;

	//objectA (ref)
	comp_.distances_.ref_.covered_ = len(alignObjectA_) - gappedBasesInA - gappedEnd;
	comp_.distances_.ref_.coverage_ = static_cast<double>(comp_.distances_.ref_.covered_)/len(objectA);
	comp_.distances_.ref_.identities_ = comp_.highQualityMatches_ + comp_.lowQualityMatches_;
	comp_.distances_.ref_.identity_ = static_cast<double>(comp_.distances_.ref_.identities_)/len(objectA);
	//objectB (query)
	comp_.distances_.query_.covered_ = len(alignObjectB_) - gappedBasesInB - gappedEnd;
	comp_.distances_.query_.coverage_ = static_cast<double>(comp_.distances_.query_.covered_)/len(objectB);
	comp_.distances_.query_.identities_ = comp_.highQualityMatches_ + comp_.lowQualityMatches_;
	comp_.distances_.query_.identity_ = static_cast<double>(comp_.distances_.query_.identities_)/len(objectB);
	//of the overlapping alignment
	comp_.distances_.basesInAln_ = len(alignObjectA_) - gappedEnd;
  comp_.distances_.percentMatch_ = (comp_.highQualityMatches_ + comp_.lowQualityMatches_)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.percentMismatch_ = (comp_.hqMismatches_ + comp_.lqMismatches_)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.percentGaps_ = (gappedBasesInA + gappedBasesInB)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.overLappingEvents_ = comp_.highQualityMatches_
			+ comp_.lowQualityMatches_ + comp_.hqMismatches_ + comp_.lqMismatches_
			+ comp_.lowKmerMismatches_ + comp_.distances_.alignmentGaps_.size();
  if(comp_.distances_.overLappingEvents_ == 0){
  		comp_.distances_.eventBasedIdentity_ = 0;
  }else{
  		comp_.distances_.eventBasedIdentity_ = (comp_.highQualityMatches_ + comp_.lowQualityMatches_)/static_cast<double>(comp_.distances_.overLappingEvents_);
  }
  comp_.setEventBaseIdentityHq();
  comp_.refName_ = objectA.name_;
  comp_.queryName_ = objectB.name_;
  comp_.alnScore_ = parts_.score_;
  return comp_;
}


size_t aligner::getAlignPosForSeqAPos(size_t seqAPos){
	return getAlnPosForRealPos(alignObjectA_.seqBase_.seq_, seqAPos);
}
size_t aligner::getAlignPosForSeqBPos(size_t seqBPos){
	return getAlnPosForRealPos(alignObjectB_.seqBase_.seq_, seqBPos);

}
size_t aligner::getSeqPosForAlnAPos(size_t alnAPos){
	return getRealPosForAlnPos(alignObjectA_.seqBase_.seq_, alnAPos);
}
size_t aligner::getSeqPosForAlnBPos(size_t alnBPos){
	return getRealPosForAlnPos(alignObjectB_.seqBase_.seq_, alnBPos);
}

const comparison & aligner::profileAlignment(const seqInfo& objectA,
		const seqInfo& objectB, bool checkKmer, bool usingQuality,
		bool doingMatchQuality, uint32_t start, uint32_t stop) {
  resetAlignmentInfo();
	uint32_t firstOffset = 0;
	uint32_t secondOffset = 0;
	uint32_t gappedBasesInA = 0;
	uint32_t gappedBasesInB = 0;
	uint32_t numberOfGaps = 0;
  //if the start is the very beginning, need to calculate the offset
  if (start != 0) {
    for (uint32_t i = 0; i < start; ++i) {
    	//gap in alignObjectA, normally reference sequence
      if (alignObjectA_.seqBase_.seq_[i] == '-') {
        ++firstOffset;
        while (alignObjectA_.seqBase_.seq_[i + 1] == '-' && (i + 1) != start) {
          ++i;
          ++firstOffset;
        }
        continue;
      }
      //gap in alignObjectB, normally query sequence
      if (alignObjectB_.seqBase_.seq_[i] == '-') {
        ++secondOffset;
        while (alignObjectB_.seqBase_.seq_[i + 1] == '-' && (i + 1) != start) {
          ++i;
          ++secondOffset;
        }
        continue;
      }
    }
  }
  //if stop not manually set this to the size of the alignment
  /**@todo consider throwing an exception or at least printing
   *  a warning if the stop is request to be greater than the sequence*/
  if (stop == 0 || stop > alignObjectA_.seqBase_.seq_.size()) {
    stop = alignObjectA_.seqBase_.seq_.size();
  }
  for (uint32_t i = start; i < stop; ++i) {
  	//gap in alignObjectA, normally reference sequence
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectB_.seqBase_.seq_.substr(i, 1),
					alignObjectB_.seqBase_.qual_[i], true);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
        ++i;
        ++firstOffset;
      }
      /**@todo for now this is kept the real stop and start for
       *  reasons this function is used for but this might change
       *   or become a new function*/
			bool endGap = (newGap.startPos_ + newGap.size_ >= len(alignObjectA_)) || 0 == newGap.startPos_;
			if (!endGap) {
				gappedBasesInA += newGap.size_;
				++numberOfGaps;
				handleGapCountingInA(newGap);
				comp_.distances_.alignmentGaps_.insert(
						std::make_pair(newGap.startPos_, newGap));
			} else if (countEndGaps_) {
				++numberOfGaps;
				handleGapCountingInA(newGap);
				comp_.distances_.alignmentGaps_.insert(
						std::make_pair(newGap.startPos_, newGap));
			}
      continue;
    }
    //gap in alignObjectB, normally query sequence
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectA_.seqBase_.seq_.substr(i, 1),
					alignObjectA_.seqBase_.qual_[i], false);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-' ) {
        newGap.gapedSequence_.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
        ++i;
        ++secondOffset;
      }
      /**@todo for now this is kept the real stop and start for
       *  reasons this function is used for but this might change
       *   or become a new function*/
      bool endGap = (newGap.startPos_ + newGap.size_ >= alignObjectB_.seqBase_.seq_.length()) || 0 == newGap.startPos_;
			if (!endGap) {
				++numberOfGaps;
				gappedBasesInB += newGap.size_;
				handleGapCountingInB(newGap);
				comp_.distances_.alignmentGaps_.insert(
										std::make_pair(newGap.startPos_, newGap));
			}else if (countEndGaps_) {
				++numberOfGaps;
				handleGapCountingInB(newGap);
				comp_.distances_.alignmentGaps_.insert(
										std::make_pair(newGap.startPos_, newGap));
			}
      continue;
    }
    if ( 0 > parts_.scoring_.mat_[alignObjectA_.seqBase_.seq_[i]]
                                  [alignObjectB_.seqBase_.seq_[i]]) {
			auto firstK = getKmerPos(i - firstOffset, kMaps_.kLength_, objectA.seq_);
			auto secondK = getKmerPos(i - secondOffset, kMaps_.kLength_, objectB.seq_);
      if (usingQuality) {
  			/*if (objectA.checkQual(i - firstOffset, qScorePars_)
  					&& objectB.checkQual(i - secondOffset, qScorePars_)) {*/
  			/*if (objectB.checkQual(i - secondOffset, qScorePars_)) {*/
  			if ( (objectA.cnt_ >= 3 || objectA.checkQual(i - firstOffset, qScorePars_) )
  					&& objectB.checkQual(i - secondOffset, qScorePars_)) {

          if (checkKmer &&
              (kMaps_.isKmerLowFreq(firstK)
  								|| kMaps_.isKmerLowFreq(secondK))) {
            ++comp_.lowKmerMismatches_;
            comp_.distances_.lowKmerMismatches_.insert(std::make_pair(
                i,
                mismatch(
                    alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                    objectA.getLeadQual(i - firstOffset,qScorePars_.qualThresWindow_),
                    objectA.getTrailQual(i - firstOffset,qScorePars_.qualThresWindow_), i - firstOffset,
                    alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                    objectB.getLeadQual(i - secondOffset,qScorePars_.qualThresWindow_),
                    objectB.getTrailQual(i - secondOffset,qScorePars_.qualThresWindow_), i - secondOffset,
										kMaps_.kmersByPos_->getKmerFreq(secondK),
										kMaps_.kmersNoPos_->getKmerFreq(secondK))));
          } else {
            ++comp_.hqMismatches_;
            comp_.distances_.mismatches_.insert(std::make_pair(
                          i,
                          mismatch(
                              alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                              objectA.getLeadQual(i - firstOffset,qScorePars_.qualThresWindow_),
                              objectA.getTrailQual(i - firstOffset,qScorePars_.qualThresWindow_), i - firstOffset,
                              alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                              objectB.getLeadQual(i - secondOffset,qScorePars_.qualThresWindow_),
                              objectB.getTrailQual(i - secondOffset,qScorePars_.qualThresWindow_), i - secondOffset,
															kMaps_.kmersByPos_->getKmerFreq(secondK),
															kMaps_.kmersNoPos_->getKmerFreq(secondK))));
          }
        } else {
        	comp_.distances_.mismatches_.insert(std::make_pair(
              i,
              mismatch(
                  alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                  objectA.getLeadQual(i - firstOffset,qScorePars_.qualThresWindow_),
                  objectA.getTrailQual(i - firstOffset,qScorePars_.qualThresWindow_), i - firstOffset,
                  alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                  objectB.getLeadQual(i - secondOffset,qScorePars_.qualThresWindow_),
                  objectB.getTrailQual(i - secondOffset,qScorePars_.qualThresWindow_), i - secondOffset,
                  kMaps_.kmersByPos_->getKmerFreq(secondK),
									kMaps_.kmersNoPos_->getKmerFreq(secondK))));
          ++comp_.lqMismatches_;
        }
      } else {
        if (checkKmer &&
            (kMaps_.isKmerLowFreq(firstK)
								|| kMaps_.isKmerLowFreq(secondK))) {
          ++comp_.lowKmerMismatches_;
          comp_.distances_.lowKmerMismatches_.insert(std::make_pair(
              i,
              mismatch(
                  alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                  objectA.getLeadQual(i - firstOffset,qScorePars_.qualThresWindow_),
                  objectA.getTrailQual(i - firstOffset,qScorePars_.qualThresWindow_), i - firstOffset,
                  alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                  objectB.getLeadQual(i - secondOffset,qScorePars_.qualThresWindow_),
                  objectB.getTrailQual(i - secondOffset,qScorePars_.qualThresWindow_), i - secondOffset,
									kMaps_.kmersByPos_->getKmerFreq(secondK),
									kMaps_.kmersNoPos_->getKmerFreq(secondK))));
        } else {
          ++comp_.hqMismatches_;
          comp_.distances_.mismatches_.insert(std::make_pair(
                        i,
                        mismatch(
                            alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                            objectA.getLeadQual(i - firstOffset,qScorePars_.qualThresWindow_),
                            objectA.getTrailQual(i - firstOffset,qScorePars_.qualThresWindow_), i - firstOffset,
                            alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                            objectB.getLeadQual(i - secondOffset,qScorePars_.qualThresWindow_),
                            objectB.getTrailQual(i - secondOffset,qScorePars_.qualThresWindow_), i - secondOffset,
														kMaps_.kmersByPos_->getKmerFreq(secondK),
														kMaps_.kmersNoPos_->getKmerFreq(secondK))));
        }
      }
    } else {
      if (usingQuality) {
        if (doingMatchQuality) {
          if (objectA.checkQual(i - firstOffset, qScorePars_) &&
              objectB.checkQual(i - secondOffset, qScorePars_)) {
            comp_.highQualityMatches_++;
          } else {
            comp_.lowQualityMatches_++;
          }
        } else {
          comp_.highQualityMatches_++;
        }
      } else {
        ++comp_.highQualityMatches_;
      }
    }
  }

	uint32_t gappedLeftEndA = 0;
	uint32_t gappedRightEndA = 0;
	if (alignObjectA_.seqBase_.seq_.front() == '-') {
		gappedLeftEndA += countBeginChar(alignObjectA_.seqBase_.seq_);
	}
	if (alignObjectA_.seqBase_.seq_.back() == '-') {
		gappedRightEndA += countEndChar(alignObjectA_.seqBase_.seq_);
	}
	uint32_t gappedEndA = gappedRightEndA + gappedLeftEndA;

	uint32_t gappedLeftEndB = 0;
	uint32_t gappedRightEndB = 0;
	if (alignObjectB_.seqBase_.seq_.front() == '-') {
		gappedLeftEndB += countBeginChar(alignObjectB_.seqBase_.seq_);
	}
	if (alignObjectB_.seqBase_.seq_.back() == '-') {
		gappedRightEndB += countEndChar(alignObjectB_.seqBase_.seq_);
	}
	uint32_t gappedEndB = gappedLeftEndB + gappedRightEndB;

	uint32_t gappedEnd = gappedEndA + gappedEndB;

	//objectA (ref)
	comp_.distances_.ref_.covered_ = len(alignObjectA_) - gappedBasesInA - gappedEnd;
	comp_.distances_.ref_.coverage_ = static_cast<double>(comp_.distances_.ref_.covered_)/len(objectA);
	comp_.distances_.ref_.identities_ = comp_.highQualityMatches_ + comp_.lowQualityMatches_;
	comp_.distances_.ref_.identity_ = static_cast<double>(comp_.distances_.ref_.identities_)/len(objectA);
	//objectB (query)
	comp_.distances_.query_.covered_ = len(alignObjectB_) - gappedBasesInB - gappedEnd;
	comp_.distances_.query_.coverage_ = static_cast<double>(comp_.distances_.query_.covered_)/len(objectB);
	comp_.distances_.query_.identities_ = comp_.highQualityMatches_ + comp_.lowQualityMatches_;
	comp_.distances_.query_.identity_ = static_cast<double>(comp_.distances_.query_.identities_)/len(objectB);
	//of the overlapping alignment
	comp_.distances_.basesInAln_ = len(alignObjectA_) - gappedEnd;
  comp_.distances_.percentMatch_ = (comp_.highQualityMatches_ + comp_.lowQualityMatches_)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.percentMismatch_ = (comp_.hqMismatches_ + comp_.lqMismatches_)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.percentGaps_ = (gappedBasesInA + gappedBasesInB)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.overLappingEvents_ = comp_.highQualityMatches_
			+ comp_.lowQualityMatches_ + comp_.hqMismatches_ + comp_.lqMismatches_
			+ comp_.lowKmerMismatches_ + comp_.distances_.alignmentGaps_.size();
  if(comp_.distances_.overLappingEvents_ == 0){
  		comp_.distances_.eventBasedIdentity_ = 0;
  }else{
  		comp_.distances_.eventBasedIdentity_ = (comp_.highQualityMatches_ + comp_.lowQualityMatches_)/static_cast<double>(comp_.distances_.overLappingEvents_);
  }
  comp_.setEventBaseIdentityHq();
  comp_.refName_ = objectA.name_;
  comp_.queryName_ = objectB.name_;
  comp_.alnScore_ = parts_.score_;
  return comp_;
}








comparison aligner::compareAlignment(
    const seqInfo& objectA, const seqInfo& objectB,
    bool checkKmers) {

	resetAlignmentInfo();;
	uint32_t firstOffset = 0;
	uint32_t secondOffset = 0;
  uint32_t gappedBasesInA  = 0;
  uint32_t gappedBasesInB  = 0;
  uint32_t numberOfGaps = 0;
  for (uint32_t i = 0; i < len(alignObjectA_); ++i) {
  	//gap in alignObjectA, normally reference sequence
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
    	++numberOfGaps;
      ++firstOffset;
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectB_.seqBase_.seq_.substr(i, 1),
					alignObjectB_.seqBase_.qual_[i], true);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        //newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
        ++i;
        ++firstOffset;
      }
      bool endGap = (newGap.startPos_ + newGap.size_ >= len(alignObjectA_)) || 0 == newGap.startPos_;
			if (!endGap) {
				gappedBasesInA += newGap.size_;
				++numberOfGaps;
				handleGapCountingInA(newGap);
				comp_.distances_.alignmentGaps_.insert(
						std::make_pair(newGap.startPos_, newGap));
			} else if (countEndGaps_) {
				++numberOfGaps;
				handleGapCountingInA(newGap);
			}
      continue;
    }
    //gap in alignObjectB, normally query sequence
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectA_.seqBase_.seq_.substr(i, 1),
					alignObjectA_.seqBase_.qual_[i], false);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        //newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
        ++i;
        ++secondOffset;
      }
      bool endGap = (newGap.startPos_ + newGap.size_ >= alignObjectB_.seqBase_.seq_.length()) || 0 == newGap.startPos_;
			if (!endGap) {
				++numberOfGaps;
				gappedBasesInB += newGap.size_;
				handleGapCountingInB(newGap);
			} else if (countEndGaps_) {
				++numberOfGaps;
				handleGapCountingInB(newGap);
			}
      continue;
    }
    /**@todo consider changing this to be scoringMatrix_[A][B] < 0
     *  instead to allow custom scoring matrix scores to determine mismatch */
		if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
			/*if (objectA.checkQual(i - firstOffset, qScorePars_)
					&& objectB.checkQual(i - secondOffset, qScorePars_)) {*/
			/*if (objectB.checkQual(i - secondOffset, qScorePars_)) {*/
			if ( (objectA.cnt_ >= 3 || objectA.checkQual(i - firstOffset, qScorePars_) )
					&& objectB.checkQual(i - secondOffset, qScorePars_)) {
				auto firstK = getKmerPos(i - firstOffset, kMaps_.kLength_,
						objectA.seq_);
				auto secondK = getKmerPos(i - secondOffset, kMaps_.kLength_,
						objectB.seq_);
				if (checkKmers
						&& (kMaps_.isKmerLowFreq(firstK)
								|| kMaps_.isKmerLowFreq(secondK))) {
					++comp_.lowKmerMismatches_;
				} else {
					++comp_.hqMismatches_;
				}
			} else {
				++comp_.lqMismatches_;
			}
		} else {
			++comp_.highQualityMatches_;
		}
  }


	uint32_t gappedLeftEndA = 0;
	uint32_t gappedRightEndA = 0;
	if (alignObjectA_.seqBase_.seq_.front() == '-') {
		gappedLeftEndA += countBeginChar(alignObjectA_.seqBase_.seq_);
	}
	if (alignObjectA_.seqBase_.seq_.back() == '-') {
		gappedRightEndA += countEndChar(alignObjectA_.seqBase_.seq_);
	}
	uint32_t gappedEndA = gappedRightEndA + gappedLeftEndA;

	uint32_t gappedLeftEndB = 0;
	uint32_t gappedRightEndB = 0;
	if (alignObjectB_.seqBase_.seq_.front() == '-') {
		gappedLeftEndB += countBeginChar(alignObjectB_.seqBase_.seq_);
	}
	if (alignObjectB_.seqBase_.seq_.back() == '-') {
		gappedRightEndB += countEndChar(alignObjectB_.seqBase_.seq_);
	}
	uint32_t gappedEndB = gappedLeftEndB + gappedRightEndB;

	uint32_t gappedEnd = gappedEndA + gappedEndB;

	//objectA (ref)
	comp_.distances_.ref_.covered_ = len(alignObjectA_) - gappedBasesInA - gappedEnd;
	comp_.distances_.ref_.coverage_ = static_cast<double>(comp_.distances_.ref_.covered_)/len(objectA);
	comp_.distances_.ref_.identities_ = comp_.highQualityMatches_ + comp_.lowQualityMatches_;
	comp_.distances_.ref_.identity_ = static_cast<double>(comp_.distances_.ref_.identities_)/len(objectA);
	//objectB (query)
	comp_.distances_.query_.covered_ = len(alignObjectB_) - gappedBasesInB - gappedEnd;
	comp_.distances_.query_.coverage_ = static_cast<double>(comp_.distances_.query_.covered_)/len(objectB);
	comp_.distances_.query_.identities_ = comp_.highQualityMatches_ + comp_.lowQualityMatches_;
	comp_.distances_.query_.identity_ = static_cast<double>(comp_.distances_.query_.identities_)/len(objectB);
	//of the overlapping alignment
	comp_.distances_.basesInAln_ = len(alignObjectA_) - gappedEnd;
  comp_.distances_.percentMatch_ = (comp_.highQualityMatches_ + comp_.lowQualityMatches_)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.percentMismatch_ = (comp_.hqMismatches_ + comp_.lqMismatches_)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.percentGaps_ = (gappedBasesInA + gappedBasesInB)/static_cast<double>(comp_.distances_.basesInAln_);
  comp_.distances_.overLappingEvents_ = comp_.highQualityMatches_
			+ comp_.lowQualityMatches_ + comp_.hqMismatches_ + comp_.lqMismatches_
			+ comp_.lowKmerMismatches_ + comp_.distances_.alignmentGaps_.size();
  if(comp_.distances_.overLappingEvents_ == 0){
  	comp_.distances_.eventBasedIdentity_ = 0;
  }else{
  	comp_.distances_.eventBasedIdentity_ = (comp_.highQualityMatches_ + comp_.lowQualityMatches_)/static_cast<double>(comp_.distances_.overLappingEvents_);
  }
  comp_.setEventBaseIdentityHq();
  comp_.refName_ = objectA.name_;
  comp_.queryName_ = objectB.name_;
  comp_.alnScore_ = parts_.score_;
  return comp_;
}

void aligner::handleGapCountingInA(gap& currentGap) {
	if (!seqUtil::isHomopolyer(currentGap.gapedSequence_) || !weighHomopolymers_) {
		if (currentGap.size_ >= 3) {
			++comp_.largeBaseIndel_;
		} else if (currentGap.size_ == 2) {
			++comp_.twoBaseIndel_;
		} else if (currentGap.size_ == 1) {
			++comp_.oneBaseIndel_;
		}
	} else {
		uint32_t alnABases = 0;
		uint32_t alnBBases = 0;
		// forwards
		uint32_t cursor = 0;
		while ((currentGap.startPos_ + currentGap.size_ + cursor)
				< alignObjectA_.seqBase_.seq_.size()
				&& alignObjectA_.seqBase_.seq_[currentGap.startPos_ + currentGap.size_
						+ cursor] == currentGap.gapedSequence_[0]) {
			++alnABases;
			++cursor;
		}
		cursor = 1;
		// backwards
		while (cursor <= currentGap.startPos_
				&& alignObjectA_.seqBase_.seq_[currentGap.startPos_ - cursor]
						== currentGap.gapedSequence_[0]) {
			++alnABases;
			++cursor;
		}
		cursor = 0;
		// forwards
		while ((currentGap.startPos_ + cursor) < alignObjectB_.seqBase_.seq_.size()
				&& alignObjectB_.seqBase_.seq_[currentGap.startPos_ + cursor]
						== currentGap.gapedSequence_[0]) {
			++alnBBases;
			++cursor;
		}
		cursor = 1;
		// backwards
		while (cursor <= currentGap.startPos_
				&& alignObjectB_.seqBase_.seq_[currentGap.startPos_ - cursor]
						== currentGap.gapedSequence_[0]) {
			++alnBBases;
			++cursor;
		}
		//if it is a whole chuck of homopolymer missing, no weighting
		if (alnBBases == 0 || alnABases == 0) {
			if (currentGap.size_ >= 3) {
				++comp_.largeBaseIndel_;
			} else if (currentGap.size_ == 2) {
				++comp_.twoBaseIndel_;
			} else if (currentGap.size_ == 1) {
				++comp_.oneBaseIndel_;
			}
		} else {
			/*double currentScore = currentGap.size_
					/ ((alnABases * alignObjectA_.seqBase_.cnt_
							+ alnBBases * alignObjectB_.seqBase_.cnt_)
							/ (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));*/
			double currentScore = currentGap.size_ / static_cast<double>(alnABases + alnBBases);
			if (currentGap.size_ >= 3) {
				if (currentScore > 1) {
					++comp_.largeBaseIndel_;
				} else {
					comp_.largeBaseIndel_ += currentScore;
				}
			} else if (currentGap.size_ == 2) {
				comp_.twoBaseIndel_ += currentScore;
			} else if (currentGap.size_ == 1) {
				comp_.oneBaseIndel_ += currentScore;
			}
		}
	}
}
void aligner::handleGapCountingInB(gap& currentGap) {
	if (!seqUtil::isHomopolyer(currentGap.gapedSequence_) || !weighHomopolymers_) {
		if (currentGap.size_ >= 3) {
			++comp_.largeBaseIndel_;
		} else if (currentGap.size_ == 2) {
			++comp_.twoBaseIndel_;
		} else if (currentGap.size_ == 1) {
			++comp_.oneBaseIndel_;
		}
	} else {
		uint32_t alnABases = 0;
		uint32_t alnBBases = 0;
		// forwards
		uint32_t cursor = 0;
		while ((currentGap.startPos_ + currentGap.size_ + cursor)
				< alignObjectB_.seqBase_.seq_.size()
				&& alignObjectB_.seqBase_.seq_[currentGap.startPos_ + currentGap.size_
						+ cursor] == currentGap.gapedSequence_[0]) {
			++alnBBases;
			++cursor;
		}
		cursor = 1;
		// backwards
		while (cursor <= currentGap.startPos_
				&& alignObjectB_.seqBase_.seq_[currentGap.startPos_ - cursor]
						== currentGap.gapedSequence_[0]) {
			++alnBBases;
			++cursor;
		}
		cursor = 0;
		// forwards
		while ((currentGap.startPos_ + cursor) < alignObjectA_.seqBase_.seq_.size()
				&& alignObjectA_.seqBase_.seq_[currentGap.startPos_ + cursor]
						== currentGap.gapedSequence_[0]) {
			++alnABases;
			++cursor;
		}
		cursor = 1;
		// backwards
		while (cursor <= currentGap.startPos_
				&& alignObjectA_.seqBase_.seq_[currentGap.startPos_ - cursor]
						== currentGap.gapedSequence_[0]) {
			++alnABases;
			++cursor;
		}
		//if it is a whole chuck of homopolymer missing, no weighting
		if (alnBBases == 0 || alnABases == 0) {
			if (currentGap.size_ >= 3) {
				++comp_.largeBaseIndel_;
			} else if (currentGap.size_ == 2) {
				++comp_.twoBaseIndel_;
			} else if (currentGap.size_ == 1) {
				++comp_.oneBaseIndel_;
			}
		} else {
			/*double currentScore = currentGap.size_
					/ ((alnABases * alignObjectA_.seqBase_.cnt_
							+ alnBBases * alignObjectB_.seqBase_.cnt_)
							/ (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));*/
			double currentScore = currentGap.size_ / static_cast<double>(alnABases + alnBBases);
			if (currentGap.size_ >= 3) {
				if (currentScore > 1) {
					++comp_.largeBaseIndel_;
				} else {
					comp_.largeBaseIndel_ += currentScore;
				}
			} else if (currentGap.size_ == 2) {
				comp_.twoBaseIndel_ += currentScore;
			} else if (currentGap.size_ == 1) {
				comp_.oneBaseIndel_ += currentScore;
			}
		}
	}
}

void aligner::outPutParameterInfo(std::ostream& out) const {
  out << "numberOfOneIndel:" << comp_.oneBaseIndel_
      << " numberOfTwoIndel:" << comp_.twoBaseIndel_
      << " numberOfLargeGaps:" << comp_.largeBaseIndel_
      << " highQualityMismatch:" << comp_.hqMismatches_
      << " lowQualityMismatch:" << comp_.lqMismatches_
      << " lowKmerMismatch:" << comp_.lowKmerMismatches_ << std::endl;
}

// check for tandem repeat gaps
bool aligner::checkForTandemRepeatGap() {
  bool check = false;
  if (comp_.largeBaseIndel_ == 1) {
    std::map<uint32_t, gap>::iterator gapIter;
    for (gapIter = comp_.distances_.alignmentGaps_.begin(); gapIter != comp_.distances_.alignmentGaps_.end();
         ++gapIter) {
      if (gapIter->second.size_ >= 3) {
        std::string search;
        std::vector<TandemRepeat> gapTand =
            findTandemRepeatsInSequence(gapIter->second.gapedSequence_);
        if (gapTand.empty()) {
          search = gapIter->second.gapedSequence_;
        } else {
          search = gapTand[0].repeat_;
        }
        bool gapWithinTandem = false;
        if (alignObjectA_.seqBase_.seq_[gapIter->second.startPos_] == '-') {
          TandemRepeat secondTandems = findTandemRepeatOfStrInSequence(
              alignObjectB_.seqBase_.seq_, search);
          if (gapIter->second.startPos_ >= secondTandems.startPos_ &&
              gapIter->second.startPos_ + gapIter->second.size_ - 1 <=
                  secondTandems.stopPos_) {
            gapWithinTandem = true;
          }
          if (gapWithinTandem) {
            check = true;
          }
        } else if (alignObjectB_.seqBase_.seq_[gapIter->second.startPos_] ==
                   '-') {
          TandemRepeat secondTandems = findTandemRepeatOfStrInSequence(
              alignObjectA_.seqBase_.seq_, search);
          if (gapIter->second.startPos_ >= secondTandems.startPos_ &&
              gapIter->second.startPos_ + gapIter->second.size_ - 1 <=
                  secondTandems.stopPos_) {
            gapWithinTandem = true;
          }
          if (gapWithinTandem) {
            check = true;
          }
        } else {
          std::cout << "Error, gap start pos is not a gap in either sequences"
                    << std::endl;
        }
      }
    }
  }
  return check;
}

std::vector<TandemRepeat> aligner::findTandemRepeatsInSequence(
    const std::string& str, int match, int mismatch, int gap,
    int minimumAlignScore) {
  uint32_t sizeChecker = 2;
  std::vector<TandemRepeat> repeats;
  std::vector<TandemRepeat>::iterator repIter;
  bool foundTandem = false;
  int startPos = 0;
  int pos = 0;
  int numberOfRepeats = 0;
  int stopPos = pos;
  std::string tandem = "";
  while (sizeChecker < str.size() / 2 && !foundTandem) {
    pos = 0;
    bool keepSearching = true;
    numberOfRepeats = 0;
    int alignScore = 0;
    while (keepSearching && (pos + sizeChecker) < str.size()) {
      tandem = str.substr(pos, sizeChecker);
      bool homopolymer = true;
      for (uint32_t i = 0; i < tandem.size() - 1; ++i) {
        if (tandem[i] != tandem[i + 1]) {
          homopolymer = false;
          break;
        }
      }
      if (homopolymer) {
        ++pos;
        continue;
      }
      numberOfRepeats = 0;
      while (tandem == str.substr(pos + sizeChecker, sizeChecker)) {
        if (numberOfRepeats == 0) {
          startPos = pos;
        }
        pos += sizeChecker;
        if (numberOfRepeats == 0) {
          numberOfRepeats = 2;
        } else {
          ++numberOfRepeats;
        }
      }
      alignScore = (int)tandem.size() * (numberOfRepeats) * match;
      if (alignScore >= minimumAlignScore) {
        foundTandem = true;
      } else {
        foundTandem = false;
      }
      if (foundTandem) {
        stopPos = pos + sizeChecker - 1;
        bool alreadySmallerRepeat = false;
        for (repIter = repeats.begin(); repIter != repeats.end(); ++repIter) {
          TandemRepeat tempRep =
              findTandemRepeatOfStrInSequence(tandem, repIter->repeat_);
          if (tempRep.numberOfRepeats_ != 0) {
            alreadySmallerRepeat = true;
            break;
          }
        }
        if (!alreadySmallerRepeat) {
					repeats.emplace_back(tandem,
							numberOfRepeats,
							alignScore,
							startPos,
							startPos + numberOfRepeats * tandem.size());
        }
        pos = stopPos;
        foundTandem = false;
      }
      ++pos;
    }
    ++sizeChecker;
  }
  return repeats;
}

TandemRepeat aligner::findTandemRepeatOfStrInSequence(const std::string & str,
                                                      const std::string & tandem,
                                                      int match,
																									 int mismatch,
                                                      int gap,
                                                      int minimumAlignScore) {
  size_t position = str.find(tandem);
  if(std::string::npos != position){
    uint32_t startPosition = position;
    uint32_t numberOfRepeats = 1;
    while(position + tandem.size() < str.size() + 1 - tandem.size() &&
    		std::equal(tandem.begin(), tandem.end(), str.begin() + position + tandem.size()) ){
    		++numberOfRepeats;
    		position+=tandem.size();
    }
    int alignScore = tandem.size() * match * numberOfRepeats;
    if(alignScore >=minimumAlignScore){
    		return TandemRepeat(tandem,
    				numberOfRepeats,
						alignScore,
						startPosition,
						startPosition + tandem.size() * numberOfRepeats);
    }
  }
  return TandemRepeat("", 0, 0, 0, 0);
}

TandemRepeat aligner::findTandemRepeatOfStrInSequenceDegen(
    std::string str, std::string tandem, int match, int mismatch, int gap,
    int minimumAlignScore) {
  size_t sizeChecker = tandem.size();
  bool foundTandem = false;
  int startPos = 0;
  int pos = -(int)sizeChecker;
  int numberOfRepeats = 0;
  int alignScore = 0;
  int stopPos = 0;
  bool keepSearching = true;
  numberOfRepeats = 0;
  while (keepSearching && (pos + sizeChecker) < str.size()) {
    bool homopolymer = true;
    for (uint32_t i = 0; i < tandem.size() - 1; ++i) {
      if (tandem[i] != tandem[i + 1]) {
        homopolymer = false;
        break;
      }
    }
    if (homopolymer) {
      ++pos;
      continue;
    }
    numberOfRepeats = 0;
    while (tandem == str.substr(pos + sizeChecker, sizeChecker)) {
      if (numberOfRepeats == 1) {
        startPos = pos;
      }
      ++numberOfRepeats;
      pos += sizeChecker;
    }
    alignScore = tandem.size() * (numberOfRepeats) * match;
    if (alignScore >= minimumAlignScore) {
      foundTandem = true;
    } else {
      foundTandem = false;
    }
    if (foundTandem) {
      stopPos = pos + sizeChecker - 1;
      keepSearching = false;
    }
    ++pos;
  }
  if (!foundTandem) {
    return (TandemRepeat("", numberOfRepeats, alignScore, 0, 0));
  } else {
    return (
        TandemRepeat(tandem, numberOfRepeats, alignScore, startPos, stopPos));
  }
}
// limited right now for kmer checking
bool aligner::checkTwoEqualSeqs(const std::string& seq1,
                                const std::string& seq2,
                                int allowableMismatches) {
  int currentMismatches = 0;
  for (auto i : iter::range(seq1.size())) {
    if (seq1[i] != seq2[i]) {
      ++currentMismatches;
      if (currentMismatches > allowableMismatches) {
        return false;
      }
    }
  }
  return true;
}

bool aligner::checkTwoStringsDegen(
    const std::string& str1, const std::string& str2, int allowableMismatches,
    const substituteMatrix& scoringArrayIn) {
  if (str1.size() != str2.size()) {
    std::cout << "Strings should be same size, return false" << std::endl;
    return false;
  }
  int mismatches_ = 0;
  for (auto i : iter::range(str1.size())) {
    if (scoringArrayIn.mat_[str1[i]][str2[i]] < 0) {
      ++mismatches_;
    }
  }
  if (mismatches_ > allowableMismatches) {
    return false;
  } else {
    return true;
  }
}

void aligner::scoreAlignment(bool editTheSame) {
  parts_.score_ = 0;
  //editDistance_ = 0;
  for (uint32_t i = 0; i < len(alignObjectA_); ++i) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectB_.seqBase_.seq_.substr(i, 1),
					alignObjectB_.seqBase_.qual_[i], true);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
        ++i;
      }
      if (newGap.startPos_ + newGap.size_ >=
          alignObjectB_.seqBase_.seq_.length()) {
        if (editTheSame) {
          //editDistance_ += gapScores_.gapRightExtend_ +
                          // gapScores_.gapRightExtend_ * (newGap.size_ - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size_;
        }
        parts_.score_ -= parts_.gapScores_.gapRightRefOpen_;
        parts_.score_ -= parts_.gapScores_.gapRightRefExtend_ * (newGap.size_ - 1);


      } else if (newGap.startPos_ == 0) {
       	parts_.score_ -= parts_.gapScores_.gapLeftRefOpen_;
       	parts_.score_ -= parts_.gapScores_.gapLeftRefExtend_ * (newGap.size_ - 1);

        if (editTheSame) {
          //editDistance_ += gapScores_.gapLeftOpen_ +
                          // gapScores_.gapLeftExtend_ * (newGap.size_ - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size_;
        }
      } else {
      	parts_.score_ -= parts_.gapScores_.gapOpen_;
      	parts_.score_ -= parts_.gapScores_.gapExtend_ * (newGap.size_ - 1);
        if (editTheSame) {
          //editDistance_ +=
            //  gapScores_.gapOpen_ + gapScores_.gapExtend_ * (newGap.size_ - 1);
        } else {
          // std::cout << "adding gap" << std::endl;
          // std::cout << "gapSize: " << newGap.size_ << std::endl;
          //editDistance_ += newGap.size_;
        }
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
			gap newGap = gap(i, getSeqPosForAlnAPos(i), getSeqPosForAlnBPos(i),
					alignObjectA_.seqBase_.seq_.substr(i, 1),
					alignObjectA_.seqBase_.qual_[i], false);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
        ++i;
      }
      if (newGap.startPos_ + newGap.size_ >=
          alignObjectB_.seqBase_.seq_.length()) {
        if (editTheSame) {
          //editDistance_ += gapScores_.gapRightOpen_ +
                          // gapScores_.gapRightExtend_ * (newGap.size_ - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size_;
        }
        parts_.score_ -= parts_.gapScores_.gapRightQueryOpen_;
        parts_.score_ -= parts_.gapScores_.gapRightQueryExtend_ * (newGap.size_ - 1);
      } else if (newGap.startPos_ == 0) {
      	parts_.score_ -= parts_.gapScores_.gapLeftQueryOpen_;
      	parts_.score_ -= parts_.gapScores_.gapLeftQueryExtend_ * (newGap.size_ - 1);
        if (editTheSame) {
          //editDistance_ += gapScores_.gapLeftOpen_ +
                          // gapScores_.gapLeftExtend_ * (newGap.size_ - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size_;
        }
      } else {
      	parts_.score_ -= parts_.gapScores_.gapOpen_;
      	parts_.score_ -= parts_.gapScores_.gapExtend_ * (newGap.size_ - 1);
        if (editTheSame) {
          //editDistance_ +=
             // gapScores_.gapOpen_ + gapScores_.gapExtend_ * (newGap.size_ - 1);
        } else {
          // std::cout << "adding gap" << std::endl;
          // std::cout << "gapSize: " << newGap.size_ << std::endl;
          //editDistance_ += newGap.size_;
        }
      }
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      if (editTheSame) {
       // editDistance_ = -scoringArray_[alignObjectA_.seqBase_.seq_[i]]
         //                             [alignObjectB_.seqBase_.seq_[i]];
      } else {
        // std::cout << "adding mismatch" << std::endl;
        //++editDistance_;
      }
    }
    parts_.score_ += parts_.scoring_.mat_[alignObjectA_.seqBase_.seq_[i]]
                           [alignObjectB_.seqBase_.seq_[i]];
  }
  comp_.alnScore_ = parts_.score_;
}

void aligner::noAlignSetAndScore(const seqInfo& objectA,
		const seqInfo& objectB) {
	alignObjectA_.seqBase_ = objectA;
	alignObjectB_.seqBase_ = objectB;

	if (alignObjectA_.seqBase_.seq_.size() < alignObjectB_.seqBase_.seq_.size()) {
		alignObjectA_.seqBase_.append(
				std::string('-',
						alignObjectB_.seqBase_.seq_.size()
								- alignObjectA_.seqBase_.seq_.size()), 0);
	} else if (alignObjectA_.seqBase_.seq_.size()
			> alignObjectB_.seqBase_.seq_.size()) {
		alignObjectB_.seqBase_.append(
				std::string('-',
						alignObjectA_.seqBase_.seq_.size()
								- alignObjectB_.seqBase_.seq_.size()), 0);
	}

	scoreAlignment(false);

}


void aligner::processAlnInfoInputNoCheck(const std::string& alnInfoDirName, bool verbose){
	//std::cout << __FILE__ << ' '<< __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if (alnInfoDirName != "") {
		if(bib::files::bfs::exists(alnInfoDirName)){
			//std::cout << __FILE__ << ' '<< __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			alnHolder_.read(alnInfoDirName, verbose);
		}
	}
}

void aligner::processAlnInfoOutputNoCheck(const std::string& outAlnInfoDirName, bool verbose){
	if ("" != outAlnInfoDirName) {
		if(!bib::files::bfs::exists(outAlnInfoDirName)){
			bib::files::makeDir(bib::files::MkdirPar(outAlnInfoDirName));
		}
		alnHolder_.write(outAlnInfoDirName, verbose);
	}
}

void aligner::processAlnInfoInput(const std::string& alnInfoDirName, bool verbose) {
	if (alnInfoDirName != "") {
		alignment::alnCacheDirSearchLock.lock();
		auto fullPath = bib::files::normalize(alnInfoDirName).string();
		auto search = alignment::alnCacheDirLocks.find(fullPath);
		if (search == alignment::alnCacheDirLocks.end()) {
			alignment::alnCacheDirLocks.emplace(fullPath, std::make_unique<std::shared_timed_mutex>());
		}
		auto realSearch = alignment::alnCacheDirLocks.find(fullPath);
		{
			std::shared_lock<std::shared_timed_mutex> lock(*(realSearch->second));
			alignment::alnCacheDirSearchLock.unlock();
			processAlnInfoInputNoCheck(alnInfoDirName, verbose);
		}
	}
}

void aligner::processAlnInfoOutput(const std::string& outAlnInfoDirName, bool verbose) {
	if (outAlnInfoDirName != "") {
		alignment::alnCacheDirSearchLock.lock();
		auto fullPath = bib::files::normalize(outAlnInfoDirName).string();
		auto search = alignment::alnCacheDirLocks.find(fullPath);
		if (search == alignment::alnCacheDirLocks.end()) {
			alignment::alnCacheDirLocks.emplace(fullPath, std::make_unique<std::shared_timed_mutex>());
		}
		auto realSearch = alignment::alnCacheDirLocks.find(fullPath);
		{
			std::unique_lock<std::shared_timed_mutex> lock(*(realSearch->second));
			alignment::alnCacheDirSearchLock.unlock();
			processAlnInfoOutputNoCheck(outAlnInfoDirName, verbose);
		}
	}
}

void aligner::setDefaultQualities() {
	qScorePars_.primaryQual_ = 20;
	qScorePars_.secondaryQual_ = 15;
	qScorePars_.qualThresWindow_ = 5;
};

void aligner::setGeneralScorring(int32_t generalMatch, int32_t generalMismatch){
	parts_.scoring_ = substituteMatrix(generalMatch, generalMismatch);
}


}  // namespace bib
