#include "aligner.hpp"
//
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

namespace bibseq {
struct kmerPos{
	std::string kmer_;
	uint32_t kPos_;
};

inline kmerPos getKmerPos(uint32_t realPos, uint32_t kLength, std::string str) {
	std::string kmer = "";
	uint32_t kPos = realPos;
	if ((kLength / 2) > realPos) {
		kPos = 0;
		kmer = str.substr(0, kLength);
	} else if (realPos + kLength / 2 >= str.size()) {
		kmer = str.substr(str.size() - kLength, kLength);
		kPos = str.size() - kLength;
	} else {
		kmer = str.substr(realPos - kLength / 2, kLength);
		kPos = realPos - kLength / 2;
	}
	return kmerPos { kmer, kPos };
}

void aligner::resetCounts() {
  comp_.resetCounts();
}

void aligner::resetAlignmentInfo() {
  // the mismatches_
  mismatches_.clear();
  lowKmerMismatches_.clear();
  // the indels
  alignmentGaps_.clear();
  resetCounts();
}
void aligner::setStartingParamters(const kmerMaps& inKmaps, uint32_t primaryQuality,
		uint32_t secondaryQuality, uint32_t qualThresholdWindow,
                          bool countEndGaps) {
  setDefaultQualities();
  countEndGaps_ = countEndGaps;
  setQual(primaryQuality, secondaryQuality, qualThresholdWindow);
  kMaps_ = inKmaps;
}

void aligner::setQual(int pQual, int sQual, int qualThresholdWindow) {
  primaryQual_ = pQual;
  secondaryQual_ = sQual;
  qualThresWindow_ = qualThresholdWindow;
}


void aligner::profilePrimerAlignment(const seqInfo& objectA,
                            const seqInfo& objectB,
                            bool weighHomopolymers){
  resetAlignmentInfo();
  int firstOffset = 0;
  int secondOffset = 0;
  int sumOfGappedEndsInA = 0;
  int sumOfRegularGappedInA = 0;
  int sumOfGappedEndsInB = 0;
  int sumOfRegularGappedInB = 0;
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); ++i) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i], true);
      //newGap.ref = true;
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
        ++i;
        ++firstOffset;
      }
      if ((((newGap.startPos_ + newGap.size_) >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos_ == 0) &&
          !countEndGaps_) {
        sumOfGappedEndsInA += newGap.size_;
      } else {
        handleGapCountingInA(newGap, weighHomopolymers);
        sumOfRegularGappedInA += newGap.size_;
        alignmentGaps_.insert(std::make_pair(newGap.startPos_, newGap));
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i], false);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
        ++i;
        ++secondOffset;
      }
      if (((newGap.startPos_ + newGap.size_) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos_ == 0) &&
          !countEndGaps_) {
        sumOfGappedEndsInB += newGap.size_;
      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
        sumOfRegularGappedInB += newGap.size_;
        alignmentGaps_.insert(std::make_pair(newGap.startPos_, newGap));
      }
      continue;
    }
    // alignObjectA_.seqBase_.seq_[i]!=alignObjectB_.seqBase_.seq_[i]
    if (0 > parts_.scoring_.mat_[alignObjectA_.seqBase_.seq_[i]]
                         [alignObjectB_.seqBase_.seq_[i]]) {
      ++comp_.hqMismatches_;
      mismatches_.insert(std::make_pair(
                  i,
                  mismatch(
                      alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                      objectA.getLeadQual(i - firstOffset),
                      objectA.getTrailQual(i - firstOffset), i - firstOffset,
                      alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                      objectB.getLeadQual(i - secondOffset),
                      objectB.getTrailQual(i - secondOffset), i - secondOffset,
                      0,0)));
    } else {
      ++comp_.highQualityMatches_;
    }
  }

  double coverageArea = alignObjectA_.seqBase_.seq_.size() - sumOfGappedEndsInA - sumOfGappedEndsInB;
  //std::cout << "cov area"  << coverageArea << std::endl;
  //std::cout << "sumOfGappedEndsInA area"  << sumOfGappedEndsInA << std::endl;
  //std::cout << "sumOfGappedEndsInB area"  << sumOfGappedEndsInB << std::endl;
  comp_.distances_.refCoverage_ = (coverageArea - sumOfRegularGappedInA)/objectA.seq_.size();
  comp_.distances_.queryCoverage_ = (coverageArea - sumOfRegularGappedInB)/objectB.seq_.size();
  comp_.distances_.percentIdentity_ = (comp_.highQualityMatches_ + comp_.lowQualityMatches_) / coverageArea;
  comp_.distances_.percentMismatch_ = (comp_.hqMismatches_ + comp_.lowKmerMismatches_ + comp_.lqMismatches_)/coverageArea;
  comp_.distances_.percentageGaps_ = (sumOfRegularGappedInA + sumOfRegularGappedInB)/ coverageArea;
  comp_.distances_.identity_ = comp_.distances_.percentIdentity_;
  comp_.distances_.ownDistance_ = comp_.hqMismatches_;
  comp_.distances_.ownGapDistance_ = 0;
  for (const auto& g : alignmentGaps_) {
  	comp_.distances_.ownGapDistance_ += g.second.homoploymerScore_ * g.second.size_;
  	comp_.distances_.ownDistance_ += g.second.homoploymerScore_ * g.second.size_;
  }
  comp_.distances_.eventBasedIdentity_ = 1.0 - ((alignmentGaps_.size() + comp_.hqMismatches_)/
  		static_cast<double>(comp_.hqMismatches_ + comp_.highQualityMatches_ + comp_.hqMismatches_));
  return;
}
// profile primer alignment
void aligner::profilePrimerAlignment(const baseReadObject& objectA,
                                     const baseReadObject& objectB,
                                     bool weighHomopolymers) {
	profilePrimerAlignment(objectA.seqBase_, objectB.seqBase_, weighHomopolymers);
}

uint32_t aligner::getAlignPosForSeqAPos(uint32_t seqAPos) {
  uint32_t offSet = 0;
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.size(); ++i) {
    if ((i - offSet) == seqAPos) {
      return i;
    }
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++offSet;
    }
  }
  return alignObjectA_.seqBase_.seq_.size();
}
uint32_t aligner::getAlignPosForSeqBPos(uint32_t seqBPos) {
  uint32_t offSet = 0;
  for (uint32_t i = 0; i < alignObjectB_.seqBase_.seq_.size(); ++i) {
    if ((i - offSet) == seqBPos) {
      return i;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++offSet;
    }
  }
  return alignObjectB_.seqBase_.seq_.size();
}
void aligner::profileAlignment(const baseReadObject& objectA,
                               const baseReadObject& objectB, int kLength,
                               bool kmersByPosition, bool checkKmer,
                               bool usingQuality, bool doingMatchQuality,
                               bool weighHomopolyer, uint32_t start,
                               uint32_t stop) {
  profileAlignment(objectA.seqBase_, objectB.seqBase_, kLength, kmersByPosition,
                   checkKmer, usingQuality, doingMatchQuality, weighHomopolyer,
                   start, stop);
}

void aligner::profileAlignment(const seqInfo& objectA, const seqInfo& objectB,
                               int kLength, bool kmersByPosition,
                               bool checkKmer, bool usingQuality,
                               bool doingMatchQuality, bool weighHomopolyer,
                               uint32_t start, uint32_t stop) {
  resetAlignmentInfo();
  //editDistance_ = 0;
  int firstOffset = 0;
  int secondOffset = 0;
  int sumOfGappedEndsInA = 0;
  int sumOfRegularGappedInA = 0;
  // std::cout << "start: " << start << std::endl;
  if (start != 0) {
    for (uint32_t i = 0; i < start; ++i) {
      // std::cout << i <<"/" << alignObjectA_.seqBase_.seq_.length() <<
      // std::endl;
      if (alignObjectA_.seqBase_.seq_[i] == '-') {
        ++firstOffset;
        gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                         alignObjectB_.seqBase_.qual_[i], true);
        while (alignObjectA_.seqBase_.seq_[i + 1] == '-' && (i + 1) != start) {
          newGap.gapedSequence_.append(
              alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
          ++newGap.size_;
          newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
          ++i;
          ++firstOffset;
        }
        // firstOffset += newGap.size_;
        if (((newGap.startPos_ + newGap.size_ >=
              alignObjectA_.seqBase_.seq_.length()) ||
             newGap.startPos_ == 0) &&
            !countEndGaps_) {
          sumOfGappedEndsInA += newGap.size_;
        } else {
          sumOfRegularGappedInA += newGap.size_;
        }
        continue;
      }
      if (alignObjectB_.seqBase_.seq_[i] == '-') {
        ++secondOffset;
        gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                         alignObjectA_.seqBase_.qual_[i], false);
        while (alignObjectB_.seqBase_.seq_[i + 1] == '-' && (i + 1) != start) {
          newGap.gapedSequence_.append(
              alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
          ++newGap.size_;
          newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
          ++i;
          ++secondOffset;
        }
        if (((newGap.startPos_ + newGap.size_) >=
                 alignObjectB_.seqBase_.seq_.length() ||
             newGap.startPos_ == 0) &&
            !countEndGaps_) {

        } else {
        }
        continue;
      }
    }
  }
  if (stop == 0) {
    stop = alignObjectA_.seqBase_.seq_.size();
  }
  for (uint32_t i = start; i < stop; ++i) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i],true);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
        ++i;
        ++firstOffset;
      }
      // firstOffset += newGap.size_;
      if (((newGap.startPos_ + newGap.size_ >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos_ == 0) &&
          !countEndGaps_) {
        sumOfGappedEndsInA += newGap.size_;
      } else {
        sumOfRegularGappedInA += newGap.size_;
        handleGapCountingInA(newGap, weighHomopolyer);
        alignmentGaps_.insert(std::make_pair(newGap.startPos_, newGap));
        //editDistance_ += newGap.size_;
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i],false);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size_;
        newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
        ++i;
        ++secondOffset;
      }
      if (((newGap.startPos_ + newGap.size_) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos_ == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolyer);
        alignmentGaps_.insert(std::make_pair(newGap.startPos_, newGap));
        //editDistance_ += newGap.size_;
      }
      continue;
    }

    if ( 0 > parts_.scoring_.mat_[alignObjectA_.seqBase_.seq_[i]]
                                  [alignObjectB_.seqBase_.seq_[i]]) {
			auto firstK = getKmerPos(i - firstOffset, kMaps_.kLength_, objectA.seq_);
			auto secondK = getKmerPos(i - secondOffset, kMaps_.kLength_, objectB.seq_);

      if (usingQuality) {
        if (objectA.checkQual(i - firstOffset, primaryQual_, secondaryQual_,
                              qualThresWindow_) &&
            objectB.checkQual(i - secondOffset, primaryQual_, secondaryQual_,
                              qualThresWindow_)) {
          if (checkKmer &&
              (kMaps_.isKmerLowFrequency(firstK.kmer_, firstK.kPos_, kmersByPosition,
                                         kMaps_.runCutOff_) ||
               kMaps_.isKmerLowFrequency(secondK.kmer_, secondK.kPos_, kmersByPosition,
                                         kMaps_.runCutOff_))) {
            ++comp_.lowKmerMismatches_;
            lowKmerMismatches_.insert(std::make_pair(
                i,
                mismatch(
                    alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                    objectA.getLeadQual(i - firstOffset),
                    objectA.getTrailQual(i - firstOffset), i - firstOffset,
                    alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                    objectB.getLeadQual(i - secondOffset),
                    objectB.getTrailQual(i - secondOffset), i - secondOffset,
                    kMaps_.kmersByPos_[secondK.kPos_][secondK.kmer_].readCnt_,
                    kMaps_.kmersAnyWhere_[secondK.kmer_].readCnt_)));
          } else {
            ++comp_.hqMismatches_;
            mismatches_.insert(std::make_pair(
                          i,
                          mismatch(
                              alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                              objectA.getLeadQual(i - firstOffset),
                              objectA.getTrailQual(i - firstOffset), i - firstOffset,
                              alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                              objectB.getLeadQual(i - secondOffset),
                              objectB.getTrailQual(i - secondOffset), i - secondOffset,
                              kMaps_.kmersByPos_[secondK.kPos_][secondK.kmer_].readCnt_,
                              kMaps_.kmersAnyWhere_[secondK.kmer_].readCnt_)));
          }

        } else {
          mismatches_.insert(std::make_pair(
              i,
              mismatch(
                  alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                  objectA.getLeadQual(i - firstOffset),
                  objectA.getTrailQual(i - firstOffset), i - firstOffset,
                  alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                  objectB.getLeadQual(i - secondOffset),
                  objectB.getTrailQual(i - secondOffset), i - secondOffset,
                  kMaps_.kmersByPos_[secondK.kPos_][secondK.kmer_].readCnt_,
                  kMaps_.kmersAnyWhere_[secondK.kmer_].readCnt_)));
          ++comp_.lqMismatches_;
        }
      } else {
        if (checkKmer &&
            (kMaps_.isKmerLowFrequency(firstK.kmer_, firstK.kPos_, kmersByPosition,
                                       kMaps_.runCutOff_) ||
             kMaps_.isKmerLowFrequency(secondK.kmer_, secondK.kPos_, kmersByPosition,
                                       kMaps_.runCutOff_))) {
          ++comp_.lowKmerMismatches_;
          lowKmerMismatches_.insert(std::make_pair(
              i,
              mismatch(
                  alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                  objectA.getLeadQual(i - firstOffset),
                  objectA.getTrailQual(i - firstOffset), i - firstOffset,
                  alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                  objectB.getLeadQual(i - secondOffset),
                  objectB.getTrailQual(i - secondOffset), i - secondOffset,
                  kMaps_.kmersByPos_[secondK.kPos_][secondK.kmer_].readCnt_,
                  kMaps_.kmersAnyWhere_[secondK.kmer_].readCnt_)));
        } else {
          ++comp_.hqMismatches_;
          mismatches_.insert(std::make_pair(
                        i,
                        mismatch(
                            alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                            objectA.getLeadQual(i - firstOffset),
                            objectA.getTrailQual(i - firstOffset), i - firstOffset,
                            alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                            objectB.getLeadQual(i - secondOffset),
                            objectB.getTrailQual(i - secondOffset), i - secondOffset,
                            kMaps_.kmersByPos_[secondK.kPos_][secondK.kmer_].readCnt_,
                            kMaps_.kmersAnyWhere_[secondK.kmer_].readCnt_)));
        }
      }

    }
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      if (usingQuality) {
        if (doingMatchQuality) {
          if (objectA.checkQual(i - firstOffset, primaryQual_, secondaryQual_,
                                qualThresWindow_) &&
              objectB.checkQual(i - secondOffset, primaryQual_, secondaryQual_,
                                qualThresWindow_)) {
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

	uint32_t overLappingBases = comp_.highQualityMatches_
			+ comp_.lowQualityMatches_ + comp_.hqMismatches_ + comp_.lqMismatches_
			+ comp_.lowKmerMismatches_;
	uint32_t gapedEndA = 0;
	if(alignObjectA_.seqBase_.seq_.back() == '-'){
		gapedEndA+=countEndChar(alignObjectA_.seqBase_.seq_);
	}
	if(alignObjectA_.seqBase_.seq_.front() == '-'){
		gapedEndA+=countBeginChar(alignObjectA_.seqBase_.seq_);
	}
	uint32_t gapedEndB = 0;
	if(alignObjectB_.seqBase_.seq_.back() == '-'){
		gapedEndB+=countEndChar(alignObjectB_.seqBase_.seq_);
	}
	if(alignObjectB_.seqBase_.seq_.front() == '-'){
		gapedEndB+=countBeginChar(alignObjectB_.seqBase_.seq_);
	}

	comp_.distances_.percentIdentity_ = static_cast<double>(comp_.highQualityMatches_
			+ comp_.lowQualityMatches_)/ (static_cast<double>(objectB.seq_.size()) - gapedEndA);
	comp_.distances_.queryCoverage_ = static_cast<double> (objectB.seq_.size() - gapedEndA) / objectB.seq_.size();
	comp_.distances_.percentageGaps_ =static_cast<double>(alignObjectB_.seqBase_.seq_.size() - objectB.seq_.size() - gapedEndA)/ alignObjectB_.seqBase_.seq_.size();
  comp_.distances_.eventBasedIdentity_ = 1.0 - ((comp_.hqMismatches_ + comp_.lqMismatches_ + comp_.lowKmerMismatches_ + alignmentGaps_.size())/
  		static_cast<double>(overLappingBases));
}






comparison aligner::compareAlignment(
    const baseReadObject& objectA, const baseReadObject& objectB,
    const runningParameters& runParams, bool checkKmers, bool kmersByPosition,
    bool weighHomopolymers) {
	return compareAlignment(objectA.seqBase_, objectB.seqBase_, runParams, checkKmers, kmersByPosition, weighHomopolymers);
}

comparison aligner::compareAlignment(
    const seqInfo& objectA, const seqInfo& objectB,
    const runningParameters& runParams, bool checkKmers, bool kmersByPosition,
    bool weighHomopolymers) {

	resetAlignmentInfo();;
  int firstOffset = 0;
  int secondOffset = 0;

  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); ++i) {
  	//std::cout << i << std::endl;
    // std::cout <<"cab2.1" << std::endl;
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      firstOffset++;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i],true);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size_++;
        newGap.qualities_.emplace_back(alignObjectB_.seqBase_.qual_[i + 1]);
        ++i;
      }
      if (((newGap.startPos_ + newGap.size_ >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos_ == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInA(newGap, weighHomopolymers);
      }
      // std::cout <<"cab2.2" << std::endl;
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      // std::cout <<"cab3.1" << std::endl;
      secondOffset++;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i],false);
      // std::cout <<"cab3.2" << std::endl;
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence_.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size_++;
        newGap.qualities_.emplace_back(alignObjectA_.seqBase_.qual_[i + 1]);
        ++i;
      }
      if (((newGap.startPos_ + newGap.size_) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos_ == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
      }
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
			if (objectA.checkQual(i - firstOffset, primaryQual_,
					secondaryQual_, qualThresWindow_)
					&& objectB.checkQual(i - secondOffset, primaryQual_,
							secondaryQual_, qualThresWindow_)) {
				auto firstK = getKmerPos(i - firstOffset, kMaps_.kLength_, objectA.seq_);
				auto secondK = getKmerPos(i - secondOffset, kMaps_.kLength_, objectB.seq_);
				if (checkKmers
						&& (kMaps_.isKmerLowFrequency(firstK.kmer_, firstK.kPos_, kmersByPosition,
								kMaps_.runCutOff_)
								|| kMaps_.isKmerLowFrequency(secondK.kmer_, secondK.kPos_,
										kmersByPosition, kMaps_.runCutOff_))) {
					++comp_.lowKmerMismatches_;
				} else {
					comp_.hqMismatches_++;
				}
			} else {
				comp_.lqMismatches_++;
			}
    }
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      ++comp_.highQualityMatches_;
    }
  }

	uint32_t overLappingBases = comp_.highQualityMatches_
			+ comp_.lowQualityMatches_ + comp_.hqMismatches_ + comp_.lqMismatches_
			+ comp_.lowKmerMismatches_;
	uint32_t gapedEndA = 0;
	if(alignObjectA_.seqBase_.seq_.back() == '-'){
		gapedEndA+=countEndChar(alignObjectA_.seqBase_.seq_);
	}
	if(alignObjectA_.seqBase_.seq_.front() == '-'){
		gapedEndA+=countBeginChar(alignObjectA_.seqBase_.seq_);
	}
	uint32_t gapedEndB = 0;
	if(alignObjectB_.seqBase_.seq_.back() == '-'){
		gapedEndB+=countEndChar(alignObjectB_.seqBase_.seq_);
	}
	if(alignObjectB_.seqBase_.seq_.front() == '-'){
		gapedEndB+=countBeginChar(alignObjectB_.seqBase_.seq_);
	}

	comp_.distances_.percentIdentity_ = static_cast<double>(comp_.highQualityMatches_
			+ comp_.lowQualityMatches_)/ (static_cast<double>(objectB.seq_.size()) - gapedEndA);
	comp_.distances_.queryCoverage_ = static_cast<double> (objectB.seq_.size() - gapedEndA) / objectB.seq_.size();
	comp_.distances_.percentageGaps_ =static_cast<double>(alignObjectB_.seqBase_.seq_.size() - objectB.seq_.size() - gapedEndA)/ alignObjectB_.seqBase_.seq_.size();
  comp_.distances_.eventBasedIdentity_ = 1.0 - ((comp_.hqMismatches_ + comp_.lqMismatches_ + comp_.lowKmerMismatches_ + alignmentGaps_.size())/
  		static_cast<double>(overLappingBases));
  return comp_;
}



void aligner::handleGapCountingInA(gap& currentGap, bool weighHomopolymers) {

  if (!seqUtil::isHomopolyer(currentGap.gapedSequence_) || !weighHomopolymers) {
    if (currentGap.size_ >= 3) {
      ++comp_.largeBaseIndel_;
      currentGap.homoploymerScore_ = 1;
    } else if (currentGap.size_ == 2) {
      ++comp_.twoBaseIndel_;
      currentGap.homoploymerScore_ = 1;
    } else if (currentGap.size_ == 1) {
      ++comp_.oneBaseIndel_;
      currentGap.homoploymerScore_ = 1;
    }
  } else {
    double firstBases = 0.00;
    double secondBases = 0.00;
    // bases++;
    // forwards
    int cursor = 0;
    while ((currentGap.startPos_ + currentGap.size_ + cursor) <
               alignObjectA_.seqBase_.seq_.size() &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos_ + currentGap.size_ +
                                       cursor] == currentGap.gapedSequence_[0]) {
      ++firstBases;
      ++cursor;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos_ - cursor) >= 0 &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos_ - cursor] ==
               currentGap.gapedSequence_[0]) {
      ++firstBases;
      ++cursor;
    }
    cursor = 0;
    // forwards
    while ((currentGap.startPos_ + cursor) <
               alignObjectB_.seqBase_.seq_.size() &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos_ + cursor] ==
               currentGap.gapedSequence_[0]) {
      ++secondBases;
      ++cursor;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos_ - cursor) >= 0 &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos_ - cursor] ==
               currentGap.gapedSequence_[0]) {
      ++secondBases;
      ++cursor;
    }
    // std::cout<<"secondBases: "<<secondBases<<" firstBases:
    // "<<firstBases<<std::endl;
    // std::cout<<"secondBasesSize: "<<alignObjectB_.seqBase_.cnt_<<"
    // firstBasesSize: "<<alignObjectA_.seqBase_.cnt_<<std::endl;
    // std::cout<<"size: "<<currentGap.size_<<std::endl<<std::endl;
    if (secondBases == 0 || firstBases == 0) {
      if (currentGap.size_ >= 3) {
        ++comp_.largeBaseIndel_;
        currentGap.homoploymerScore_ = 1;
      } else if (currentGap.size_ == 2) {
        ++comp_.twoBaseIndel_;
        currentGap.homoploymerScore_ = 1;
      } else if (currentGap.size_ == 1) {
        ++comp_.oneBaseIndel_;
        currentGap.homoploymerScore_ = 1;
      }
      // //std::cout<<"mark 1"<<std::endl;
    } else {
      // //std::cout<<"mark 2"<<std::endl;

      if (secondBases < firstBases) {
        // //std::cout<<"mark 3"<<std::endl;
        if (currentGap.size_ >= 3) {
          double currentScore =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++comp_.largeBaseIndel_;
            currentGap.homoploymerScore_ = 1;
          } else {
            comp_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore_ = currentScore;
          }

        } else if (currentGap.size_ == 2) {
          comp_.twoBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size_ == 1) {
          comp_.oneBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
      } else {
        // //std::cout<<"mark 4"<<std::endl;
        if (currentGap.size_ >= 3) {
          double currentScore =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++comp_.largeBaseIndel_;
            currentGap.homoploymerScore_ = currentScore;
          } else {
            comp_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore_ = currentScore;
          }

        } else if (currentGap.size_ == 2) {
          comp_.twoBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size_ == 1) {
          comp_.oneBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
      }
    }
  }
}
void aligner::handleGapCountingInB(gap& currentGap, bool weighHomopolymers) {

  if (!seqUtil::isHomopolyer(currentGap.gapedSequence_) || !weighHomopolymers) {
    if (currentGap.size_ >= 3) {
      ++comp_.largeBaseIndel_;
      currentGap.homoploymerScore_ = 1;
    } else if (currentGap.size_ == 2) {
      ++comp_.twoBaseIndel_;
      currentGap.homoploymerScore_ = 1;
    } else if (currentGap.size_ == 1) {
      ++comp_.oneBaseIndel_;
      currentGap.homoploymerScore_ = 1;
    }
  } else {
    double firstBases = 0.00;
    double secondBases = 0.00;
    // bases++;
    // forwards
    int cursor = 0;
    while ((currentGap.startPos_ + cursor) <
               alignObjectA_.seqBase_.seq_.size() &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos_ + cursor] ==
               currentGap.gapedSequence_[0]) {
      ++firstBases;
      ++cursor;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos_ - cursor) >= 0 &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos_ - cursor] ==
               currentGap.gapedSequence_[0]) {
      ++firstBases;
      ++cursor;
    }
    // forwards
    cursor = 0;
    while ((currentGap.startPos_ + currentGap.size_ + cursor) <
               alignObjectB_.seqBase_.seq_.size() &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos_ + currentGap.size_ +
                                       cursor] == currentGap.gapedSequence_[0]) {
      ++secondBases;
      ++cursor;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos_ - cursor) >= 0 &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos_ - cursor] ==
               currentGap.gapedSequence_[0]) {
      ++secondBases;
      ++cursor;
    }

    // std::cout<<"secondBases: "<<secondBases<<" firstBases:
    // "<<firstBases<<std::endl;
    // std::cout<<"secondBasesSize: "<<alignObjectB_.seqBase_.cnt_<<"
    // firstBasesSize: "<<alignObjectA_.seqBase_.cnt_<<std::endl;
    // std::cout<<"size: "<<currentGap.size_<<std::endl<<std::endl;

    if (secondBases == 0 || firstBases == 0) {
      if (currentGap.size_ >= 3) {
        ++comp_.largeBaseIndel_;
        currentGap.homoploymerScore_ = 1;
      } else if (currentGap.size_ == 2) {
        ++comp_.twoBaseIndel_;
        currentGap.homoploymerScore_ = 1;
      } else if (currentGap.size_ == 1) {
        ++comp_.oneBaseIndel_;
        currentGap.homoploymerScore_ = 1;
      }
      // //std::cout<<"mark 1"<<std::endl;
    } else {
      // //std::cout<<"mark 2"<<std::endl;

      if (secondBases < firstBases) {
        // //std::cout<<"mark 3"<<std::endl;
        if (currentGap.size_ >= 3) {
          double currentScore =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++comp_.largeBaseIndel_;
            currentGap.homoploymerScore_ = 1;
          } else {
            comp_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore_ = currentScore;
          }
        } else if (currentGap.size_ == 2) {
          comp_.twoBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size_ == 1) {
          comp_.oneBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
      } else {
        // //std::cout<<"mark 4"<<std::endl;
        if (currentGap.size_ >= 3) {
          double currentScore =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++comp_.largeBaseIndel_;
            currentGap.homoploymerScore_ = 1;
          } else {
            comp_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore_ = currentScore;
          }
        } else if (currentGap.size_ == 2) {
          comp_.twoBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size_ == 1) {
          comp_.oneBaseIndel_ +=
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore_ =
              currentGap.size_ /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
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
    for (gapIter = alignmentGaps_.begin(); gapIter != alignmentGaps_.end();
         ++gapIter) {
      if (gapIter->second.size_ >= 3) {
        std::string search;
        std::vector<tandemRepeat> gapTand =
            findTandemRepeatsInSequence(gapIter->second.gapedSequence_);
        if (gapTand.empty()) {
          search = gapIter->second.gapedSequence_;
        } else {
          search = gapTand[0].repeat;
        }
        bool gapWithinTandem = false;
        if (alignObjectA_.seqBase_.seq_[gapIter->second.startPos_] == '-') {
          tandemRepeat secondTandems = findTandemRepeatOfStrInSequence(
              alignObjectB_.seqBase_.seq_, search);
          if ((int)gapIter->second.startPos_ >= secondTandems.startPos &&
              (int)gapIter->second.startPos_ + (int)gapIter->second.size_ - 1 <=
                  secondTandems.stopPos) {
            gapWithinTandem = true;
          }
          if (gapWithinTandem) {
            check = true;
          }
        } else if (alignObjectB_.seqBase_.seq_[gapIter->second.startPos_] ==
                   '-') {
          tandemRepeat secondTandems = findTandemRepeatOfStrInSequence(
              alignObjectA_.seqBase_.seq_, search);
          if ((int)gapIter->second.startPos_ >= secondTandems.startPos &&
              (int)gapIter->second.startPos_ + (int)gapIter->second.size_ - 1 <=
                  secondTandems.stopPos) {
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

std::vector<tandemRepeat> aligner::findTandemRepeatsInSequence(
    const std::string& str, int match, int mismatch, int gap,
    int minimumAlignScore) {
  uint32_t sizeChecker = 2;
  std::vector<tandemRepeat> repeats;
  std::vector<tandemRepeat>::iterator repIter;
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
          tandemRepeat tempRep =
              findTandemRepeatOfStrInSequence(tandem, repIter->repeat);
          if (tempRep.numberOfRepeats != 0) {
            alreadySmallerRepeat = true;
            break;
          }
        }
        if (!alreadySmallerRepeat) {
          repeats.push_back(tandemRepeat(tandem, numberOfRepeats, alignScore,
                                         startPos, stopPos));
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

tandemRepeat aligner::findTandemRepeatOfStrInSequence(std::string str,
                                                      std::string tandem,
                                                      int match, int mismatch,
                                                      int gap,
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
    alignScore = (int)tandem.size() * (numberOfRepeats) * match;
    if (alignScore >= minimumAlignScore) {
      foundTandem = true;
    } else {
      foundTandem = false;
    }
    if (foundTandem) {
      stopPos = pos + (int)sizeChecker - 1;
      keepSearching = false;
    }
    ++pos;
  }
  if (!foundTandem) {
    return (tandemRepeat("", numberOfRepeats, alignScore, 0, 0));
  } else {
    return (
        tandemRepeat(tandem, numberOfRepeats, alignScore, startPos, stopPos));
  }
}

tandemRepeat aligner::findTandemRepeatOfStrInSequenceDegen(
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
    alignScore = (int)tandem.size() * (numberOfRepeats) * match;
    if (alignScore >= minimumAlignScore) {
      foundTandem = true;
    } else {
      foundTandem = false;
    }
    if (foundTandem) {
      stopPos = pos + (int)sizeChecker - 1;
      keepSearching = false;
    }
    ++pos;
  }
  if (!foundTandem) {
    return (tandemRepeat("", numberOfRepeats, alignScore, 0, 0));
  } else {
    return (
        tandemRepeat(tandem, numberOfRepeats, alignScore, startPos, stopPos));
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
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); ++i) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
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
        parts_.score_ -= parts_.gapScores_.gapRightOpen_;
        parts_.score_ -= parts_.gapScores_.gapRightExtend_ * (newGap.size_ - 1);
      } else if (newGap.startPos_ == 0) {
      	parts_.score_ -= parts_.gapScores_.gapLeftOpen_;
      	parts_.score_ -= parts_.gapScores_.gapLeftExtend_ * (newGap.size_ - 1);
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
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
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
        parts_.score_ -= parts_.gapScores_.gapRightOpen_;
        parts_.score_ -= parts_.gapScores_.gapRightExtend_ * (newGap.size_ - 1);
      } else if (newGap.startPos_ == 0) {
      	parts_.score_ -= parts_.gapScores_.gapLeftOpen_;
      	parts_.score_ -= parts_.gapScores_.gapLeftExtend_ * (newGap.size_ - 1);
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
}

void aligner::noAlignSetAndScore(const baseReadObject& objectA,
		const baseReadObject& objectB) {
	noAlignSetAndScore(objectA.seqBase_, objectB.seqBase_);
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


}  // namespace bib
