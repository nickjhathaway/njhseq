//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "aligner.hpp"

namespace bibseq {


void aligner::resetCounts() {
  errors_.resetCounts();
  highQualityMatch_ = 0;
  lowQualityMatch_ = 0;
  numberOfSmallGaps_ = 0;
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
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      newGap.ref = true;
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
        ++firstOffset;
      }
      if ((((newGap.startPos + newGap.size) >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
          !countEndGaps_) {
        sumOfGappedEndsInA += newGap.size;
      } else {
        handleGapCountingInA(newGap, weighHomopolymers);
        sumOfRegularGappedInA += newGap.size;
        alignmentGaps_.insert(std::make_pair(newGap.startPos, newGap));
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      newGap.ref = false;
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
        ++secondOffset;
      }
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {
        sumOfGappedEndsInB += newGap.size;
      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
        sumOfRegularGappedInB += newGap.size;
        alignmentGaps_.insert(std::make_pair(newGap.startPos, newGap));
      }
      continue;
    }
    // alignObjectA_.seqBase_.seq_[i]!=alignObjectB_.seqBase_.seq_[i]
    if (0 > parts_.scoring_.mat_[alignObjectA_.seqBase_.seq_[i]]
                         [alignObjectB_.seqBase_.seq_[i]]) {
      ++errors_.hqMismatches_;
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
      ++highQualityMatch_;
    }
  }

  double coverageArea = len(alignObjectA_.seqBase_.seq_) - sumOfGappedEndsInA - sumOfGappedEndsInB;
  //std::cout << "cov area"  << coverageArea << std::endl;
  //std::cout << "sumOfGappedEndsInA area"  << sumOfGappedEndsInA << std::endl;
  //std::cout << "sumOfGappedEndsInB area"  << sumOfGappedEndsInB << std::endl;
  distances_.refCoverage_ = (coverageArea - sumOfRegularGappedInA)/objectA.seq_.size();
  distances_.queryCoverage_ = (coverageArea - sumOfRegularGappedInB)/objectB.seq_.size();
  distances_.percentIdentity_ = (highQualityMatch_ + lowQualityMatch_) / coverageArea;
  distances_.percentMismatch_ = (errors_.hqMismatches_ + errors_.lowKmerMismatches_ + errors_.lqMismatches_)/coverageArea;
  distances_.percentageGaps_ = (sumOfRegularGappedInA + sumOfRegularGappedInB)/ coverageArea;
  distances_.identity_ = distances_.percentIdentity_;
  distances_.ownDistance_ = errors_.hqMismatches_;
  distances_.ownGapDistance_ = 0;
  for (const auto& g : alignmentGaps_) {
  	distances_.ownGapDistance_ += g.second.homoploymerScore * g.second.size;
  	distances_.ownDistance_ += g.second.homoploymerScore * g.second.size;
  }
  distances_.eventBasedIdentity_ = 1.0 - ((alignmentGaps_.size() + errors_.hqMismatches_)/
  		static_cast<double>(errors_.hqMismatches_ + highQualityMatch_ + errors_.hqMismatches_));
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
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.size(); i++) {
    if ((i - offSet) == seqAPos) {
      return i;
    }
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++offSet;
    }
  }
  return len(alignObjectA_.seqBase_.seq_);
}
uint32_t aligner::getAlignPosForSeqBPos(uint32_t seqBPos) {
  uint32_t offSet = 0;
  for (uint32_t i = 0; i < alignObjectB_.seqBase_.seq_.size(); i++) {
    if ((i - offSet) == seqBPos) {
      return i;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++offSet;
    }
  }
  return len(alignObjectB_.seqBase_.seq_);
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
    for (uint32_t i = 0; i < start; i++) {
      // std::cout << i <<"/" << alignObjectA_.seqBase_.seq_.length() <<
      // std::endl;
      if (alignObjectA_.seqBase_.seq_[i] == '-') {
        ++firstOffset;
        gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                         alignObjectB_.seqBase_.qual_[i]);
        newGap.ref = true;
        while (alignObjectA_.seqBase_.seq_[i + 1] == '-' && (i + 1) != start) {
          newGap.gapedSequence.append(
              alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
          newGap.size++;
          newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
          i++;
          ++firstOffset;
        }
        // firstOffset += newGap.size;
        if (((newGap.startPos + newGap.size >=
              alignObjectA_.seqBase_.seq_.length()) ||
             newGap.startPos == 0) &&
            !countEndGaps_) {
          sumOfGappedEndsInA += newGap.size;
        } else {
          sumOfRegularGappedInA += newGap.size;
        }
        continue;
      }
      if (alignObjectB_.seqBase_.seq_[i] == '-') {
        ++secondOffset;
        gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                         alignObjectA_.seqBase_.qual_[i]);
        newGap.ref = false;
        while (alignObjectB_.seqBase_.seq_[i + 1] == '-' && (i + 1) != start) {
          newGap.gapedSequence.append(
              alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
          newGap.size++;
          newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
          i++;
          ++secondOffset;
        }
        if (((newGap.startPos + newGap.size) >=
                 alignObjectB_.seqBase_.seq_.length() ||
             newGap.startPos == 0) &&
            !countEndGaps_) {

        } else {
        }
        continue;
      }
    }
  }
  if (stop == 0) {
    stop = len(alignObjectA_.seqBase_.seq_);
  }
  for (uint32_t i = start; i < stop; i++) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      newGap.ref = true;
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
        ++firstOffset;
      }
      // firstOffset += newGap.size;
      if (((newGap.startPos + newGap.size >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
          !countEndGaps_) {
        sumOfGappedEndsInA += newGap.size;
      } else {
        sumOfRegularGappedInA += newGap.size;
        handleGapCountingInA(newGap, weighHomopolyer);
        alignmentGaps_.insert(std::make_pair(newGap.startPos, newGap));
        //editDistance_ += newGap.size;
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      newGap.ref = false;
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
        ++secondOffset;
      }
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolyer);
        alignmentGaps_.insert(std::make_pair(newGap.startPos, newGap));
        //editDistance_ += newGap.size;
      }
      continue;
    }

    if ( 0 > parts_.scoring_.mat_[alignObjectA_.seqBase_.seq_[i]]
                                  [alignObjectB_.seqBase_.seq_[i]]) {
      //++editDistance_;
      std::string firstKmer = "";
      std::string secondKmer = "";
      uint32_t firstPos = i;
      uint32_t secondPos = i;
      if ((int)i - firstOffset - kLength / 2 < 0) {
        firstPos = 0;
        firstKmer = objectA.seq_.substr(0, kLength);
      } else if (i - firstOffset + kLength / 2 >= objectA.seq_.size()) {
        firstKmer =
            objectA.seq_.substr(objectA.seq_.size() - kLength, kLength);
        firstPos = (int)objectA.seq_.size()  - kLength;
      } else {
        firstKmer = objectA.seq_.substr(i - firstOffset - kLength / 2, kLength);
        firstPos = i - firstOffset - kLength / 2;
      }

      if (i - secondOffset + kLength / 2 >= objectB.seq_.size()) {
        secondKmer =
            objectB.seq_.substr(objectB.seq_.size()  - kLength, kLength);
        secondPos = (int)objectB.seq_.size() - kLength;
      } else if ((int)i - secondOffset - kLength / 2 < 0) {
        secondKmer = objectB.seq_.substr(0, kLength);
        secondPos = 0;
      } else {
        secondKmer =
            objectB.seq_.substr(i - secondOffset - kLength / 2, kLength);
        secondPos = i - secondOffset - kLength / 2;
      }
      if (checkKmer &&
          (kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                     kMaps_.runCutOff_) ||
           kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                     kMaps_.runCutOff_))) {
        ++errors_.lowKmerMismatches_;
        lowKmerMismatches_.insert(std::make_pair(
            i,
            mismatch(
                alignObjectA_.seqBase_.seq_[i], alignObjectA_.seqBase_.qual_[i],
                objectA.getLeadQual(i - firstOffset),
                objectA.getTrailQual(i - firstOffset), i - firstOffset,
                alignObjectB_.seqBase_.seq_[i], alignObjectB_.seqBase_.qual_[i],
                objectB.getLeadQual(i - secondOffset),
                objectB.getTrailQual(i - secondOffset), i - secondOffset,
                kMaps_.kmersByPos_[secondPos][secondKmer].readCnt_,
                kMaps_.kmersAnyWhere_[secondKmer].readCnt_)));
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
                kMaps_.kmersByPos_[secondPos][secondKmer].readCnt_,
                kMaps_.kmersAnyWhere_[secondKmer].readCnt_)));
        if (usingQuality) {
          if (objectA.checkQual(i - firstOffset, primaryQual_, secondaryQual_,
                                qualThresWindow_) &&
              objectB.checkQual(i - secondOffset, primaryQual_, secondaryQual_,
                                qualThresWindow_)) {
            ++errors_.hqMismatches_;
          } else {
            ++errors_.lqMismatches_;
          }
        } else {
          ++errors_.hqMismatches_;
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
            highQualityMatch_++;
          } else {
            lowQualityMatch_++;
          }
        } else {
          highQualityMatch_++;
        }
      } else {
        ++highQualityMatch_;
      }
    }
  }
  std::string noGapsAlignB =
      seqUtil::removeGapsReturn(alignObjectB_.seqBase_.seq_);
  distances_.percentIdentity_ = (double)(highQualityMatch_ + lowQualityMatch_) /
                     ((double)noGapsAlignB.size() - sumOfGappedEndsInA -
                      sumOfRegularGappedInA);
  distances_.queryCoverage_ =
      (double)(noGapsAlignB.size() - sumOfGappedEndsInA) / objectB.seq_.size();
  distances_.percentageGaps_ =
      (double)(alignObjectB_.seqBase_.seq_.size() - noGapsAlignB.size()) /
      alignObjectB_.seqBase_.seq_.size();
  distances_.eventBasedIdentity_ = 1.0 - ((errors_.hqMismatches_ + errors_.lqMismatches_ + errors_.lowKmerMismatches_ + alignmentGaps_.size())/
  		static_cast<double>(highQualityMatch_ + lowQualityMatch_ + errors_.hqMismatches_ + errors_.lqMismatches_ + errors_.lowKmerMismatches_));
}
// no kmer checking
bool aligner::checkAlignmentBool(const baseReadObject& objectA,
                                 const baseReadObject& objectB,
                                 const runningParameters& runParams,
                                 bool usingQuality, bool doingMatchQuality,
                                 bool weighHomopolymer) {
	resetAlignmentInfo();;
  int firstOffset = 0;
  int secondOffset = 0;
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); ++i) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        ++newGap.size;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        ++i;
      }
      if (((newGap.startPos + newGap.size >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInA(newGap, weighHomopolymer);
        if (errors_.largeBaseIndel_ > runParams.errors_.largeBaseIndel_) {
          return false;
        }
        if (errors_.twoBaseIndel_ > runParams.errors_.twoBaseIndel_) {
          return false;
        }
        if (errors_.oneBaseIndel_ > runParams.errors_.oneBaseIndel_) {
          return false;
        }

        // firstAlignmentGaps.insert(std::make_pair(newGap.startPos, newGap));
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      secondOffset++;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        ++i;
      }
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolymer);
        if (errors_.largeBaseIndel_ > runParams.errors_.largeBaseIndel_) {
          return false;
        }
        if (errors_.twoBaseIndel_ > runParams.errors_.twoBaseIndel_) {
          return false;
        }
        if (errors_.oneBaseIndel_ > runParams.errors_.oneBaseIndel_) {
          return false;
        }
        // secondAlignmentGaps.insert(std::make_pair(newGap.startPos, newGap));
      }
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      if (usingQuality) {
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_)) {
          ++errors_.hqMismatches_;
          if (errors_.hqMismatches_ > runParams.errors_.hqMismatches_) {
            return false;
          }
        } else {
          ++errors_.lqMismatches_;
          if (errors_.lqMismatches_ > runParams.errors_.lqMismatches_) {
            return false;
          }
        }
      } else {
        ++errors_.hqMismatches_;
        if (errors_.hqMismatches_ > runParams.errors_.hqMismatches_) {
          return false;
        }
      }
    }
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      if (usingQuality) {
        if (doingMatchQuality) {
          if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                         secondaryQual_, qualThresWindow_) &&
              objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                         secondaryQual_, qualThresWindow_)) {
            highQualityMatch_++;
          } else {
            lowQualityMatch_++;
          }
        } else {
          highQualityMatch_++;
        }
      } else {
        ++highQualityMatch_;
      }
    }
  }
  return true;
}
errorProfile aligner::checkAlignment(const baseReadObject& objectA,
                                     const baseReadObject& objectB,
                                     const runningParameters& runParams,
                                     bool checkKmers, bool kmersByPosition,
                                     bool weighHomopolymers) {
	resetAlignmentInfo();;
  int firstOffset = 0;
  int secondOffset = 0;
  /*if (score<0) {
   return false;
   }*/
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      firstOffset++;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInA(newGap, weighHomopolymers);
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {

      secondOffset++;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
      }
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      std::string firstKmer = "";
      std::string secondKmer = "";
      uint32_t firstPos = i;
      uint32_t secondPos = i;
      if ((int)i - firstOffset - kMaps_.kLength_ / 2 < 0) {
        firstPos = 0;
        firstKmer = objectA.seqBase_.seq_.substr(0, kMaps_.kLength_);
      } else if (i - firstOffset + kMaps_.kLength_ / 2 >=
                 objectA.seqBase_.seq_.size()) {
        firstKmer = objectA.seqBase_.seq_.substr(
            objectA.seqBase_.seq_.size() - 1 - kMaps_.kLength_,
            kMaps_.kLength_);
        firstPos = (int)objectA.seqBase_.seq_.size() - 1 - kMaps_.kLength_;
      } else {
        firstKmer = objectA.seqBase_.seq_.substr(
            i - firstOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        firstPos = i - firstOffset - kMaps_.kLength_ / 2;
      }
      if (i - secondOffset + kMaps_.kLength_ / 2 >=
          objectB.seqBase_.seq_.size()) {
        secondKmer = objectB.seqBase_.seq_.substr(
            objectB.seqBase_.seq_.size() - 1 - kMaps_.kLength_,
            kMaps_.kLength_);
        secondPos = (int)objectB.seqBase_.seq_.size() - 1 - kMaps_.kLength_;
      } else if ((int)i - secondOffset - kMaps_.kLength_ / 2 < 0) {
        secondKmer = objectB.seqBase_.seq_.substr(0, kMaps_.kLength_);
        secondPos = 0;
      } else {
        secondKmer = objectB.seqBase_.seq_.substr(
            i - secondOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        secondPos = i - secondOffset - kMaps_.kLength_ / 2;
      }
      if ((kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                     kMaps_.runCutOff_) ||
           kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                     kMaps_.runCutOff_)) &&
          checkKmers) {
        ++errors_.lowKmerMismatches_;
      } else {
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_)) {
          errors_.hqMismatches_++;
        } else {
          errors_.lqMismatches_++;
        }
      }
    }
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      highQualityMatch_++;
    }
  }
  return errors_;
}
errorProfile aligner::checkAlignmentLowKmerQual(
    const baseReadObject& objectA, const baseReadObject& objectB,
    const runningParameters& runParams, bool checkKmers, bool kmersByPosition,
    bool weighHomopolymers) {
	resetAlignmentInfo();;
  int firstOffset = 0;
  int secondOffset = 0;
  /*if (score<0) {
   return false;
   }*/
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      firstOffset++;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInA(newGap, weighHomopolymers);
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {

      secondOffset++;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
      }
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      std::string firstKmer = "";
      std::string secondKmer = "";
      uint32_t firstPos = i;
      uint32_t secondPos = i;
      if ((int)i - firstOffset - kMaps_.kLength_ / 2 < 0) {
        firstPos = 0;
        firstKmer = objectA.seqBase_.seq_.substr(0, kMaps_.kLength_);
      } else if (i - firstOffset + kMaps_.kLength_ / 2 >=
                 objectA.seqBase_.seq_.size()) {
        firstKmer = objectA.seqBase_.seq_.substr(
            objectA.seqBase_.seq_.size() - 1 - kMaps_.kLength_,
            kMaps_.kLength_);
        firstPos = (int)objectA.seqBase_.seq_.size() - 1 - kMaps_.kLength_;
      } else {
        firstKmer = objectA.seqBase_.seq_.substr(
            i - firstOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        firstPos = i - firstOffset - kMaps_.kLength_ / 2;
      }
      if (i - secondOffset + kMaps_.kLength_ / 2 >=
          objectB.seqBase_.seq_.size()) {
        secondKmer = objectB.seqBase_.seq_.substr(
            objectB.seqBase_.seq_.size() - 1 - kMaps_.kLength_,
            kMaps_.kLength_);
        secondPos = (int)objectB.seqBase_.seq_.size() - 1 - kMaps_.kLength_;
      } else if ((int)i - secondOffset - kMaps_.kLength_ / 2 < 0) {
        secondKmer = objectB.seqBase_.seq_.substr(0, kMaps_.kLength_);
        secondPos = 0;
      } else {
        secondKmer = objectB.seqBase_.seq_.substr(
            i - secondOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        secondPos = i - secondOffset - kMaps_.kLength_ / 2;
      }
      if ((kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                     kMaps_.runCutOff_) ||
           kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                     kMaps_.runCutOff_)) &&
          checkKmers) {
        ++errors_.lowKmerMismatches_;
        // std::cout<<"primaryQualLowKmer_: "<<primaryQualLowKmer_<<std::endl;
        // std::cout<<"secondaryQualLowKmer_"<<secondaryQualLowKmer_<<std::endl;
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQualLowKmer_,
                                       secondaryQualLowKmer_,
                                       qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQualLowKmer_,
                                       secondaryQualLowKmer_,
                                       qualThresWindow_)) {
          errors_.hqMismatches_++;
        } else {
          errors_.lqMismatches_++;
        }
      } else {
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_)) {
          errors_.hqMismatches_++;
        } else {
          errors_.lqMismatches_++;
        }
      }
    }
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      highQualityMatch_++;
    }
  }
  return errors_;
}

errorProfile aligner::checkAlignmentBothQualKmer(
    const baseReadObject& objectA, const baseReadObject& objectB,
    const runningParameters& runParams, bool checkKmers, bool kmersByPosition,
    bool weighHomopolymers) {
	resetAlignmentInfo();;
  int firstOffset = 0;
  int secondOffset = 0;
  /*if (score<0) {
   return false;
   }*/
  // std::cout <<"cab1" << std::endl;
  // std::cout << alignObjectA_.seqBase_.seq_.size() << ":"
  //<< alignObjectA_.seqBase_.seq_ << std::endl;
  // std::cout << alignObjectB_.seqBase_.seq_.size() << ":"
  //<< alignObjectB_.seqBase_.seq_ << std::endl;

  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
    // std::cout <<"cab2.1" << std::endl;
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      firstOffset++;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
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
                       alignObjectA_.seqBase_.qual_[i]);
      // std::cout <<"cab3.2" << std::endl;
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
      }
      // std::cout <<"cab3.3" << std::endl;
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
      }
      // std::cout <<"cab3.4" << std::endl;
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      // std::cout <<"cab4.1" << std::endl;
      std::string firstKmer = "";
      std::string secondKmer = "";
      uint32_t firstPos = i;
      uint32_t secondPos = i;
      if ((int)i - firstOffset - kMaps_.kLength_ / 2 < 0) {
        firstPos = 0;
        firstKmer = objectA.seqBase_.seq_.substr(0, kMaps_.kLength_);
      } else if (i - firstOffset + kMaps_.kLength_ / 2 >=
                 objectA.seqBase_.seq_.size()) {
        firstKmer = objectA.seqBase_.seq_.substr(
            objectA.seqBase_.seq_.size() + 1 - kMaps_.kLength_,
            kMaps_.kLength_);
        firstPos = (int)objectA.seqBase_.seq_.size() + 1 - kMaps_.kLength_;
      } else {
        firstKmer = objectA.seqBase_.seq_.substr(
            i - firstOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        firstPos = i - firstOffset - kMaps_.kLength_ / 2;
      }
      // std::cout <<"cab4.2" << std::endl;
      if (i - secondOffset + kMaps_.kLength_ / 2 >=
          objectB.seqBase_.seq_.size()) {
        secondKmer = objectB.seqBase_.seq_.substr(
            objectB.seqBase_.seq_.size() + 1 - kMaps_.kLength_,
            kMaps_.kLength_);
        secondPos = (int)objectB.seqBase_.seq_.size() + 1 - kMaps_.kLength_;
      } else if ((int)i - secondOffset - kMaps_.kLength_ / 2 < 0) {
        secondKmer = objectB.seqBase_.seq_.substr(0, kMaps_.kLength_);
        secondPos = 0;
      } else {
        secondKmer = objectB.seqBase_.seq_.substr(
            i - secondOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        secondPos = i - secondOffset - kMaps_.kLength_ / 2;
      }
      // std::cout <<"cab4.3" << std::endl;
      /*if ((kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                     kMaps_.qualRunCutOff_) ||
           kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                     kMaps_.qualRunCutOff_)) &&
          checkKmers)*/
      if (checkKmers &&
      		(kMaps_.isKmerLowFreqByQual(firstKmer, firstPos, kmersByPosition,
                                      alignObjectA_.seqBase_.qual_[i]) ||
           kMaps_.isKmerLowFreqByQual(secondKmer, secondPos, kmersByPosition,
                                      alignObjectA_.seqBase_.qual_[i])) ) {
        ++errors_.lowKmerMismatches_;
        /*
       //++errors_.lowKmerMismatches_;
       // std::cout<<"primaryQualLowKmer_: "<<primaryQualLowKmer_<<std::endl;
       // std::cout<<"secondaryQualLowKmer_"<<secondaryQualLowKmer_<<std::endl;
       if (objectA.checkQual(i - firstOffset, primaryQualLowKmer_,
                             secondaryQualLowKmer_, qualThresWindow_) &&
           objectB.checkQual(i - secondOffset, primaryQualLowKmer_,
                             secondaryQualLowKmer_, qualThresWindow_)) {
         if ((kMaps_.isKmerLowFreqByQual(firstKmer, firstPos, kmersByPosition,
       alignObjectA_.seqBase_.qual_[i] ) ||
              kMaps_.isKmerLowFreqByQual(secondKmer, secondPos, kmersByPosition,
                        alignObjectA_.seqBase_.qual_[i])) &&
             checkKmers) {
         /STARif ((  kMaps_.isKmerLowFrequency(firstKmer, firstPos,
       kmersByPosition,
                                        kMaps_.runCutOff_) ||
              kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                        kMaps_.runCutOff_)) &&
             checkKmers) {STAR/
           ++errors_.lowKmerMismatches_;
         } else {
           errors_.hqMismatches_++;
         }
       } else {
         errors_.lqMismatches_++;
       }*/
        // std::cout <<"cab4.4" << std::endl;
      } else {
        // std::cout <<"cab4.5" << std::endl;
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_)) {
          errors_.hqMismatches_++;
        } else {
          errors_.lqMismatches_++;
        }
      }
    }
    // std::cout <<"cab5.1" << std::endl;
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      highQualityMatch_++;
    }
    // std::cout <<"cab6.1" << std::endl;
  }
  // std::cout <<"cab7.1" << std::endl;
  return errors_;
}

errorProfile aligner::checkAlignmentBothRegKmer(
    const baseReadObject& objectA, const baseReadObject& objectB,
    const runningParameters& runParams, bool checkKmers, bool kmersByPosition,
    bool weighHomopolymers) {
	//std::cout << "begin" << std::endl;
	//alignObjectA_.seqBase_.outPutSeq(std::cout);
	//alignObjectB_.seqBase_.outPutSeq(std::cout);
	//std::cout << alignObjectA_.seqBase_.seq_.size() << std::endl;;
	//std::cout << alignObjectB_.seqBase_.seq_.size() << std::endl;

	resetAlignmentInfo();;
  int firstOffset = 0;
  int secondOffset = 0;
  /*if (score<0) {
   return false;
   }*/
  // std::cout <<"cab1" << std::endl;
  // std::cout << alignObjectA_.seqBase_.seq_.size() << ":"
  //<< alignObjectA_.seqBase_.seq_ << std::endl;
  // std::cout << alignObjectB_.seqBase_.seq_.size() << ":"
  //<< alignObjectB_.seqBase_.seq_ << std::endl;

  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
  	//std::cout << i << std::endl;
    // std::cout <<"cab2.1" << std::endl;
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      firstOffset++;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
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
                       alignObjectA_.seqBase_.qual_[i]);
      // std::cout <<"cab3.2" << std::endl;
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
      }
      // std::cout <<"cab3.3" << std::endl;
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
      }
      // std::cout <<"cab3.4" << std::endl;
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      //std::cout <<"cab4.1" << std::endl;
      std::string firstKmer = "";
      std::string secondKmer = "";
      uint32_t firstPos = i;
      uint32_t secondPos = i;

      if ((int)i - firstOffset - kMaps_.kLength_ / 2 < 0) {
        firstPos = 0;
        firstKmer = objectA.seqBase_.seq_.substr(0, kMaps_.kLength_);
      } else if (i - firstOffset + kMaps_.kLength_ / 2 >=
                 objectA.seqBase_.seq_.size()) {
        firstKmer = objectA.seqBase_.seq_.substr(
            objectA.seqBase_.seq_.size() - kMaps_.kLength_,
            kMaps_.kLength_);
        firstPos = (int)objectA.seqBase_.seq_.size() - kMaps_.kLength_;
      } else {
        firstKmer = objectA.seqBase_.seq_.substr(
            i - firstOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        firstPos = i - firstOffset - kMaps_.kLength_ / 2;
      }
      //std::cout <<"cab4.2" << std::endl;
      if (i - secondOffset + kMaps_.kLength_ / 2 >=
          objectB.seqBase_.seq_.size()) {
        secondKmer = objectB.seqBase_.seq_.substr(
            objectB.seqBase_.seq_.size()  - kMaps_.kLength_,
            kMaps_.kLength_);
        secondPos = (int)objectB.seqBase_.seq_.size() - kMaps_.kLength_;
      } else if ((int)i - secondOffset - kMaps_.kLength_ / 2 < 0) {
        secondKmer = objectB.seqBase_.seq_.substr(0, kMaps_.kLength_);
        secondPos = 0;
      } else {
        secondKmer = objectB.seqBase_.seq_.substr(
            i - secondOffset - kMaps_.kLength_ / 2, kMaps_.kLength_);
        secondPos = i - secondOffset - kMaps_.kLength_ / 2;
      }

      if (checkKmers &&
      		(kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                     kMaps_.runCutOff_) ||
           kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                     kMaps_.runCutOff_)) ) {
/*      	std::cout << "firstKmer: " << firstKmer << std::endl;
      	std::cout << "freq" << kMaps_.kmersByPos_[firstPos][firstKmer].readCnt_  << " " << kMaps_.runCutOff_ << std::endl;
      	std::cout << "secondKmer: " << secondKmer << std::endl;
      	std::cout << "freq" << kMaps_.kmersByPos_[secondPos][secondKmer].readCnt_ << " " << kMaps_.runCutOff_ << std::endl << std::endl;
*/
        ++errors_.lowKmerMismatches_;
      } else {
        // std::cout <<"cab4.5" << std::endl;
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_)) {
          errors_.hqMismatches_++;
        } else {
          errors_.lqMismatches_++;
        }
      }
      /*
      if ((kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                     kMaps_.qualRunCutOff_) ||
           kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                     kMaps_.qualRunCutOff_)) &&
          checkKmers){
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQualLowKmer_,
                                       secondaryQualLowKmer_,
                                       qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQualLowKmer_,
                                       secondaryQualLowKmer_,
                                       qualThresWindow_)) {
          if ((kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                         kMaps_.runCutOff_) ||
               kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                         kMaps_.runCutOff_)) &&
              checkKmers) {
            ++errors_.lowKmerMismatches_;
          } else {
            errors_.hqMismatches_++;
          }
        } else {
          errors_.lqMismatches_++;
        }
        // std::cout <<"cab4.4" << std::endl;
      } else {
        // std::cout <<"cab4.5" << std::endl;
        if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_) &&
            objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                       secondaryQual_, qualThresWindow_)) {
          errors_.hqMismatches_++;
        } else {
          errors_.lqMismatches_++;
        }
      }*/
      //std::cout <<"cab4.end" << std::endl;
    }
    // std::cout <<"cab5.1" << std::endl;
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      highQualityMatch_++;
    }
    // std::cout <<"cab6.1" << std::endl;
  }
  // std::cout <<"cab7.1" << std::endl;
  //std::cout << "end" << std::endl;
  return errors_;
}

bool aligner::checkAlignmentBool(const baseReadObject& objectA,
                                 const baseReadObject& objectB,
                                 const runningParameters& runParams,
                                 int kLength, bool kmersByPosition,
                                 bool usingQuality, bool doingMatchQuality,
                                 bool weighHomopolymers) {
	resetAlignmentInfo();;
  int firstOffset = 0;
  int secondOffset = 0;
  /*if (score<0) {
      return false;
  }*/
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      firstOffset++;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size >=
            alignObjectA_.seqBase_.seq_.length()) ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInA(newGap, weighHomopolymers);
        if (errors_.largeBaseIndel_ > runParams.errors_.largeBaseIndel_) {
          return false;
        }
        if (errors_.twoBaseIndel_ > runParams.errors_.twoBaseIndel_) {
          return false;
        }
        if (errors_.oneBaseIndel_ > runParams.errors_.oneBaseIndel_) {
          return false;
        }
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {

      secondOffset++;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
      }
      if (((newGap.startPos + newGap.size) >=
               alignObjectB_.seqBase_.seq_.length() ||
           newGap.startPos == 0) &&
          !countEndGaps_) {

      } else {
        handleGapCountingInB(newGap, weighHomopolymers);
        if (errors_.largeBaseIndel_ > runParams.errors_.largeBaseIndel_) {
          return false;
        }
        if (errors_.twoBaseIndel_ > runParams.errors_.twoBaseIndel_) {
          return false;
        }
        if (errors_.oneBaseIndel_ > runParams.errors_.oneBaseIndel_) {
          return false;
        }
      }
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      std::string firstKmer = "";
      std::string secondKmer = "";
      uint32_t firstPos = i;
      uint32_t secondPos = i;
      if ((int)i - firstOffset - kLength / 2 < 0) {
        firstPos = 0;
        firstKmer = objectA.seqBase_.seq_.substr(0, kLength);
      } else if (i - firstOffset + kLength / 2 >=
                 objectA.seqBase_.seq_.size()) {
        firstKmer = objectA.seqBase_.seq_.substr(
            objectA.seqBase_.seq_.size() - kLength, kLength);
        firstPos = (int)objectA.seqBase_.seq_.size() - kLength;
      } else {
        firstKmer = objectA.seqBase_.seq_.substr(i - firstOffset - kLength / 2,
                                                 kLength);
        firstPos = i - firstOffset - kLength / 2;
      }
      if (i - secondOffset + kLength / 2 >= objectB.seqBase_.seq_.size()) {
        secondKmer = objectB.seqBase_.seq_.substr(
            objectB.seqBase_.seq_.size() - kLength, kLength);
        secondPos = (int)objectB.seqBase_.seq_.size() - kLength;
      } else if ((int)i - secondOffset - kLength / 2 < 0) {
        secondKmer = objectB.seqBase_.seq_.substr(0, kLength);
        secondPos = 0;
      } else {
        secondKmer = objectB.seqBase_.seq_.substr(
            i - secondOffset - kLength / 2, kLength);
        secondPos = i - secondOffset - kLength / 2;
      }

      if (kMaps_.isKmerLowFrequency(firstKmer, firstPos, kmersByPosition,
                                    kMaps_.runCutOff_) ||
          kMaps_.isKmerLowFrequency(secondKmer, secondPos, kmersByPosition,
                                    kMaps_.runCutOff_)) {
        ++errors_.lowKmerMismatches_;
      } else {
        if (usingQuality) {
          if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                         secondaryQual_, qualThresWindow_) &&
              objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                         secondaryQual_, qualThresWindow_)) {
            errors_.hqMismatches_++;
            if (errors_.hqMismatches_ > runParams.errors_.hqMismatches_) {
              return false;
            }
          } else {
            errors_.lqMismatches_++;
            if (errors_.lqMismatches_ > runParams.errors_.lqMismatches_) {
              return false;
            }
          }
        } else {
          errors_.hqMismatches_++;
          if (errors_.hqMismatches_ > runParams.errors_.hqMismatches_) {
            return false;
          }
        }
      }
    }
    if (alignObjectA_.seqBase_.seq_[i] == alignObjectB_.seqBase_.seq_[i]) {
      if (usingQuality) {
        if (doingMatchQuality) {
          if (objectA.seqBase_.checkQual(i - firstOffset, primaryQual_,
                                         secondaryQual_, qualThresWindow_) &&
              objectB.seqBase_.checkQual(i - secondOffset, primaryQual_,
                                         secondaryQual_, qualThresWindow_)) {
            highQualityMatch_++;
          } else {
            lowQualityMatch_++;
          }
        } else {
          highQualityMatch_++;
        }
      } else {
        highQualityMatch_++;
      }
    }
  }
  return true;
}

void aligner::simpleProfileAlignment(const std::string& firstSeq,
                                     const std::string& secondSeq,
                                     int& forwardOverHangSizeA,
                                     int& backOverHangSizeA,
                                     int& forwardOverHangSizeB,
                                     int& backOverHangSizeB) {
  resetAlignmentInfo();
  int firstOffset = 0;
  int secondOffset = 0;
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      ++firstOffset;
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      newGap.ref = true;
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
      }
      firstOffset += newGap.size;
      if (newGap.startPos + newGap.size >=
          alignObjectA_.seqBase_.seq_.length()) {
        backOverHangSizeA = newGap.size;
      } else if (newGap.startPos == 0) {
        forwardOverHangSizeA = newGap.size;
      } else {
        alignmentGaps_.insert({newGap.startPos, newGap});
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      ++secondOffset;
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      newGap.ref = false;
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
      }
      secondOffset += newGap.size;
      if (newGap.startPos + newGap.size >=
          alignObjectB_.seqBase_.seq_.length()) {
        backOverHangSizeB = newGap.size;
      } else if (newGap.startPos == 0) {
        forwardOverHangSizeB = newGap.size;
      } else {
        alignmentGaps_.insert({newGap.startPos, newGap});
      }
      continue;
    }
    if (alignObjectA_.seqBase_.seq_[i] != alignObjectB_.seqBase_.seq_[i]) {
      ++errors_.hqMismatches_;
    } else if (alignObjectA_.seqBase_.seq_[i] ==
               alignObjectB_.seqBase_.seq_[i]) {
      ++highQualityMatch_;
    } else {
      std::cout << "well this shouldn't be happening " << std::endl;
      std::cout << "A: " << alignObjectA_.seqBase_.seq_[i]
                << " B:" << alignObjectB_.seqBase_.seq_[i] << std::endl;
    }
  }
  return;
}
void aligner::handleGapCountingInA(gap& currentGap, bool weighHomopolymers) {

  if (!seqUtil::isHomopolyer(currentGap.gapedSequence) || !weighHomopolymers) {
    if (currentGap.size >= 3) {
      ++errors_.largeBaseIndel_;
      currentGap.homoploymerScore = 1;
    } else if (currentGap.size == 2) {
      ++errors_.twoBaseIndel_;
      currentGap.homoploymerScore = 1;
    } else if (currentGap.size == 1) {
      ++errors_.oneBaseIndel_;
      currentGap.homoploymerScore = 1;
    }
  } else {
    double firstBases = 0.00;
    double secondBases = 0.00;
    // bases++;
    // forwards
    int cursor = 0;
    while ((currentGap.startPos + currentGap.size + cursor) <
               alignObjectA_.seqBase_.seq_.size() &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos + currentGap.size +
                                       cursor] == currentGap.gapedSequence[0]) {
      firstBases++;
      cursor++;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos - cursor) >= 0 &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos - cursor] ==
               currentGap.gapedSequence[0]) {
      firstBases++;
      cursor++;
    }
    cursor = 0;
    // forwards
    while ((currentGap.startPos + cursor) <
               alignObjectB_.seqBase_.seq_.size() &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos + cursor] ==
               currentGap.gapedSequence[0]) {
      secondBases++;
      cursor++;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos - cursor) >= 0 &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos - cursor] ==
               currentGap.gapedSequence[0]) {
      secondBases++;
      cursor++;
    }
    // std::cout<<"secondBases: "<<secondBases<<" firstBases:
    // "<<firstBases<<std::endl;
    // std::cout<<"secondBasesSize: "<<alignObjectB_.seqBase_.cnt_<<"
    // firstBasesSize: "<<alignObjectA_.seqBase_.cnt_<<std::endl;
    // std::cout<<"size: "<<currentGap.size<<std::endl<<std::endl;
    if (secondBases == 0 || firstBases == 0) {
      if (currentGap.size >= 3) {
        ++errors_.largeBaseIndel_;
        currentGap.homoploymerScore = 1;
      } else if (currentGap.size == 2) {
        ++errors_.twoBaseIndel_;
        currentGap.homoploymerScore = 1;
      } else if (currentGap.size == 1) {
        ++errors_.oneBaseIndel_;
        currentGap.homoploymerScore = 1;
      }
      // //std::cout<<"mark 1"<<std::endl;
    } else {
      // //std::cout<<"mark 2"<<std::endl;

      if (secondBases < firstBases) {
        // //std::cout<<"mark 3"<<std::endl;
        if (currentGap.size >= 3) {
          double currentScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++errors_.largeBaseIndel_;
            currentGap.homoploymerScore = 1;
          } else {
            errors_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore = currentScore;
          }

        } else if (currentGap.size == 2) {
          errors_.twoBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size == 1) {
          errors_.oneBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
      } else {
        // //std::cout<<"mark 4"<<std::endl;
        if (currentGap.size >= 3) {
          double currentScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++errors_.largeBaseIndel_;
            currentGap.homoploymerScore = currentScore;
          } else {
            errors_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore = currentScore;
          }

        } else if (currentGap.size == 2) {
          errors_.twoBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size == 1) {
          errors_.oneBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
      }
    }
  }
}
void aligner::handleGapCountingInB(gap& currentGap, bool weighHomopolymers) {

  if (!seqUtil::isHomopolyer(currentGap.gapedSequence) || !weighHomopolymers) {
    if (currentGap.size >= 3) {
      ++errors_.largeBaseIndel_;
      currentGap.homoploymerScore = 1;
    } else if (currentGap.size == 2) {
      ++errors_.twoBaseIndel_;
      currentGap.homoploymerScore = 1;
    } else if (currentGap.size == 1) {
      ++errors_.oneBaseIndel_;
      currentGap.homoploymerScore = 1;
    }
  } else {
    double firstBases = 0.00;
    double secondBases = 0.00;
    // bases++;
    // forwards
    int cursor = 0;
    while ((currentGap.startPos + cursor) <
               alignObjectA_.seqBase_.seq_.size() &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos + cursor] ==
               currentGap.gapedSequence[0]) {
      firstBases++;
      cursor++;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos - cursor) >= 0 &&
           alignObjectA_.seqBase_.seq_[currentGap.startPos - cursor] ==
               currentGap.gapedSequence[0]) {
      firstBases++;
      cursor++;
    }
    // forwards
    cursor = 0;
    while ((currentGap.startPos + currentGap.size + cursor) <
               alignObjectB_.seqBase_.seq_.size() &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos + currentGap.size +
                                       cursor] == currentGap.gapedSequence[0]) {
      secondBases++;
      cursor++;
    }
    cursor = 1;
    // backwards
    while (((int)currentGap.startPos - cursor) >= 0 &&
           alignObjectB_.seqBase_.seq_[currentGap.startPos - cursor] ==
               currentGap.gapedSequence[0]) {
      secondBases++;
      cursor++;
    }

    // std::cout<<"secondBases: "<<secondBases<<" firstBases:
    // "<<firstBases<<std::endl;
    // std::cout<<"secondBasesSize: "<<alignObjectB_.seqBase_.cnt_<<"
    // firstBasesSize: "<<alignObjectA_.seqBase_.cnt_<<std::endl;
    // std::cout<<"size: "<<currentGap.size<<std::endl<<std::endl;

    if (secondBases == 0 || firstBases == 0) {
      if (currentGap.size >= 3) {
        ++errors_.largeBaseIndel_;
        currentGap.homoploymerScore = 1;
      } else if (currentGap.size == 2) {
        ++errors_.twoBaseIndel_;
        currentGap.homoploymerScore = 1;
      } else if (currentGap.size == 1) {
        ++errors_.oneBaseIndel_;
        currentGap.homoploymerScore = 1;
      }
      // //std::cout<<"mark 1"<<std::endl;
    } else {
      // //std::cout<<"mark 2"<<std::endl;

      if (secondBases < firstBases) {
        // //std::cout<<"mark 3"<<std::endl;
        if (currentGap.size >= 3) {
          double currentScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++errors_.largeBaseIndel_;
            currentGap.homoploymerScore = 1;
          } else {
            errors_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore = currentScore;
          }
        } else if (currentGap.size == 2) {
          errors_.twoBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size == 1) {
          errors_.oneBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
      } else {
        // //std::cout<<"mark 4"<<std::endl;
        if (currentGap.size >= 3) {
          double currentScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          if (currentScore > 1) {
            ++errors_.largeBaseIndel_;
            currentGap.homoploymerScore = 1;
          } else {
            errors_.largeBaseIndel_ += currentScore;
            currentGap.homoploymerScore = currentScore;
          }
        } else if (currentGap.size == 2) {
          errors_.twoBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        } else if (currentGap.size == 1) {
          errors_.oneBaseIndel_ +=
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
          currentGap.homoploymerScore =
              currentGap.size /
              ((firstBases * alignObjectA_.seqBase_.cnt_ +
                secondBases * alignObjectB_.seqBase_.cnt_) /
               (alignObjectA_.seqBase_.cnt_ + alignObjectB_.seqBase_.cnt_));
        }
      }
    }
  }
}

void aligner::outPutParameterInfo(std::ostream& out) const {
  out << "numberOfOneIndel:" << errors_.oneBaseIndel_
      << " numberOfTwoIndel:" << errors_.twoBaseIndel_
      << " numberOfLargeGaps:" << errors_.largeBaseIndel_
      << " highQualityMismatch:" << errors_.hqMismatches_
      << " lowQualityMismatch:" << errors_.lqMismatches_
      << " lowKmerMismatch:" << errors_.lowKmerMismatches_ << std::endl;
}

// check for tandem repeat gaps
bool aligner::checkForTandemRepeatGap() {
  bool check = false;
  if (errors_.largeBaseIndel_ == 1) {
    std::map<uint32_t, gap>::iterator gapIter;
    for (gapIter = alignmentGaps_.begin(); gapIter != alignmentGaps_.end();
         ++gapIter) {
      if (gapIter->second.size >= 3) {
        std::string search;
        std::vector<tandemRepeat> gapTand =
            findTandemRepeatsInSequence(gapIter->second.gapedSequence);
        if (gapTand.empty()) {
          search = gapIter->second.gapedSequence;
        } else {
          search = gapTand[0].repeat;
        }
        bool gapWithinTandem = false;
        if (alignObjectA_.seqBase_.seq_[gapIter->second.startPos] == '-') {
          tandemRepeat secondTandems = findTandemRepeatOfStrInSequence(
              alignObjectB_.seqBase_.seq_, search);
          if ((int)gapIter->second.startPos >= secondTandems.startPos &&
              (int)gapIter->second.startPos + (int)gapIter->second.size - 1 <=
                  secondTandems.stopPos) {
            gapWithinTandem = true;
          }
          if (gapWithinTandem) {
            check = true;
          }
        } else if (alignObjectB_.seqBase_.seq_[gapIter->second.startPos] ==
                   '-') {
          tandemRepeat secondTandems = findTandemRepeatOfStrInSequence(
              alignObjectA_.seqBase_.seq_, search);
          if ((int)gapIter->second.startPos >= secondTandems.startPos &&
              (int)gapIter->second.startPos + (int)gapIter->second.size - 1 <=
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
  for (uint32_t i = 0; i < alignObjectA_.seqBase_.seq_.length(); i++) {
    if (alignObjectA_.seqBase_.seq_[i] == '-') {
      gap newGap = gap(i, alignObjectB_.seqBase_.seq_.substr(i, 1),
                       alignObjectB_.seqBase_.qual_[i]);
      newGap.ref = true;
      while (alignObjectA_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectB_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectB_.seqBase_.qual_[i + 1];
        i++;
      }
      if (newGap.startPos + newGap.size >=
          alignObjectB_.seqBase_.seq_.length()) {
        if (editTheSame) {
          //editDistance_ += gapScores_.gapRightExtend_ +
                          // gapScores_.gapRightExtend_ * (newGap.size - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size;
        }
        parts_.score_ -= parts_.gapScores_.gapRightOpen_;
        parts_.score_ -= parts_.gapScores_.gapRightExtend_ * (newGap.size - 1);
      } else if (newGap.startPos == 0) {
      	parts_.score_ -= parts_.gapScores_.gapLeftOpen_;
      	parts_.score_ -= parts_.gapScores_.gapLeftExtend_ * (newGap.size - 1);
        if (editTheSame) {
          //editDistance_ += gapScores_.gapLeftOpen_ +
                          // gapScores_.gapLeftExtend_ * (newGap.size - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size;
        }
      } else {
      	parts_.score_ -= parts_.gapScores_.gapOpen_;
      	parts_.score_ -= parts_.gapScores_.gapExtend_ * (newGap.size - 1);
        if (editTheSame) {
          //editDistance_ +=
            //  gapScores_.gapOpen_ + gapScores_.gapExtend_ * (newGap.size - 1);
        } else {
          // std::cout << "adding gap" << std::endl;
          // std::cout << "gapSize: " << newGap.size << std::endl;
          //editDistance_ += newGap.size;
        }
      }
      continue;
    }
    if (alignObjectB_.seqBase_.seq_[i] == '-') {
      gap newGap = gap(i, alignObjectA_.seqBase_.seq_.substr(i, 1),
                       alignObjectA_.seqBase_.qual_[i]);
      newGap.ref = false;
      while (alignObjectB_.seqBase_.seq_[i + 1] == '-') {
        newGap.gapedSequence.append(
            alignObjectA_.seqBase_.seq_.substr(i + 1, 1));
        newGap.size++;
        newGap.summedQuality += alignObjectA_.seqBase_.qual_[i + 1];
        i++;
      }
      if (newGap.startPos + newGap.size >=
          alignObjectB_.seqBase_.seq_.length()) {
        if (editTheSame) {
          //editDistance_ += gapScores_.gapRightOpen_ +
                          // gapScores_.gapRightExtend_ * (newGap.size - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size;
        }
        parts_.score_ -= parts_.gapScores_.gapRightOpen_;
        parts_.score_ -= parts_.gapScores_.gapRightExtend_ * (newGap.size - 1);
      } else if (newGap.startPos == 0) {
      	parts_.score_ -= parts_.gapScores_.gapLeftOpen_;
      	parts_.score_ -= parts_.gapScores_.gapLeftExtend_ * (newGap.size - 1);
        if (editTheSame) {
          //editDistance_ += gapScores_.gapLeftOpen_ +
                          // gapScores_.gapLeftExtend_ * (newGap.size - 1);
        } else if (countEndGaps_) {
          //editDistance_ += newGap.size;
        }
      } else {
      	parts_.score_ -= parts_.gapScores_.gapOpen_;
      	parts_.score_ -= parts_.gapScores_.gapExtend_ * (newGap.size - 1);
        if (editTheSame) {
          //editDistance_ +=
             // gapScores_.gapOpen_ + gapScores_.gapExtend_ * (newGap.size - 1);
        } else {
          // std::cout << "adding gap" << std::endl;
          // std::cout << "gapSize: " << newGap.size << std::endl;
          //editDistance_ += newGap.size;
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
}  // namespace bib
