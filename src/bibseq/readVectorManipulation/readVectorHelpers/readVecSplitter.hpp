#pragma once
//
//  trimmers.hpp
//
//  Created by Nick Hathaway on 2/3/13.
//
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
#include "njhseq/utils.h"
#include "njhseq/readVectorManipulation/readVectorOperations.h"
#include "njhseq/seqToolsUtils/seqToolsUtils.hpp"


namespace njhseq {

class readVecSplitter {

 public:
  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorOnReadLengthOutliers(const std::vector<T>& vec) {
    std::vector<uint32_t> lengths;
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal =
        splitVectorOnReadLengthOutliersAdd(vec, outliers, splitCount, false);
    return {normal, outliers};
  };
  template <class T>
  static std::vector<T> splitVectorOnReadLengthOutliersAdd(
      const std::vector<T>& vec, std::vector<T>& badReads, uint32_t& splitCount,
      bool mark = true) {
    std::vector<uint32_t> lengths;
    for (const auto& read : vec) {
      lengths.push_back(getSeqBase(read).seq_.length());
    }
    double meanLength = vectorMean(lengths);
    double sdLength = vectorStandardDeviationSamp(lengths);
    std::vector<T> normal;
    for (const auto& iter : vec) {
      if (fabs(iter.seqBase_.seq_.length() - meanLength) < 2 * sdLength) {
        normal.push_back(iter);
      } else {
        badReads.push_back(iter);
        if (mark) {
          getSeqBase(badReads.back()).name_.append("_len>2sdFromMean");
        }
        ++splitCount;
      }
    }
    return normal;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorWithinBasesOfMean(
      const std::vector<T>& vec, int basesWithin) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal = splitVectorWithinBasesOfMeanAdd(
        vec, basesWithin, outliers, splitCount, false);
    return {normal, outliers};
  };

  template <class T>
  static std::vector<T> splitVectorWithinBasesOfMeanAdd(
      const std::vector<T>& vec, int basesWithin, std::vector<T>& badReads,
      uint32_t& splitCount, bool mark = true) {
    std::vector<uint32_t> lengths;
    for (const auto& read : vec) {
      lengths.push_back(getSeqBase(read).seq_.length());
    }
    double meanLength = vectorMean(lengths);
    std::vector<T> normal;
    for (const auto& read : vec) {
      if (fabs(getSeqBase(read).seq_.length() - meanLength) <= basesWithin) {
        normal.push_back(read);
      } else {
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append(
              "_lenMoreThan" + std::to_string(basesWithin) + "FromMean");
        }
        ++splitCount;
      }
    }
    return normal;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorWithinBasesOfGiven(const std::vector<T>& vec, int basesWithin,
                                int given) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal = splitVectorWithinBasesOfGivenAdd(
        vec, basesWithin, given, outliers, splitCount);
    return {normal, outliers};
  };

  template <class T>
  static std::vector<T> splitVectorWithinBasesOfGivenAdd(
      const std::vector<T>& vec, int basesWithin, int given,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    std::vector<T> normal;
    for (const auto& read : vec) {
      if (std::abs(static_cast<int32_t>(getSeqBase(read).seq_.length()) - given) <= basesWithin) {
        normal.push_back(read);
      } else {
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append("_lenMoreThan" +
                                                std::to_string(basesWithin) +
                                                "From" + std::to_string(given));
        }
        ++splitCount;
      }
    }
    return normal;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorBellowLength(
      const std::vector<T>& vec, uint32_t minLength) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal =
        splitVectorBellowLengthAdd(vec, minLength, outliers, splitCount, false);
    return {normal, outliers};
  };

  template <class T>
  static std::vector<T> splitVectorBellowLengthAdd(const std::vector<T>& vec,
                                                   uint32_t minLength,
                                                   std::vector<T>& badReads,
                                                   uint32_t& splitCount,
                                                   bool mark = true) {
    std::vector<T> goodReads;
    for (const auto& read : vec) {
      if (getSeqBase(read).seq_.length() >= minLength) {
        goodReads.push_back(read);
      } else {
        ++splitCount;
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append("_length<" +
                                                std::to_string(minLength));
        }
      }
    }
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnQualCheck(
      const std::vector<T>& vec, uint32_t qualCheck, double cutOff) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;

    std::vector<T> normal =
    		splitVectorOnQualCheckAdd(vec, qualCheck, cutOff, outliers, splitCount, false);
    return {normal, outliers};
  };

  template <class T>
  static std::vector<T> splitVectorOnQualCheckAdd(const std::vector<T>& vec,
  																								uint32_t qualCheck, double cutOff,
                                                   std::vector<T>& badReads,
                                                   uint32_t& splitCount,
                                                   bool mark = true) {
    std::vector<T> goodReads;
    for (const auto& read : vec) {
      if (getSeqBase(read).getQualCheck(qualCheck) >= cutOff) {
        goodReads.push_back(read);
      } else {
        ++splitCount;
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append("_q" + estd::to_string(qualCheck) + "<" +
          		estd::to_string(cutOff));
        }
      }
    }
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorAboveLength(
      const std::vector<T>& vec, uint32_t maxLength) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal =
        splitVectorAboveLengthAdd(vec, maxLength, outliers, splitCount, false);
    return {normal, outliers};
  };
  template <class T>
  static std::vector<T> splitVectorAboveLengthAdd(const std::vector<T>& vec,
                                                  uint32_t maxLength,
                                                  std::vector<T>& badReads,
                                                  uint32_t& splitCount,
                                                  bool mark = true) {
    std::vector<T> normal;
    for (const auto& read : vec) {
      if (getSeqBase(read).seq_.length() <= maxLength) {
        normal.push_back(read);
      } else {
        ++splitCount;
        badReads.push_back(read);
        if (mark) {
          badReads.back()
              .seqBase_.name_.append("_length>" + std::to_string(maxLength));
        }
      }
    }
    return normal;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorBetweenLengths(
      const std::vector<T>& vec, uint32_t minLength, uint32_t maxLength) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal = splitVectorBetweenLengthAdd(
        vec, minLength, maxLength, outliers, splitCount, false);
    return {normal, outliers};
  };
  template <class T>
  static std::vector<T> splitVectorBetweenLengthAdd(
      const std::vector<T>& vec, uint32_t minLength, uint32_t maxLength,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    std::vector<T> goodReads =
        splitVectorAboveLengthAdd(vec, maxLength, badReads, splitCount, mark);
    goodReads = splitVectorBellowLengthAdd(goodReads, minLength, badReads,
                                           splitCount, mark);
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnReadCount(
      const std::vector<T>& reads, uint32_t runCutoff) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans =
        splitVectorOnReadCountAdd(reads, runCutoff, other, splitCount, false);
    return {ans, other};
  };

  template <class T>
  static std::vector<T> splitVectorOnReadCountAdd(const std::vector<T>& reads,
  		uint32_t runCutoff,
                                                  std::vector<T>& badReads,
                                                  uint32_t& splitCount,
                                                  bool mark = true) {
    std::vector<T> ans;
    for (const auto& read : reads) {
      if (getSeqBase(read).cnt_ > runCutoff) {
        ans.push_back(read);
      } else {
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append("_clusterSize<=" +
                                                std::to_string(runCutoff));
        }
        ++splitCount;
      }
    }
    return ans;
  };
  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnReadFraction(
      const std::vector<T>& reads, double fractionCutOff) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorOnReadCountAdd(reads, fractionCutOff, other,
                                                   splitCount, false);
    return {ans, other};
  };

  template <class T>
  static std::vector<T> splitVectorOnReadFractionAdd(
      const std::vector<T>& reads, double fractionCutOff,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    std::vector<T> ans;
    for (const auto& read : reads) {
      if (getSeqBase(read).frac_ >= fractionCutOff) {
        ans.push_back(read);
      } else {
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append("_fraction<=" +
                                                std::to_string(fractionCutOff));
        }
        ++splitCount;
      }
    }
    return ans;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnRemove(
      const std::vector<T>& reads) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans =
        splitVectorOnRemoveAdd(reads, other, splitCount, "nothing", false);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorOnRemoveAdd(
      const std::vector<T>& reads, std::vector<T>& badReads, uint32_t& splitCount,
      const std::string& markWith = "_remove", bool mark = true) {
    std::vector<T> ans;
    for (const auto& read : reads) {
      if (!read.remove) {
        ans.push_back(read);
      } else {
        ++splitCount;
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append(markWith);
        }
      }
    }
    return ans;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorOnNucleotideComposition(const std::vector<T>& reads) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans =
        splitVectorOnNucleotideCompositionAdd(reads, other, splitCount);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorOnNucleotideCompositionAdd(
      const std::vector<T>& reads, std::vector<T>& badReads,
      uint32_t& splitCount) {
    std::vector<T> ans;
    charCounter mainCounter;
    for (const auto& read : reads) {
    	charCounter counter;
    	counter.increaseCountByString(getSeqBase(read).seq_);
      for (const auto& letPos : iter::range(counter.fractions_.size())) {
        mainCounter.fractions_[letPos] += counter.fractions_[letPos];
      }
    }
    for (auto& iter : mainCounter.fractions_) {
      iter = iter / reads.size();
    }
    std::vector<double> differences;
    for (const auto& read : reads) {
    	charCounter counter;
    	counter.increaseCountByString(getSeqBase(read).seq_);
      double difference = 0.00;
      for (const auto& letPos : iter::range(counter.fractions_.size())) {
        difference += fabs(mainCounter.fractions_[letPos] - counter.fractions_[letPos]);
        // std::cout<<difference<<std::endl;
      }
      differences.push_back(difference);
    }
    double stdCalc = vectorStandardDeviationSamp(differences);
    double meanCalc = vectorMean(differences);
    std::vector<T> in;
    std::vector<T> out;
    for (auto& rIter : reads) {
    	charCounter counter;
    	counter.increaseCountByString(rIter.seqBase_.seq_);
      double difference = 0.00;
      for (const auto& letPos : iter::range(counter.fractions_.size())){

      	difference += fabs(mainCounter.fractions_[letPos] - counter.fractions_[letPos]);
        // std::cout<<difference<<std::endl;
      }
      if (difference > (meanCalc + 2 * stdCalc)) {
        badReads.push_back(rIter);
        ++splitCount;
      } else {
        ans.push_back(rIter);
      }
    }
    return ans;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorOnNucleotideComposition(const std::vector<T>& reads,
                                     charCounter mainCounter) {
    std::vector<T> ans;
    std::vector<T> other;
    uint32_t splitCount = 0;
    ans = splitVectorOnNucleotideCompositionAdd(reads, other, mainCounter,
                                                splitCount);
    return {ans, other};
  };

  template <class T>
  static std::vector<T> splitVectorOnNucleotideCompositionAdd(
      const std::vector<T>& reads, std::vector<T>& badReads,
			charCounter mainCounter, uint32_t& splitCount) {
    std::vector<T> ans;
    std::vector<double> differences;
    for (const auto& read : reads) {
    	charCounter counter;
    	counter.increaseCountByString(getSeqBase(read).seq_);
      double difference = 0.00;
      for (const auto& letPos : iter::range(counter.fractions_.size())) {
        difference += fabs(mainCounter.fractions_[letPos] - counter.fractions_[letPos]);
        // std::cout<<difference<<std::endl;
      }
      differences.push_back(difference);
    }
    double stdCalc = vectorStandardDeviationSamp(differences);
    double meanCalc = vectorMean(differences);
    std::vector<T> in;
    std::vector<T> out;
    for (const auto& read : reads) {
    	charCounter counter;
    	counter.increaseCountByString(getSeqBase(read).seq_);
      double difference = 0.00;
      for (const auto& letPos : iter::range(counter.fractions_.size())){

      	difference += fabs(mainCounter.fractions_[letPos] - counter.fractions_[letPos]);
        // std::cout<<difference<<std::endl;
      }
      if (difference > (meanCalc + 2 * stdCalc)) {
        badReads.push_back(read);
        ++splitCount;
      } else {
        ans.push_back(read);
      }
    }
    return ans;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnSeqSimilarity(
      const std::vector<T>& reads, const readObject& compareObject) {
    std::vector<T> ans;
    std::vector<T> other;
    uint32_t splitCount = 0;
    ans =
        splitVectorOnSeqSimilarityAdd(reads, compareObject, other, splitCount);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorOnSeqSimilarityAdd(
      const std::vector<T>& reads, const readObject& compareObject,
      std::vector<T>& badReads, uint32_t& splitCount) {
    std::vector<T> ans;
    uint64_t maxReadLength = 0;
    readVec::getMaxLength(reads, maxReadLength);
    if (compareObject.seqBase_.seq_.size() > maxReadLength) {
      maxReadLength = compareObject.seqBase_.seq_.size();
    }
    // aligner object
    auto scoringMatrixMap = substituteMatrix::createDegenScoreMatrix(1,-1);
    gapScoringParameters gapPars(5, 1, 5, 1, 0, 0);
    aligner alignerObj = aligner(maxReadLength, gapPars,scoringMatrixMap);
    for (const auto& read : reads) {
      alignerObj.alignCache(compareObject, read, false);
      if (alignerObj.parts_.score_ < 0) {
        ++splitCount;
        badReads.push_back(read);
      } else {
        ans.push_back(read);
      }
    }
    return ans;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorSeqContaining(
      const std::vector<T>& reads, const std::string& str, uint32_t occurences) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorSeqContainingAdd(reads, str, occurences,
                                                     other, splitCount, false);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorSeqContainingAdd(
      const std::vector<T>& reads, const std::string& str, uint32_t occurences,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    std::vector<T> goodReads;
    std::string lowerCaseSearch = stringToLowerReturn(str);
    for (const auto& read : reads) {
      if (countOccurences(getSeqBase(read).seq_, str) +
              countOccurences(getSeqBase(read).seq_, lowerCaseSearch) <
          occurences) {
        goodReads.push_back(read);
      } else {
        badReads.push_back(read);
        ++splitCount;
        if (mark) {
          getSeqBase(badReads.back()).name_.append(
              "_contained>" + std::to_string(occurences - 1) + "_" + str);
        }
      }
    }
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorWithNameContaining(const std::vector<T>& reads,
                                const std::string& exclusionString) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorWithNameContainingAdd(
        reads, exclusionString, other, splitCount, false);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorWithNameContainingAdd(
      const std::vector<T>& reads, const std::string& exclusionString,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    std::vector<T> ans;
    for (const auto& read : reads) {
      if (getSeqBase(read).name_.find(exclusionString) == std::string::npos) {
        ans.push_back(read);
      } else {
        ++splitCount;
        badReads.push_back(read);
        if (mark) {
          getSeqBase(badReads.back()).name_.append("_seqNameContains:" +
                                                exclusionString);
        }
      }
    }
    return ans;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnQualityWindow(
      const std::vector<T>& reads, int qualityWindowLength,
      int qualityWindowStep, int qualityWindowThres, bool useClipped) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorOnQualityWindowAdd(
        reads, qualityWindowLength, qualityWindowStep, qualityWindowThres,
        useClipped, other, splitCount, false);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorOnQualityWindowAdd(
      const std::vector<T>& reads, int qualityWindowLength,
      int qualityWindowStep, int qualityWindowThres, bool useClipped,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    bool passed = true;
    std::vector<T> goodReads;
    for (const auto& read : reads) {
    	passed = seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
    	                                        qualityWindowStep, getSeqBase(read).qual_);
      if (passed) {
        goodReads.push_back(read);
      } else {
        badReads.push_back(read);
        ++splitCount;
        if (mark) {
          std::stringstream window;
          window << "wl:" << qualityWindowLength << ",ws:" << qualityWindowStep
                 << ",wt:" << qualityWindowThres;
          getSeqBase(badReads.back()).name_.append("_failedWindowOf_" +
                                                window.str());
        }
      }
    }
    return goodReads;
  };
  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorOnQualityWindowPlusLength(const std::vector<T>& reads,
                                       int qualityWindowLength,
                                       int qualityWindowStep,
                                       int qualityWindowThres, bool useClipped,
                                       int minLen) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorOnQualityWindowPlusLengthAdd(
        reads, qualityWindowLength, qualityWindowStep, qualityWindowThres,
        useClipped, minLen, other, splitCount, false);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorOnQualityWindowPlusLengthAdd(
      const std::vector<T>& reads, int qualityWindowLength,
      int qualityWindowStep, int qualityWindowThres, bool useClipped,
      int minLen, std::vector<T>& badReads, uint32_t& splitCount,
      bool mark = true) {
    std::vector<T> goodReads;
    for (const auto& read : reads) {
      uint32_t windowFailedPos = 0;
      windowFailedPos = seqUtil::checkQualityWindowPos(
                  qualityWindowLength, qualityWindowThres, qualityWindowStep,
									getSeqBase(read).qual_);
      if (windowFailedPos + 1 > minLen) {
        goodReads.push_back(read);
        goodReads.back().setClip(0, windowFailedPos - 1);
      } else {
        badReads.push_back(read);
        ++splitCount;
        if (mark) {
          std::stringstream window;
          window << "wl:" << qualityWindowLength << ",ws:" << qualityWindowStep
                 << ",wt:" << qualityWindowThres;
          getSeqBase(badReads.back()).name_.append("_failedWindowOf_" +
                                                window.str());
        }
      }
    }
    return goodReads;
  };


  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorOnNsPlusLength(const std::vector<T>& reads,
                                       uint32_t minLen) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorOnNsPlusLengthAdd(
        reads, minLen, other, splitCount, false);
    return {ans, other};
  };

  template <class T>
  static std::vector<T> splitVectorOnNsPlusLengthAdd(
      const std::vector<T>& reads,
      uint32_t minLen, std::vector<T>& badReads, uint32_t& splitCount,
      bool mark = true) {
    std::vector<T> goodReads;
    for (const auto& read : reads) {
      auto firstN = getSeqBase(read).seq_.find("N");
      if(firstN > minLen){
      	goodReads.emplace_back(read);
      	if(firstN != std::string::npos){
      		goodReads.back().setClip(0, firstN - 1);
      	}
      }else{
        badReads.emplace_back(read);
        ++splitCount;
        if (mark) {
          std::stringstream window;
          window << "N_before" << minLen;
          getSeqBase(badReads.back()).name_.append(window.str());
        }
      }
    }
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnFlowNoise(
      const std::vector<T>& reads, uint32_t flowCutOff) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans =
        splitVectorOnFlowNoiseAdd(reads, flowCutOff, other, splitCount, false);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorOnFlowNoiseAdd(const std::vector<T>& reads,
                                                  uint32_t flowCutOff,
                                                  std::vector<T>& badReads,
                                                  uint32_t& splitCount,
                                                  bool mark = true) {
    std::vector<T> goodReads;
    for (auto read : reads) {
      if (read.flowNoiseProcess(flowCutOff)) {
        goodReads.push_back(read);
      } else {
        if (mark) {
        	getSeqBase(read).name_.append("_failedFlowNoiseProcessing");
        }
        badReads.push_back(read);
        ++splitCount;
      }
    }
    return goodReads;
  };
  /*
  template<class T>
  static void splitVectorOnQWindowLengthNsAdd(std::vector<T>& reads,int
  minLength, int maxLength, int qualityWindowLength, int qualityWindowStep, int
  qualityWindowThres, bool useClipped, int occurences, uint32_t& badReadsLength,
  uint32_t & containsNs, uint32_t & badQualityWindow, uint32_t& flowNoise,bool
  flowFiltering, int maxFlows, std::pair<std::vector<T>, std::vector<T>> &
  alreadySplitVectors){

      alreadySplitVectors.first=splitVectorSeqContainingAdd(reads, "N",
  occurences, alreadySplitVectors.second,containsNs);
      alreadySplitVectors.first=splitVectorBetweenLengthAdd(alreadySplitVectors.first,
  minLength,maxLength, alreadySplitVectors.second,badReadsLength);
      alreadySplitVectors.first=splitVectorOnQualityWindowPlusLengthAdd(alreadySplitVectors.first,
  qualityWindowLength, qualityWindowStep, qualityWindowThres,
  useClipped,minLength, alreadySplitVectors.second,badQualityWindow);
      if (flowFiltering) {
          alreadySplitVectors.first=splitVectorOnFlowNoiseAdd(alreadySplitVectors.first,
  maxFlows, alreadySplitVectors.second, flowNoise);
      }
      alreadySplitVectors.first=splitVectorBetweenLengthAdd(alreadySplitVectors.first,
  minLength,maxLength, alreadySplitVectors.second,badReadsLength);
  };*/
};
}  // namespace njhseq


