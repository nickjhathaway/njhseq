#pragma once
//
//  trimmers.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 2/3/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"
#include "bibseq/seqToolsUtils/seqToolsUtils.hpp"

namespace bibseq {

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
    for (const auto& iter : vec) {
      lengths.push_back(iter.seqBase_.seq_.length());
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
          badReads.back().seqBase_.name_.append("_len>2sdFromMean");
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
    for (const auto& iter : vec) {
      lengths.push_back(iter.seqBase_.seq_.length());
    }
    double meanLength = vectorMean(lengths);
    std::vector<T> normal;
    for (const auto& iter : vec) {
      if (fabs(iter.seqBase_.seq_.length() - meanLength) <= basesWithin) {
        normal.push_back(iter);
      } else {
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append(
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
    for (const auto& iter : vec) {
      if (std::abs(static_cast<int32_t>(iter.seqBase_.seq_.length()) - given) <= basesWithin) {
        normal.push_back(iter);
      } else {
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append("_lenMoreThan" +
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
      const std::vector<T>& vec, int minLength) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal =
        splitVectorBellowLengthAdd(vec, minLength, outliers, splitCount, false);
    return {normal, outliers};
  };

  template <class T>
  static std::vector<T> splitVectorBellowLengthAdd(const std::vector<T>& vec,
                                                   int minLength,
                                                   std::vector<T>& badReads,
                                                   uint32_t& splitCount,
                                                   bool mark = true) {
    std::vector<T> goodReads;
    for (const auto& iter : vec) {
      if ((int)iter.seqBase_.seq_.length() >= minLength) {
        goodReads.push_back(iter);
      } else {
        ++splitCount;
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append("_length<" +
                                                std::to_string(minLength));
        }
      }
    }
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnQualCheck(
      const std::vector<T>& vec, int qualCheck, double cutOff) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;

    std::vector<T> normal =
    		splitVectorOnQualCheckAdd(vec, qualCheck, cutOff, outliers, splitCount, false);
    return {normal, outliers};
  };

  template <class T>
  static std::vector<T> splitVectorOnQualCheckAdd(const std::vector<T>& vec,
  																								int qualCheck, double cutOff,
                                                   std::vector<T>& badReads,
                                                   uint32_t& splitCount,
                                                   bool mark = true) {
    std::vector<T> goodReads;
    for (const auto& iter : vec) {
      if (iter.fractionAboveQualCheck_ >= cutOff) {
        goodReads.push_back(iter);
      } else {
        ++splitCount;
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append("_q" + to_string(qualCheck) + "<" +
                                                std::to_string(cutOff));
        }
      }
    }
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorAboveLength(
      const std::vector<T>& vec, int maxLength) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal =
        splitVectorAboveLengthAdd(vec, maxLength, outliers, splitCount, false);
    return {normal, outliers};
  };
  template <class T>
  static std::vector<T> splitVectorAboveLengthAdd(const std::vector<T>& vec,
                                                  int maxLength,
                                                  std::vector<T>& badReads,
                                                  uint32_t& splitCount,
                                                  bool mark = true) {
    std::vector<T> normal;
    for (const auto& iter : vec) {
      if ((int)iter.seqBase_.seq_.length() <= maxLength) {
        normal.push_back(iter);
      } else {
        ++splitCount;
        badReads.push_back(iter);
        if (mark) {
          badReads[badReads.size() - 1]
              .seqBase_.name_.append("_length>" + std::to_string(maxLength));
        }
      }
    }
    return normal;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorBetweenLengths(
      const std::vector<T>& vec, int minLength, int maxLength) {
    std::vector<T> outliers;
    uint32_t splitCount = 0;
    std::vector<T> normal = splitVectorBetweenLengthAdd(
        vec, minLength, maxLength, outliers, splitCount, false);
    return {normal, outliers};
  };
  template <class T>
  static std::vector<T> splitVectorBetweenLengthAdd(
      const std::vector<T>& vec, int minLength, int maxLength,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    std::vector<T> goodReads =
        splitVectorAboveLengthAdd(vec, maxLength, badReads, splitCount, mark);
    goodReads = splitVectorBellowLengthAdd(goodReads, minLength, badReads,
                                           splitCount, mark);
    return goodReads;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorOnReadCount(
      const std::vector<T>& reads, int runCutoff) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans =
        splitVectorOnReadCountAdd(reads, runCutoff, other, splitCount, false);
    return {ans, other};
  };

  template <class T>
  static std::vector<T> splitVectorOnReadCountAdd(const std::vector<T>& reads,
                                                  int runCutoff,
                                                  std::vector<T>& badReads,
                                                  uint32_t& splitCount,
                                                  bool mark = true) {
    std::vector<T> ans;
    for (const auto& iter : reads) {
      if (iter.seqBase_.cnt_ > runCutoff) {
        ans.push_back(iter);
      } else {
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append("_clusterSize<=" +
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
    for (const auto& iter : reads) {
      if (iter.seqBase_.frac_ >= fractionCutOff) {
        ans.push_back(iter);
      } else {
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append("_fraction<=" +
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
    for (const auto& iter : reads) {
      if (!iter.remove) {
        ans.push_back(iter);
      } else {
        ++splitCount;
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append(markWith);
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
    letterCounter mainCounter = letterCounter();
    for (auto& rIter : reads) {
      letterCounter counter(rIter.seqBase_.seq_);
      for (const auto& iter : counter.fractions) {
        mainCounter.fractions[iter.first] += iter.second;
      }
    }
    for (auto& iter : mainCounter.fractions) {
      iter.second = iter.second / reads.size();
    }
    std::vector<double> differences;
    for (auto& rIter : reads) {
      letterCounter counter(rIter.seqBase_.seq_);
      double difference = 0.00;
      for (const auto& fIter : counter.fractions) {

        difference += fabs(mainCounter.fractions[fIter.first] - fIter.second);
        // std::cout<<difference<<std::endl;
      }
      differences.push_back(difference);
    }
    double stdCalc = vectorStandardDeviationSamp(differences);
    double meanCalc = vectorMean(differences);
    std::vector<T> in;
    std::vector<T> out;
    for (auto& rIter : reads) {
      letterCounter counter(rIter.seqBase_.seq_);
      double difference = 0.00;

      for (const auto& fIter : counter.fractions) {

        difference += fabs(mainCounter.fractions[fIter.first] - fIter.second);
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
                                     letterCounter mainCounter) {
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
      letterCounter mainCounter, uint32_t& splitCount) {
    std::vector<T> ans;
    std::vector<double> differences;
    for (auto& rIter : reads) {
      letterCounter counter(rIter.seqBase_.seq_);
      double difference = 0.00;
      for (const auto& fIter : counter.fractions) {

        difference += fabs(mainCounter.fractions[fIter.first] - fIter.second);
        // std::cout<<difference<<std::endl;
      }
      differences.push_back(difference);
    }
    double stdCalc = vectorStandardDeviationSamp(differences);
    double meanCalc = vectorMean(differences);
    std::vector<T> in;
    std::vector<T> out;
    for (auto& rIter : reads) {
      letterCounter counter(rIter.seqBase_.seq_);
      double difference = 0.00;

      for (const auto& fIter : counter.fractions) {

        difference += fabs(mainCounter.fractions[fIter.first] - fIter.second);
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
    kmerMaps emptyMaps;
    gapScoringParameters gapPars(7, 1, 7, 1, 0, 0);
    aligner alignerObj = aligner(maxReadLength, gapPars,scoringMatrixMap, emptyMaps, 20, 15, 5,
                                 false);

    for (auto& rIter : reads) {
      alignerObj.alignVec(compareObject, rIter, false);
      if (alignerObj.parts_.score_ < 0) {
        ++splitCount;
        badReads.push_back(rIter);
      } else {
        ans.push_back(rIter);
      }
    }
    return ans;
  };

  template <class T>
  static std::pair<std::vector<T>, std::vector<T>> splitVectorSeqContaining(
      const std::vector<T>& reads, const std::string& str, int occurences) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorSeqContainingAdd(reads, str, occurences,
                                                     other, splitCount, false);
    return {ans, other};
  };
  template <class T>
  static std::vector<T> splitVectorSeqContainingAdd(
      const std::vector<T>& reads, const std::string& str, int occurences,
      std::vector<T>& badReads, uint32_t& splitCount, bool mark = true) {
    std::vector<T> goodReads;
    std::string lowerCaseSearch = stringToLowerReturn(str);
    for (const auto& iter : reads) {
      if (countOccurences(iter.seqBase_.seq_, str) +
              countOccurences(iter.seqBase_.seq_, lowerCaseSearch) <
          occurences) {
        goodReads.push_back(iter);
      } else {
        badReads.push_back(iter);
        ++splitCount;
        if (mark) {
          badReads.back().seqBase_.name_.append(
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
    for (const auto& iter : reads) {
      if (iter.seqBase_.name_.find(exclusionString) == std::string::npos) {
        ans.push_back(iter);
      } else {
        ++splitCount;
        badReads.push_back(iter);
        if (mark) {
          badReads.back().seqBase_.name_.append("_seqNameContains:" +
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
    for (const auto& iter : reads) {
      if (useClipped) {
        passed =
            seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
                                        qualityWindowStep, iter.qualityClip);
      } else {
        passed =
            seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
                                        qualityWindowStep, iter.seqBase_.qual_);
      }
      if (passed) {
        goodReads.push_back(iter);
      } else {
        badReads.push_back(iter);
        ++splitCount;
        if (mark) {
          std::stringstream window;
          window << "wl:" << qualityWindowLength << ",ws:" << qualityWindowStep
                 << ",wt:" << qualityWindowThres;
          badReads.back().seqBase_.name_.append("_failedWindowOf_" +
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
    for (const auto& iter : reads) {
      uint32_t windowFailedPos = 0;
      if (useClipped) {
        windowFailedPos = seqUtil::checkQualityWindowPos(
            qualityWindowLength, qualityWindowThres, qualityWindowStep,
            iter.qualityClip);
      } else {
        windowFailedPos = seqUtil::checkQualityWindowPos(
            qualityWindowLength, qualityWindowThres, qualityWindowStep,
            iter.seqBase_.qual_);
      }
      if (windowFailedPos + 1 > minLen) {
        goodReads.push_back(iter);
        goodReads.back().setClip(0, windowFailedPos - 1);
      } else {
        badReads.push_back(iter);
        ++splitCount;
        if (mark) {
          std::stringstream window;
          window << "wl:" << qualityWindowLength << ",ws:" << qualityWindowStep
                 << ",wt:" << qualityWindowThres;
          badReads.back().seqBase_.name_.append("_failedWindowOf_" +
                                                window.str());
        }
      }
    }
    return goodReads;
  };


  template <class T>
  static std::pair<std::vector<T>, std::vector<T>>
  splitVectorOnNsPlusLength(const std::vector<T>& reads,
                                       int minLen) {
    std::vector<T> other;
    uint32_t splitCount = 0;
    std::vector<T> ans = splitVectorOnNsPlusLengthAdd(
        reads, minLen, other, splitCount, false);
    return {ans, other};
  };

  template <class T>
  static std::vector<T> splitVectorOnNsPlusLengthAdd(
      const std::vector<T>& reads,
      int minLen, std::vector<T>& badReads, uint32_t& splitCount,
      bool mark = true) {
    std::vector<T> goodReads;
    for (const auto& read : reads) {
      auto firstN = read.seqBase_.seq_.find("N");
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
          badReads.back().seqBase_.name_.append(window.str());
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
    for (auto iter : reads) {
      if (iter.flowNoiseProcess(flowCutOff)) {
        goodReads.push_back(iter);
      } else {
        if (mark) {
          iter.seqBase_.name_.append("_failedFlowNoiseProcessing");
        }
        badReads.push_back(iter);
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
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "readVecSplitter.cpp"
#endif
