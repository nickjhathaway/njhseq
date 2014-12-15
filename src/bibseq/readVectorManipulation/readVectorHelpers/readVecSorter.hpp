#pragma once
//
//  sorting.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 2/3/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"

// template <class T>
namespace bibseq {

class readVecSorter {

 public:
  template <typename T>
  static void sort(std::vector<T>& vec, bool decending = true) {
    sortReadVector(vec, "totalCount", decending);
    return;
  }

  template <class T>
  static void sortReadVector(std::vector<T>& vec, const std::string& sortBy,
                             bool decending = true) {
    if (sortBy == "averageError" || sortBy == "totalCount") {
      sortByTotalCountAE(vec, decending);
    } else if (sortBy == "seq") {
      sortBySeq(vec, decending);
    } else if (sortBy == "seqClip") {
      sortBySeqClip(vec, decending);
    } else if (sortBy == "seqCondensed") {
      sortBySeqCondensed(vec, decending);
    } else if (sortBy == "qualCheck") {
      sortByQualCheck(vec, decending);
    } else if (sortBy == "size") {
      sortBySeqSize(vec, decending);
    } else if (sortBy == "name") {
      sortByName(vec, decending);
    } else if (sortBy == "fraction") {
      sortByFraction(vec, decending);
    } else if (sortBy == "cumulativeFraction") {
      sortByCumulativeFraction(vec, decending);
    } else if (sortBy == "normalizedFraction") {
      sortByNormalizedFraction(vec, decending);
    } else {
      std::cout << "unrecognized sort option: " << sortBy << ", not sorting"
                << std::endl;
    }
  }

 private:
  // sorting functions
  template <class T>
  static void sortByQualCheck(std::vector<T>& vec, bool decending) {
    qualCheckComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortByTotalCountAE(std::vector<T>& vec, bool decending) {
    totalCountAEComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortBySeq(std::vector<T>& vec, bool decending) {
    seqComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortBySeqClip(std::vector<T>& vec, bool decending) {
    seqClipComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortByName(std::vector<T>& vec, bool decending) {
    nameComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortBySeqCondensed(std::vector<T>& vec, bool decending) {
    readVec::allSetCondensedSeq(vec);
    seqCondensedComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortBySeqSize(std::vector<T>& vec, bool decending) {
    seqSizeComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortByFraction(std::vector<T>& vec, bool decending) {
    fractionComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }

  template <class T>
  static void sortByCumulativeFraction(std::vector<T>& vec, bool decending) {
    cumulativeFractionComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }
  template <class T>
  static void sortByNormalizedFraction(std::vector<T>& vec, bool decending) {
    normalizedFractionComparerStr<T> comparer;
    if (decending) {
      std::sort(vec.begin(), vec.end(), comparer);
    } else {
      std::sort(vec.rbegin(), vec.rend(), comparer);
    }
  }

  // structs
  template <class T>
  struct seqSizeComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.seqBase_.seq_.size() < second.seqBase_.seq_.size();
    }
  };
  template <class T>
  struct nameComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.seqBase_.name_ < second.seqBase_.name_;
    }
  };
  template <class T>
  struct seqComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.seqBase_.seq_ < second.seqBase_.seq_;
    }
  };
  template <class T>
  struct seqClipComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.seqClip < second.seqClip;
    }
  };
  template <class T>
  struct seqCondensedComparerStr {
    bool operator()(const T& first, const T& second) const {
      if (first.condensedSeq == second.condensedSeq) {
        return first.seqBase_.seq_ < second.seqBase_.seq_;
      } else {
        return first.condensedSeq < second.condensedSeq;
      }
    }
  };
  template <class T>
  struct qualCheckComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.fractionAboveQualCheck_ > second.fractionAboveQualCheck_;
    }
  };
  template <class T>
  struct fractionComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.seqBase_.frac_ > second.seqBase_.frac_;
    }
  };

  template <class T>
  struct cumulativeFractionComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.cumulativeFraction > second.cumulativeFraction;
    }
  };
  template <class T>
  struct normalizedFractionComparerStr {
    bool operator()(const T& first, const T& second) const {
      return first.normalizedFraction > second.normalizedFraction;
    }
  };

  template <class T>
  struct totalCountAEComparerStr {
    bool operator()(const T& first, const T& second) const {
      if (first.seqBase_.cnt_ == second.seqBase_.cnt_) {
        return first.averageErrorRate < second.averageErrorRate;
      } else {
        return first.seqBase_.cnt_ > second.seqBase_.cnt_;
      }
    }
  };
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "readVecSorter.cpp"
#endif
