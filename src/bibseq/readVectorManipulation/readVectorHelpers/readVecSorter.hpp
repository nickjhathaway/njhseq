#pragma once
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
//template <typename T>
class readVecSorter {

 public:
  template <typename T>
  static void sort(std::vector<T>& vec, bool decending = true) {
    sortReadVector(vec, "totalCount", decending);
    return;
  }

  template <typename T>
  static void sortReadVector(std::vector<T>& vec, const std::string& sortBy, bool decending = true) {
    if (sortBy == "averageError" || sortBy == "totalCount") {
      sortByTotalCountAE<T>(vec, decending);
    } else if (sortBy == "seq") {
      sortBySeq<T>(vec, decending);
    } else if (sortBy == "seqCondensed") {
      sortBySeqCondensed<T>(vec, decending);
    } else if (sortBy == "qualCheck") {
      sortByQualCheck<T>(vec, decending);
    } else if (sortBy == "size") {
      sortBySeqSize<T>(vec, decending);
    } else if (sortBy == "name") {
      sortByName<T>(vec, decending);
    } else if (sortBy == "fraction") {
      sortByFraction<T>(vec, decending);
    } else {
      std::cout << "unrecognized sort option: " << sortBy << ", not sorting"
                << std::endl;
    }
  }
/*
	template <typename T>
	static void sortReadVectorFunc(std::vector<T>& vec,
			std::function<bool(const T & read1, const T & read2)> func,
			bool decending = true) {
		if (decending) {
			std::sort(vec.begin(), vec.end(), func);
		} else {
			std::sort(vec.rbegin(), vec.rend(), func);
		}
	}*/
	template <typename T, typename FUNC>
	static void sortReadVectorFunc(std::vector<T>& vec,
			FUNC func,
			bool decending = true) {
		if (decending) {
			std::sort(vec.begin(), vec.end(), func);
		} else {
			std::sort(vec.rbegin(), vec.rend(), func);
		}
	}

 private:
  // sorting functions
  template <typename T>
  static void sortByQualCheck(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second) -> bool{
      return first.fractionAboveQualCheck_ > second.fractionAboveQualCheck_;
    }, decending);
  }

  template <typename T>
  static void sortByTotalCountAE(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second) -> bool{
  		/*
			if (first.seqBase_.cnt_ == second.seqBase_.cnt_) {
				if (first.averageErrorRate < second.averageErrorRate) {
					return true;
				} else {
					return false;
				}
			} else {
				return first.seqBase_.cnt_ > second.seqBase_.cnt_;
			}*/
  		if (roundDecPlaces(first.seqBase_.cnt_, 2) == roundDecPlaces(second.seqBase_.cnt_, 2) ) {
  			if (roundDecPlaces(first.averageErrorRate, 2)  < roundDecPlaces(second.averageErrorRate, 2) ) {
  				return true;
  			} else {
  				return false;
  			}
  		} else {
  			return roundDecPlaces(first.seqBase_.cnt_, 2)  > roundDecPlaces(second.seqBase_.cnt_, 2) ;
  		}
    }, decending);
  }

  template <typename T>
  static void sortBySeq(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second) {
      return first.seqBase_.seq_ < second.seqBase_.seq_;
    }, decending);
  }

  template <typename T>
  static void sortByName(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      return first.seqBase_.name_ < second.seqBase_.name_;
    }, decending);
  }
  template <typename T>
  static void sortBySeqCondensed(std::vector<T>& vec, bool decending) {
    readVec::allSetCondensedSeq(vec);
    sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      if (first.condensedSeq == second.condensedSeq) {
        return first.seqBase_.seq_ < second.seqBase_.seq_;
      } else {
        return first.condensedSeq < second.condensedSeq;
      }
    }, decending);
  }
  template <typename T>
  static void sortBySeqSize(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      return first.seqBase_.seq_.size() < second.seqBase_.seq_.size();
    }, decending);
  }
  template <typename T>
  static void sortByFraction(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      return first.seqBase_.frac_ > second.seqBase_.frac_;
    }, decending);
  }
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "readVecSorter.cpp"
#endif
