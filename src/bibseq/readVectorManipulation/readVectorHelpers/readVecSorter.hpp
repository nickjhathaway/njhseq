#pragma once
//
//  sorting.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 2/3/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//
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
#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"


namespace bibseq {

class readVecSorter {

 public:
  template <typename T>
  static void sort(std::vector<T>& vec, bool decending = true) {
    sortReadVector(vec, "totalCount", decending);
    return;
  }

	template<typename T>
	static void sortReadVector(std::vector<T>& vec, const std::string& sortBy,
			bool decending = true) {
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
		} else if (sortBy == "reverse") {
			bib::reverse(vec);
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error unrecognized sort option: "
					<< sortBy << ", not sorting" << std::endl;
			throw std::runtime_error { ss.str() };
		}
	}

  template <typename T>
  static void sortReadVector(std::vector<T>& vec,
  		const std::vector<uint32_t> & positions,
  		const std::string& sortBy,
			bool decending = true) {
    if (sortBy == "averageError" || sortBy == "totalCount") {
      sortByTotalCountAE<T>(vec,positions, decending);
    } else if (sortBy == "seq") {
      sortBySeq<T>(vec,positions, decending);
    } else if (sortBy == "seqCondensed") {
      sortBySeqCondensed<T>(vec,positions, decending);
    } else if (sortBy == "qualCheck") {
      sortByQualCheck<T>(vec,positions, decending);
    } else if (sortBy == "size") {
      sortBySeqSize<T>(vec,positions, decending);
    } else if (sortBy == "name") {
      sortByName<T>(vec,positions, decending);
    } else if (sortBy == "fraction") {
      sortByFraction<T>(vec,positions, decending);
    } else if (sortBy == "reverse") {
    	if(1 != positions.size()){
    		auto sortedPositions = positions;
    		bib::sort(sortedPositions);
				for (const auto pos : iter::range(sortedPositions.size() / 2)) {
					//std::cout << pos << ":" << sortedPositions[pos] << ":" << positions[positions.size() - 1 - pos] << std::endl;
					std::iter_swap(vec.begin() + sortedPositions[pos],
							vec.begin() + sortedPositions[sortedPositions.size() - 1 - pos]);
				}
    	}
    } else {
    	std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ": Error unrecognized sort option: " << sortBy << ", not sorting"
                << std::endl;
      throw std::runtime_error{ss.str()};
    }
  }

  template <typename T>
  static void sortReadVectorSimple(std::vector<T>& vec, const std::string& sortBy, bool decending = true) {
    if (sortBy == "totalCount") {
      sortByTotalCount<T>(vec, decending);
    } else if (sortBy == "seq") {
      sortBySeq<T>(vec, decending);
    } else if (sortBy == "size") {
      sortBySeqSize<T>(vec, decending);
    } else if (sortBy == "name") {
      sortByName<T>(vec, decending);
    } else if (sortBy == "fraction") {
      sortByFraction<T>(vec, decending);
    } else {
    	std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ": Error unrecognized sort option: " << sortBy << ", not sorting"
                << std::endl;
      throw std::runtime_error{ss.str()};
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

	template <typename T, typename FUNC>
	static void sortReadVectorFunc(std::vector<T>& vec,
			const std::vector<uint32_t> & positions,
			FUNC func,
			bool decending = true) {
		std::vector<uint32_t> sortedPositions = positions;
		if (decending) {
			std::sort(sortedPositions.begin(), sortedPositions.end(),
					[&func, &vec](uint32_t pos1, uint32_t pos2){
				return func(vec[pos1], vec[pos2]);
			});
		} else {
			std::sort(sortedPositions.rbegin(), sortedPositions.rend(),
					[&func, &vec](const uint32_t pos1, uint32_t pos2){
				return func(vec[pos1], vec[pos2]);
			});
		}
		/**@todo find an more efficient way of doing this that would avoid the copy */
		std::unordered_map<uint32_t, T> tempCpy;
		for(const auto pos : iter::range(sortedPositions.size())){
			if(positions[pos] !=  sortedPositions[pos]){
				tempCpy.emplace(pos, vec[sortedPositions[pos]]);
			}
		}
		for(const auto&  pos : tempCpy){
			vec[positions[pos.first]] = pos.second;
		}
	}

  // sorting functions
  template <typename T>
  static void sortByQualCheck(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second) -> bool{
      return getRef(first).fractionAboveQualCheck_ > getRef(second).fractionAboveQualCheck_;
    }, decending);
  }

  template<typename T>
	static void sortByQualCheck(std::vector<T>& vec,
			const std::vector<uint32_t>& positions, bool decending) {
		sortReadVectorFunc<T>(vec, positions,
				[](const T& first, const T& second) -> bool {
					return getRef(first).fractionAboveQualCheck_ > getRef(second).fractionAboveQualCheck_;
				}, decending);
	}

  template <typename T>
  static void sortByTotalCountAE(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second) -> bool{
  		/*
			if (getSeqBase(first).cnt_ == getSeqBase(second).cnt_) {
				if (getRef(first).averageErrorRate < getRef(second).averageErrorRate) {
					return true;
				} else {
					return false;
				}
			} else {
				return getSeqBase(first).cnt_ > getSeqBase(second).cnt_;
			}*/
  		if (roundDecPlaces(getSeqBase(first).cnt_, 2) == roundDecPlaces(getSeqBase(second).cnt_, 2) ) {
  			if (roundDecPlaces(getRef(first).averageErrorRate, 2)  < roundDecPlaces(getRef(second).averageErrorRate, 2) ) {
  				return true;
  			} else {
  				return false;
  			}
  		} else {
  			return roundDecPlaces(getSeqBase(first).cnt_, 2)  > roundDecPlaces(getSeqBase(second).cnt_, 2) ;
  		}
    }, decending);
  }

  template <typename T>
  static void sortByTotalCountAE(std::vector<T>& vec,
  		const std::vector<uint32_t>& positions, bool decending) {
  	sortReadVectorFunc<T>(vec, positions, [](const T& first, const T& second) -> bool{
  		/*
			if (getSeqBase(first).cnt_ == getSeqBase(second).cnt_) {
				if (getRef(first).averageErrorRate < getRef(second).averageErrorRate) {
					return true;
				} else {
					return false;
				}
			} else {
				return getSeqBase(first).cnt_ > getSeqBase(second).cnt_;
			}*/
  		if (roundDecPlaces(getSeqBase(first).cnt_, 2) == roundDecPlaces(getSeqBase(second).cnt_, 2) ) {
  			if (roundDecPlaces(getSeqBase(first).getAverageErrorRate(), 2)  < roundDecPlaces(getSeqBase(second).getAverageErrorRate(), 2) ) {
  				return true;
  			} else {
  				return false;
  			}
  		} else {
  			return roundDecPlaces(getSeqBase(first).cnt_, 2)  > roundDecPlaces(getSeqBase(second).cnt_, 2) ;
  		}
    }, decending);
  }

  template <typename T>
  static void sortByTotalCount(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second) -> bool{
  		return roundDecPlaces(getSeqBase(first).cnt_, 2)  > roundDecPlaces(getSeqBase(second).cnt_, 2) ;
    }, decending);
  }

  template <typename T>
  static void sortByTotalCount(std::vector<T>& vec,
  		const std::vector<uint32_t>& positions, bool decending) {
  	sortReadVectorFunc<T>(vec,positions, [](const T& first, const T& second) -> bool{
  		return roundDecPlaces(getSeqBase(first).cnt_, 2)  > roundDecPlaces(getSeqBase(second).cnt_, 2) ;
    }, decending);
  }

  template <typename T>
  static void sortBySeq(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second) {
      return getSeqBase(first).seq_ < getSeqBase(second).seq_;
    }, decending);
  }

  template <typename T>
  static void sortBySeq(std::vector<T>& vec,
  		const std::vector<uint32_t>& positions, bool decending) {
  	sortReadVectorFunc<T>(vec,positions, [](const T& first, const T& second) {
      return getSeqBase(first).seq_ < getSeqBase(second).seq_;
    }, decending);
  }

  template <typename T>
  static void sortByName(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      return getSeqBase(first).name_ < getSeqBase(second).name_;
    }, decending);
  }

  template <typename T>
  static void sortByName(std::vector<T>& vec,
  		const std::vector<uint32_t>& positions, bool decending) {
  	sortReadVectorFunc<T>(vec,positions, [](const T& first, const T& second)  {
      return getSeqBase(first).name_ < getSeqBase(second).name_;
    }, decending);
  }

  template <typename T>
  static void sortBySeqCondensed(std::vector<T>& vec, bool decending) {
    readVec::allSetCondensedSeq(vec);
    sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      if (getRef(first).condensedSeq == getRef(second).condensedSeq) {
        return getSeqBase(first).seq_ < getSeqBase(second).seq_;
      } else {
        return getRef(first).condensedSeq < getRef(second).condensedSeq;
      }
    }, decending);
  }

  template <typename T>
  static void sortBySeqCondensed(std::vector<T>& vec,
  		const std::vector<uint32_t>& positions, bool decending) {
    sortReadVectorFunc<T>(vec,positions, [](const T& first, const T& second)  {
      if (condenseSeqSimple(getSeqBase(first).seq_) == condenseSeqSimple(getSeqBase(second).seq_)) {
        return getSeqBase(first).seq_ < getSeqBase(second).seq_;
      } else {
        return condenseSeqSimple(getSeqBase(first).seq_) < condenseSeqSimple(getSeqBase(second).seq_);
      }
    }, decending);
  }

  template <typename T>
  static void sortBySeqSize(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      return getSeqBase(first).seq_.size() < getSeqBase(second).seq_.size();
    }, decending);
  }

  template <typename T>
  static void sortBySeqSize(std::vector<T>& vec,
  		const std::vector<uint32_t>& positions, bool decending) {
  	sortReadVectorFunc<T>(vec,positions, [](const T& first, const T& second)  {
      return getSeqBase(first).seq_.size() < getSeqBase(second).seq_.size();
    }, decending);
  }

  template <typename T>
  static void sortByFraction(std::vector<T>& vec, bool decending) {
  	sortReadVectorFunc<T>(vec, [](const T& first, const T& second)  {
      return getSeqBase(first).frac_ > getSeqBase(second).frac_;
    }, decending);
  }

  template <typename T>
  static void sortByFraction(std::vector<T>& vec,
  		const std::vector<uint32_t>& positions, bool decending) {
  	sortReadVectorFunc<T>(vec,positions, [](const T& first, const T& second)  {
      return getSeqBase(first).frac_ > getSeqBase(second).frac_;
    }, decending);
  }
};
}  // namespace bibseq

