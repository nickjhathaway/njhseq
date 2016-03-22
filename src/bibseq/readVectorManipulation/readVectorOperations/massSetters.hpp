#pragma once
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
//
//  massSetters.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/7/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//
//#include "bibseq/objects/seqObjects/readObject.hpp"
// mass setters
namespace bibseq {
template<typename T>
const T & getRef(const T & val) {
	return val;
}

template<typename T>
const T & getRef(const std::shared_ptr<const T> & val) {
	return getRef(*val);
}

template<typename T>
const T & getRef(const std::unique_ptr<const T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(T & val) {
	return val;
}

template<typename T>
T & getRef(const std::shared_ptr<T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(const std::unique_ptr<T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(std::shared_ptr<T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(std::unique_ptr<T> & val) {
	return getRef(*val);
}

namespace readVec {



template <class T>
void allSetCondensedSeq(std::vector<T>& reads) {
  for_each(reads, [](T& read) { getRef(read).createCondensedSeq(); });
}

template <class T>
void allSetCondensedSeqCount(std::vector<T>& reads) {
  for_each(reads, [](T& read) { getRef(read).setCondensedCounter(); });
}

template <class T>
void allSetLetterCount(std::vector<T>& reads) {
  for_each(reads, [](T& read) {
  	getRef(read).setLetterCount();
  	getRef(read).counter_.setFractions();
  });
}
template <class T>
void allSetLetterCount(std::vector<T>& reads, const std::vector<char> & alph) {
  for_each(reads, [&](T& read) {
  	getRef(read).setLetterCount(alph);
  	getRef(read).counter_.setFractions(alph);
  });
}

template <class T>
void allUpdateName(std::vector<T>& reads) {
  for_each(reads, [](T& read) { getRef(read).updateName(); });
}

template <class T>
void allSetFractionByTotalCount(std::vector<T>& reads) {
  int count = getTotalReadCount(reads);
  for_each(reads, [&](T& read) { getRef(read).setFractionByCount(count); });
}

template <class T>
void allSetQualCheck(std::vector<T>& reads, int qualCheck) {
  for_each(reads, [&](T& read) { getRef(read).setBaseCountOnQualCheck(qualCheck); });
}



template <typename T>
std::unordered_map<std::string, std::vector<T>> organizeByCondensedSeq(
    std::vector<T>& reads) {
  std::unordered_map<std::string, std::vector<T>> ans;
  readVec::allSetCondensedSeq(reads);
  size_t smallestLength = -1;
  for (const auto& read : reads) {
    if (read.condensedSeq.size() < smallestLength) {
      smallestLength = read.condensedSeq.size();
    }
  }
  for (const auto& read : reads) {
    ans[read.condensedSeq.substr(0, smallestLength)] = {};
  }
  std::string currentShortenCondensedSeq = "";
  for (const auto& read : reads) {
  	currentShortenCondensedSeq = read.condensedSeq.substr(0, smallestLength);
  	auto currentCheck = ans.find(currentShortenCondensedSeq);
    if (currentCheck == ans.end()) {
      ans.insert({currentShortenCondensedSeq, {read}});
    } else {
    	currentCheck->second.emplace_back(read);
    }
  }
  return ans;
}

template <typename T>
std::unordered_map<std::string, std::vector<uint64_t>> organizeByCondensedSeqPositions(
    std::vector<T>& reads) {
  std::unordered_map<std::string, std::vector<uint64_t>> ans;
  readVec::allSetCondensedSeq(reads);

  size_t smallestLength = -1;
  for (const auto& read : reads) {
    if (read.condensedSeq.size() < smallestLength) {
      smallestLength = read.condensedSeq.size();
    }
  }
  for (const auto& read : reads) {
    ans[read.condensedSeq.substr(0, smallestLength)] = {};
  }
  std::string currentShortenCondensedSeq = "";
  for (const auto& readPos : iter::range(reads.size())) {
  	currentShortenCondensedSeq = reads[readPos].condensedSeq.substr(0, smallestLength);
  	auto currentCheck = ans.find(currentShortenCondensedSeq);
    if (currentCheck == ans.end()) {
      ans.insert({currentShortenCondensedSeq, {readPos}});
    } else {
    	currentCheck->second.emplace_back(readPos);
    }
  }
  return ans;
}

template <typename T>
std::unordered_map<std::string, std::vector<uint64_t>> organizeByCondensedSeqPositions(
    std::vector<T>& reads, const std::vector<uint64_t> & positions) {
  std::unordered_map<std::string, std::vector<uint64_t>> ans;
  //readVec::allSetCondensedSeq(reads);
  for(const auto & pos : positions){
  	reads[pos].createCondensedSeq();
  }
  size_t smallestLength = -1;
  for (const auto& readPos : positions) {
    if (reads[readPos].condensedSeq.size() < smallestLength) {
      smallestLength = reads[readPos].condensedSeq.size();
    }
  }
  for (const auto& readPos : positions) {
    ans[reads[readPos].condensedSeq.substr(0, smallestLength)] = {};
  }
  std::string currentShortenCondensedSeq = "";
  for (const auto& readPos : positions) {
  	currentShortenCondensedSeq = reads[readPos].condensedSeq.substr(0, smallestLength);
  	auto currentCheck = ans.find(currentShortenCondensedSeq);
    if (currentCheck == ans.end()) {
      ans.insert({currentShortenCondensedSeq, {readPos}});
    } else {
    	currentCheck->second.emplace_back(readPos);
    }
  }
  return ans;
}
}  // namsespace readVec
}  // namespace bib
