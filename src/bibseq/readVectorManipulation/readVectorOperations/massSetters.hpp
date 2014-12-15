#pragma once
//
//  massSetters.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/7/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

// mass setters
namespace bibseq {
namespace readVec {
template <class T>
void allSetCondensedSeq(std::vector<T>& reads) {
  for_each(reads, [](T& read) { read.createCondensedSeq(); });
}

template <class T>
void allSetCondensedSeqCount(std::vector<T>& reads) {
  for_each(reads, [](T& read) { read.setCondensedCounter(); });
}

template <class T>
void allSetLetterCount(std::vector<T>& reads) {
  for_each(reads, [](T& read) {
    read.setLetterCount();
    read.counter_.setFractions();
  });
}
template <class T>
void allSetLetterCount(std::vector<T>& reads, const std::vector<char> & alph) {
  for_each(reads, [&](T& read) {
    read.setLetterCount(alph);
    read.counter_.setFractions(alph);
  });
}

template <class T>
void allUpdateName(std::vector<T>& reads) {
  for_each(reads, [](T& read) { read.updateName(); });
}

template <class T>
void allSetFractionByTotalCount(std::vector<T>& reads) {
  int count = getTotalReadCount(reads);
  for_each(reads, [&](T& read) { read.setFractionByCount(count); });
}

template <class T>
void allSetQualCheck(std::vector<T>& reads, int qualCheck) {
  for_each(reads, [&](T& read) { read.setBaseCountOnQualCheck(qualCheck); });
}

template <class T>
void allSetFractionName(std::vector<T>& reads) {
  for_each(reads, [](T& read) { read.setFractionName(); });
}
template <class T>
void allSetCumulativeFractionName(std::vector<T>& reads) {
  for_each(reads, [](T& read) { read.setCumulativeFractionName(); });
}

template <class T>
void allSetNormalizedFractionName(std::vector<T>& reads) {
  for_each(reads, [](T& read) { read.setNormalizedFractionName(); });
}

template <class T>
void allSetNormalizedFraction(std::vector<T>& reads,
                              size_t totalNumberOfSamples) {
  for (auto& read : reads) {
    read.setNormalizedFraction(totalNumberOfSamples);
  }
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
