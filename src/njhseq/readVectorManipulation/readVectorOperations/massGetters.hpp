#pragma once
//
//  massGetters.hpp
//
//  Created by Nicholas Hathaway on 11/18/13.
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
namespace njhseq {
namespace readVec {
// getters
template <typename T>
int getTotalReadCount(const std::vector<T>& reads, bool countRemove = false) {
  int count = 0;
  for (const auto& read : reads) {
    if (getRef(read).remove) {
      if (countRemove) {
        count += getSeqBase(read).cnt_;
      }
    } else {
      count += getSeqBase(read).cnt_;
    }
  }
  return count;
}
template<typename T>
void allPrintSeqs(const std::vector<T> & reads, std::ostream & out = std::cout){
	njh::for_each(reads, [&](const T & read){ getSeqBase(read).outPutSeq(out);});
}

template<typename T>
std::vector<T> getSeqsWithNames(const std::vector<T> & reads, const VecStr & names){
	std::vector<T> ans;
	auto checker = [&](const T & read){
		if(njh::in(getSeqBase(read).name_, names)){
			ans.emplace_back(read);
		}
	};
	njh::for_each(reads, checker);
	return ans;
}

template <typename T>
uint32_t getReadVectorSize(const std::vector<T>& reads, bool countRemove = false) {
  uint32_t count = 0;
  if (countRemove) {
    return reads.size();
  } else {
    for (const auto& read : reads) {
      if (!read.remove) {
        count++;
      }
    }
  }
  return count;
}

template <typename T>
double getAvgLength(const std::vector<T>& reads) {
	uint64_t sum = 0;
	for (const auto& read : reads) {
		sum += getSeqBase(read).seq_.size();
	}
	return static_cast<double>(sum)/static_cast<double>(reads.size());
}


template <typename T>
void getMaxLength(const T& read, uint64_t & compare) {
	if (getSeqBase(read).seq_.length() > compare) {
		compare = getSeqBase(read).seq_.length();
	}
}
template <typename T>
void getMaxLength(const std::vector<T>& reads, uint64_t& compare) {
  for (const auto& read : reads) {
  	getMaxLength(read, compare);
  }
}

template <typename T>
uint64_t getMaxLength(const std::vector<T>& reads) {
	uint64_t maxLen = 0;
	getMaxLength(reads, maxLen);
	return maxLen;
}

template <typename T>
void getMinLength(const T& read, uint64_t & compare) {
	if (getSeqBase(read).seq_.length() < compare) {
		compare = getSeqBase(read).seq_.length();
	}
}

template <typename T>
void getMinLength(const std::vector<T>& reads, uint64_t& compare) {
  for (const auto& read : reads) {
  	getMinLength(read, compare);
  }
}

template <typename T>
uint64_t getMinLength(const std::vector<T>& reads) {
	uint64_t compare = std::numeric_limits<uint64_t>::max();
  for (const auto& read : reads) {
  	getMinLength(read, compare);
  }
  return compare;
}

template <typename T>
std::vector<uint64_t> getLengths(const std::vector<T>& seqs) {
	std::vector<uint64_t>  ret;
  for (const auto& seq : seqs) {
  	ret.emplace_back(len(getRef(seq)));
  }
  return ret;
}


template <typename T>
VecStr getNames(const std::vector<T>& reads) {
  VecStr names;
  for (const auto& read : reads) {
    names.emplace_back(getSeqBase(read).name_);
  }
  return names;
}

template<typename T>
VecStr getSeqs(const std::vector<T> & reads){
	VecStr ret;
	for(const auto & seq : reads){
		ret.emplace_back(getSeqBase(seq).seq_);
	}
	return ret;
}

template <typename T>
size_t getReadIndexByName(const std::vector<T>& reads,
                          const std::string& name) {
  size_t index = 0;
  for (const auto& read : reads) {
    if (getSeqBase(read).name_ == name) {
      return index;
    }
    ++index;
  }
  return std::string::npos;
}

template <typename T>
T& getReadByName(std::vector<T>& reads, const std::string& name) {
  size_t index = getReadIndexByName(reads, name);
  if (index != std::string::npos) {
    return reads[index];
  } else {
    return *reads.end();
  }
}
template <typename T>
const T& getReadByName(const std::vector<T>& reads, const std::string& name) {
  size_t index = getReadIndexByName(reads, name);
  if (index != std::string::npos) {
    return reads[index];
  } else {
    return *reads.end();
  }
}

template <typename T>
bool readVectorContainsReadWithName(std::vector<T>& reads,
                                    const std::string& name) {
  return (getReadIndexByName(reads, name) != std::string::npos);
}

template <typename T>
void getCountOfReadNameContaining(const std::vector<T>& reads,
                                  const std::string& contains, int& count) {
  for (const auto& read : reads) {
    if (getSeqBase(read).name_.find(contains) != std::string::npos) {
      ++count;
    }
  }
  return;
}

template <typename T>
void getReadCountOfReadNameContaining(const std::vector<T>& reads,
                                      const std::string& contains, int& count) {
  for (const auto& read : reads) {
    if (getSeqBase(read).name_.find(contains) != std::string::npos) {
      count += getSeqBase(read).cnt_;
    }
  }
  return;
}

template <typename T>
void getCountOfReadSeqContainingExact(const std::vector<T>& reads,
                                      const std::string& contains, int& count) {
  for (const auto& read : reads) {
    if (getSeqBase(read).seq_.find(contains) != std::string::npos) {
      ++count;
    }
  }
  return;
}

template <typename T>
void getReadCountOfReadSeqContainingExact(const std::vector<T>& reads,
                                          const std::string& contains,
                                          int& count) {
  for (const auto& read : reads) {
    if (getSeqBase(read).seq_.find(contains) != std::string::npos) {
      count += getSeqBase(read).cnt_;
    }
  }
  return;
}

/*
template <typename T>
charCounterArray getAverageLetterFraction(const std::vector<T>& reads,
                                       bool letterCounterSet = false) {
  charCounterArray mainCounter;
  for (auto& read : reads) {
    for (const auto& letter : read.counter_.alphabet_) {
      mainCounter.chars_[letter] += read.counter_.chars_[letter];
    }
  }
  mainCounter.resetAlphabet(false);
  mainCounter.setFractions();
  return mainCounter;
}*/

template <typename T>
VecStr translateAllRet(const std::vector<T>& reads,bool complement, bool reverse, size_t start = 0) {
  VecStr ans;
  for (const auto& read : reads) {
    ans.emplace_back(read.translateRet(complement, reverse, start).seqBase_.seq_);
  }
  return ans;
}

template <typename T>
void updateQaulCountsMultiple(const std::vector<T>& reads,
                              std::map<uint32_t, uint32_t>& qualCounts) {
  njh::for_each(reads, [&](const T& read) { read.updateQualCounts(qualCounts); });
  return;
}
template <typename T>
void updateQualCountsMultiple(
    const std::vector<T> reads,
    std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize, std::array<double, 100> errorLookUp) {
  njh::for_each(reads, [&](const T& read) {
    read.updateQualCounts(counts, qualWindowSize, errorLookUp);
  });
  return;
}


template<typename T>
bool checkIfReadVecsAreSame(const std::vector<T> & reads1,
		const std::vector<T> & reads2){
	if(reads1.size() != reads2.size() ){
		std::cout << "sizes don't match, can't be same " << std::endl;
		std::cout << "size of 1: " << reads1.size() << std::endl;
		std::cout << "size of 2: " << reads2.size() << std::endl;
		return false;
	}
	for(const auto readPos : iter::range(len(reads1))){
		if(reads1[readPos].seqBase_.seq_ != reads2[readPos].seqBase_.seq_){
			std::cout << "failed seq_ on read " << readPos << std::endl;
			std::cout << "seq_ of 1: " << reads1[readPos].seqBase_.seq_ << std::endl;
			std::cout << "seq_ of 2: " << reads2[readPos].seqBase_.seq_ << std::endl;
			return false;
		}
		if(reads1[readPos].seqBase_.name_ != reads2[readPos].seqBase_.name_){
			std::cout << "failed name on read " << readPos << std::endl;
			std::cout << "name_ of 1: " << reads1[readPos].seqBase_.name_ << std::endl;
			std::cout << "name_ of 2: " << reads2[readPos].seqBase_.name_ << std::endl;
			return false;
		}
		if(reads1[readPos].seqBase_.qual_ != reads2[readPos].seqBase_.qual_){
			std::cout << "failed qual_ on read " << readPos << std::endl;
			std::cout << "qual_ of 1: " << njh::conToStr(reads1[readPos].seqBase_.qual_,",") << std::endl;
			std::cout << "qual_ of 2: " << njh::conToStr(reads2[readPos].seqBase_.qual_,",") << std::endl;
			return false;
		}
		if(reads1[readPos].seqBase_.frac_ != reads2[readPos].seqBase_.frac_){
			std::cout << "failed frac_ on read " << readPos << std::endl;
			std::cout << "frac_ of 1: " << reads1[readPos].seqBase_.frac_ << std::endl;
			std::cout << "frac_ of 2: " << reads2[readPos].seqBase_.frac_ << std::endl;
			return false;
		}
		if(reads1[readPos].seqBase_.cnt_ != reads2[readPos].seqBase_.cnt_){
			std::cout << "failed cnt_ on read " << readPos << std::endl;
			std::cout << "cnt_ of 1: " << reads1[readPos].seqBase_.cnt_ << std::endl;
			std::cout << "cnt_ of 2: " << reads2[readPos].seqBase_.cnt_ << std::endl;
			return false;
		}
	}
	return true;
}




}  // readVec
}  // namespace njh
