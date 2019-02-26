#include "vectorUtils.hpp"
#include <njhcpp/utils/lexical_cast.hpp>
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

std::string getSubVector(const std::string& vec, uint32_t start,
		uint32_t size) {
	if(start > vec.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "start: " << start << " greater than string size: " << vec.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	return std::numeric_limits<uint32_t>::max() != size ? std::string(vec.begin() + start, size + start > vec.size() ? vec.end() : vec.begin() + size + start): std::string(vec.begin() + start, vec.end());
}

std::string getSubVector(const std::string& vec, uint32_t start) {
	if(start > vec.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "start: " << start << " greater than vector size: " << vec.size() << "\n";
		throw std::runtime_error{ss.str()};
	}
	return std::string(vec.begin() + start, vec.end());;
}

const VecStr fastPermuteVectorOneLength(std::string vec) {
  VecStr ans;
  int numOfPermutes = Factorial((int)vec.size());
  ans.reserve(numOfPermutes);
  do {
    ans.push_back(vec);
  } while (std::next_permutation(vec.begin(), vec.end()));
  return ans;
}

std::vector<std::vector<int>> permuteLengthN(int N) {
  std::vector<int> seq(N);
  std::iota(seq.begin(), seq.end(), 1);
  // int start = (int)time(NULL);
  std::vector<std::vector<int>> combinations = fastPermuteVectorOneLength(seq);
  // std::cout << "Permute " << getDuration(start) << std::endl;
  std::cout << combinations.size() << std::endl;
  return combinations;
}

std::vector<uint32_t> getPositionsOfTarget(const VecStr& vec,
                                           const std::string& target) {
  uint32_t pos = 0;
  std::vector<uint32_t> positions;
  for (const auto& iter : vec) {
    if (iter == target) {
      positions.push_back(pos);
    }
    ++pos;
  }
  return positions;
}

std::vector<uint32_t> getPositionsContainingPattern(const VecStr& vec,
		const std::regex & pattern){
	std::vector<uint32_t> positions;
	for(const auto pos : iter::range(vec.size())){
		std::smatch match;
		if(std::regex_search(vec[pos], match, pattern) ){
			positions.emplace_back(pos);
		}
	}
	return positions;
}

std::vector<uint32_t> getPositionsMatchingPattern(const VecStr& vec,
		const std::regex & pattern){
	std::vector<uint32_t> positions;
	for(const auto pos : iter::range(vec.size())){
		std::smatch match;
		if(std::regex_match(vec[pos], match, pattern) ){
			positions.emplace_back(pos);
		}
	}
	return positions;
}



std::vector<uint32_t> getPositionsOfSubStrTarget(const VecStr& vec,
		const std::string& target) {
	uint32_t pos = 0;
	std::vector<uint32_t> positions;
	for (const auto& iter : vec) {
		if (njh::containsSubString(iter, target)) {
			positions.push_back(pos);
		}
		++pos;
	}
	return positions;
}

std::vector<uint32_t> getPositionsOfTargetStartsWith(
    const VecStr& vec, const std::string& target) {
  uint32_t pos = 0;
  std::vector<uint32_t> positions;
  for (const auto& iter : vec) {
    if (njh::beginsWith(iter, target)) {
      positions.push_back(pos);
    }
    ++pos;
  }
  return positions;
}

VecStr getUniqueStrings(const VecStr& vec) {
  std::set<std::string> ret;
  for (const auto& str : vec) {
  	ret.emplace(str);
  }
  return VecStr{ret.begin(), ret.end()};
}

VecStr getStringsContains(const VecStr& vec, const std::string& contains) {
  VecStr ans;
  for (const auto& iter : vec) {
    if (iter.find(contains) != std::string::npos) {
      ans.push_back(iter);
    }
  }
  return ans;
}


double getMeanFromVecStr(const VecStr & strNums){
	auto converted = njh::lexical_cast_con<VecStr, std::vector<double>>(strNums);
	return vectorMean(converted);
}
}  // namespace njhseq
