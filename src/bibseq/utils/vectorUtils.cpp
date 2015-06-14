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
#include "vectorUtils.hpp"
#include <bibcpp/utils/lexical_cast.hpp>
namespace bibseq {

std::string getSubVector(const std::string& vec, uint32_t start,
                            uint32_t size) {
  std::string ans(vec.begin() + start, vec.begin() + size + start);
  return ans;
}

std::string getSubVector(const std::string& vec, uint32_t start) {
  std::string ans(vec.begin() + start, vec.end());
  return ans;
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
std::vector<uint32_t> getPositionsOfSubStrTarget(const VecStr& vec,
                                                 const std::string& target) {
  uint32_t pos = 0;
  std::vector<uint32_t> positions;
  for (const auto& iter : vec) {
    if (containsSubString(iter, target)) {
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
    if (beginsWith(iter, target)) {
      positions.push_back(pos);
    }
    ++pos;
  }
  return positions;
}

VecStr getUniqueStrings(const VecStr& vec) {
  VecStr ans;
  for (const auto& iter : vec) {
    if (!vectorContains(ans, iter)) {
      ans.push_back(iter);
    }
  }
  return ans;
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
void setDelim(const std::string& newDelim) { bibseq::outDelim = newDelim; }

double getMeanFromVecStr(const VecStr & strNums){
	auto converted = bib::lexical_cast_con<VecStr, std::vector<double>>(strNums);
	return vectorMean(converted);
}
}  // namespace bib
