/*
 * randStrGen.cpp
 *
 *  Created on: Jul 27, 2014
 *      Author: nickhathaway
 */

#include "randStrGen.hpp"

namespace bibseq {
std::string randStrGen::rStr(uint64_t size){
  std::string ans;
  ans.reserve(size);
  while (ans.size() < size) {
  	ans.push_back(charGen_.genObj());
  }
  return ans;
}

VecStr randStrGen::rStrs(uint64_t size, uint32_t strNum){
  std::vector<std::string> ans(strNum);
  std::generate_n(ans.begin(), strNum,
                  [&]() { return rStr(size) ; });
  return ans;
}

std::string randStrGen::rStr(uint64_t minSize, uint64_t maxSize){
	//
	return rStr(rGen_.unifRand(minSize, maxSize + 1));
}

VecStr randStrGen::rStrs(uint64_t minSize, uint64_t maxSize, uint32_t num){
  std::vector<std::string> ans(num);
  std::generate_n(ans.begin(), num,
                  [&]() { return rStr(rGen_.unifRand(minSize, maxSize + 1)) ; });
  return ans;
}
} /* namespace bib */
