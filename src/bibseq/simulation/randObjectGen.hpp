//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#pragma once
/*
 * randObjectGen.hpp
 *
 *  Created on: Jul 26, 2014
 *      Author: nickhathaway
 */

#include "bibseq/common/allSystemIncludes.h"
#include "bibseq/simulation/randomGenerator.hpp"

namespace bibseq {

template<typename T, typename N>
class randObjectGen {
public:
	//constructor
	randObjectGen(const std::vector<T> & objs):
	objs_(objs), objCounts_(std::vector<N>(objs.size(),1)),
	likelihoods_(createLikelihood(objs_,objCounts_)){
		std::random_device rd;
		mtGen_.seed(rd());
	}
	randObjectGen(const std::vector<T> & objs,
			const std::vector<N> & counts):
					objs_(objs), objCounts_(counts),
	likelihoods_(createLikelihood(objs_, objCounts_)){
		std::random_device rd;
		mtGen_.seed(rd());
	}
private:
	//members
	std::mt19937_64 mtGen_;
	std::vector<T> objs_;
	std::vector<N> objCounts_;
	std::multimap<uint64_t, T, std::less<uint64_t>> likelihoods_;
public:
	//functions
	T genObj(){
	  uint64_t sum = 0;
	  uint64_t rando = mtGen_();
	  for (const auto &likelihood : likelihoods_) {
	    sum += likelihood.first;
	    if (sum > rando) {
	      return likelihood.second;
	    }
	  }
	  return likelihoods_.rbegin()->second;
	}
	std::vector<T> genObjs(uint32_t num){
		std::vector<T> ans(num);
		std::generate_n(ans.begin(), num, [&]() { return genObj();});
		return ans;
	}
	static std::multimap<uint64_t, T, std::less<uint64_t>> createLikelihood(
	    const std::vector<T> &objs, const std::vector<N> &counts){
	  if (counts.size() != objs.size()) {
	    std::cout << "Error in createLikelihood(const std::vector<T> &objs,"
	                 " const std::vector<uint32_t> & counts)" << std::endl;
	    std::cout << "Size of counts differs from size of letters" << std::endl;
	    std::cout << "Size of counts: " << counts.size() << std::endl;
	    std::cout << "Size of objs: " << objs.size() << std::endl;
	    exit(1);
	  }
	  long double countsSum = std::accumulate(counts.begin(), counts.end(), 0);
	  std::multimap<uint64_t, char, std::less<uint64_t>> likelihoods;
	  for (const auto &pos : iter::range(objs.size())) {
	    likelihoods.emplace((std::numeric_limits<uint64_t>::max() / countsSum) * counts[pos], objs[pos]);
	  }
	  return likelihoods;
	}
};



} /* namespace bib */


