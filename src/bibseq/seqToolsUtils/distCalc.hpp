#pragma once
/*
 * distCalc.hpp
 *
 *  Created on: May 25, 2015
 *      Author: nick
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

namespace bibseq {


template<typename T, typename RET, typename... Args>
void paritialDis(const std::vector<T> & vec,
		std::vector<std::pair<uint32_t, uint32_t>> inds,
		std::vector<std::vector<RET>> & ret,
		std::function<RET(const T & e1, const T& e2, const Args &... )> func,
		const Args&... args){
  for(const auto & i : inds){
  	ret[i.first][i.second] = func(vec[i.first], vec[i.second], args...);
  }
};

template<typename T, typename RET, typename... Args>
std::vector<std::vector<RET>> getDistance(const std::vector<T> & vec,
		uint32_t numThreads, std::function<RET(const T & e1, const T& e2, Args... )> func,
		const Args&... args){
	std::vector<std::vector<RET>> ret;
	std::vector<std::pair<uint32_t, uint32_t>> indices;
  for(const auto & pos : iter::range(vec.size())){
  	ret.emplace_back(std::vector<RET>(pos));
  	for(const auto & secondPos : iter::range(pos)){
  		indices.emplace_back(pos, secondPos);
  	}
  }
  if(numThreads < 2 || numThreads >= vec.size()){
  	paritialDis(vec, indices, ret, func, args...);
  }else{
  	std::vector<std::thread> threads;
  	uint32_t step = std::round(indices.size()/static_cast<double>(numThreads));
  	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> indsSplit;
  	for(const auto & tNum : iter::range(numThreads - 1)){
  		std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + tNum * step,
  			indices.begin() + (tNum + 1)*step};
  		indsSplit.emplace_back(temp);
  	}
  	std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + (numThreads - 1) * step,
  	  			indices.end()};
  	indsSplit.emplace_back(temp);
  	for(const auto & tNum : iter::range(numThreads)){
  		threads.push_back(std::thread(paritialDis<T,RET, Args...>, std::cref(vec),
    			indsSplit[tNum], std::ref(ret), func,
    			std::cref(args)...));
  	}
  	for(auto & t : threads){
  		t.join();
  	}
  }
  return ret;
}

template<typename T, typename RET, typename... Args>
void paritialDisNonConst(const std::vector<T> & vec,
		std::vector<std::pair<uint32_t, uint32_t>> inds,
		std::vector<std::vector<RET>> & ret,
		std::function<RET(const T & e1, const T& e2, Args&... )> func,
		Args&... args){
  for(const auto & i : inds){
  	ret[i.first][i.second] = func(vec[i.first], vec[i.second], args...);
  }
};

template<typename T, typename RET, typename... Args>
std::vector<std::vector<RET>> getDistanceNonConst(const std::vector<T> & vec,
		uint32_t numThreads, std::function<RET(const T & e1, const T& e2, Args&... )> func,
		Args&... args){
	std::vector<std::vector<RET>> ret;
	std::vector<std::pair<uint32_t, uint32_t>> indices;
  for(const auto & pos : iter::range(vec.size())){
  	ret.emplace_back(std::vector<RET>(pos));
  	for(const auto & secondPos : iter::range(pos)){
  		indices.emplace_back(pos, secondPos);
  	}
  }
  if(numThreads < 2 || numThreads >= vec.size()){
  	paritialDisNonConst(vec, indices, ret, func, args...);
  }else{
  	std::vector<std::thread> ts;
  	uint32_t step = std::round(indices.size()/static_cast<double>(numThreads));
  	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> indsSplit;
  	for(const auto & tNum : iter::range(numThreads - 1)){
  		std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + tNum * step,
  			indices.begin() + (tNum + 1)*step};
  		indsSplit.emplace_back(temp);
  	}
  	std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + (numThreads - 1) * step,
  	  			indices.end()};
  	indsSplit.emplace_back(temp);
  	for(const auto & tNum : iter::range(numThreads)){
    	ts.push_back(std::thread(paritialDisNonConst<T,RET, Args...>, std::cref(vec),
    			indsSplit[tNum], std::ref(ret), func,
    			std::ref(args)...));
  	}
  	for(auto & t : ts){
  		t.join();
  	}
  }
  return ret;
}

template<typename T, typename RET, typename... Args>
void paritialDisCopy(const std::vector<T> & vec,
		std::vector<std::pair<uint32_t, uint32_t>> inds,
		std::vector<std::vector<RET>> & ret,
		std::function<RET(const T & e1, const T& e2, Args... )> func,
		Args... args){
  for(const auto & i : inds){
  	ret[i.first][i.second] = func(vec[i.first], vec[i.second], args...);
  }
};

template<typename T, typename RET, typename... Args>
std::vector<std::vector<RET>> getDistanceCopy(const std::vector<T> & vec,
		uint32_t numThreads, std::function<RET(const T & e1, const T& e2, Args... )> func,
		Args... args){
	std::vector<std::vector<RET>> ret;
	std::vector<std::pair<uint32_t, uint32_t>> indices;
  for(const auto & pos : iter::range(vec.size())){
  	ret.emplace_back(std::vector<RET>(pos));
  	for(const auto & secondPos : iter::range(pos)){
  		indices.emplace_back(pos, secondPos);
  	}
  }
  if(numThreads < 2 || numThreads >= vec.size()){
  	paritialDisCopy(vec, indices, ret, func, args...);
  }else{
  	std::vector<std::thread> ts;
  	uint32_t step = std::round(indices.size()/static_cast<double>(numThreads));
  	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> indsSplit;
  	for(const auto & tNum : iter::range(numThreads - 1)){
  		std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + tNum * step,
  			indices.begin() + (tNum + 1)*step};
  		indsSplit.emplace_back(temp);
  	}
  	std::vector<std::pair<uint32_t, uint32_t>> temp {indices.begin() + (numThreads - 1) * step,
  	  			indices.end()};
  	indsSplit.emplace_back(temp);
  	for(const auto & tNum : iter::range(numThreads)){
    	ts.push_back(std::thread(paritialDisCopy<T,RET, Args...>, std::cref(vec),
    			indsSplit[tNum], std::ref(ret), func,
    			args...));
  	}
  	for(auto & t : ts){
  		t.join();
  	}
  }
  return ret;
}


} /* namespace bibseq */


