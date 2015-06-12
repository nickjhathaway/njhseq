#pragma once
/*
 * distCalc.hpp
 *
 *  Created on: May 25, 2015
 *      Author: nick
 */

#include "bibseq/utils.h"

namespace bibseq {


template<typename T, typename RET, typename... Args>
void paritialDis(const std::vector<T> & vec,
		std::vector<std::pair<uint32_t, uint32_t>> inds,
		std::vector<std::vector<RET>> & ret,
		std::function<RET(const T & e1, const T& e2, Args... )> func,
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
  if(numThreads < 2){
  	paritialDis(vec, indices, ret, func, args...);
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
    	ts.push_back(std::thread(paritialDis<T,RET, Args...>, std::cref(vec),
    			indsSplit[tNum], std::ref(ret), func,
    			std::cref(args)...));
  	}
  	for(auto & t : ts){
  		t.join();
  	}
  }
  return ret;
}

template<typename T, typename RET, typename... Args>
void paritialDisNonConst(const std::vector<T> & vec,
		std::vector<std::pair<uint32_t, uint32_t>> inds,
		std::vector<std::vector<RET>> & ret,
		std::function<RET(const T & e1, const T& e2, Args... )> func,
		Args&... args){
  for(const auto & i : inds){
  	ret[i.first][i.second] = func(vec[i.first], vec[i.second], args...);
  }
};

template<typename T, typename RET, typename... Args>
std::vector<std::vector<RET>> getDistanceNonConst(const std::vector<T> & vec,
		uint32_t numThreads, std::function<RET(const T & e1, const T& e2, Args... )> func,
		Args&... args){
	std::vector<std::vector<RET>> ret;
	std::vector<std::pair<uint32_t, uint32_t>> indices;
  for(const auto & pos : iter::range(vec.size())){
  	ret.emplace_back(std::vector<RET>(pos));
  	for(const auto & secondPos : iter::range(pos)){
  		indices.emplace_back(pos, secondPos);
  	}
  }
  if(numThreads < 2){
  	paritialDis(vec, indices, ret, func, args...);
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
  if(numThreads < 2){
  	paritialDis(vec, indices, ret, func, args...);
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


