#pragma once
//
//  randomGenerator.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
#include "bibseq/utils.h"
#include <cfloat>
namespace bibseq {

class randomGenerator {

 public:
  randomGenerator() {
    seed();
  }
  randomGenerator(uint64_t givenSeed) {
  	seedNum(givenSeed);
  }
  // Members
  std::mt19937_64 mtGen_;
  // Generate a random number between 0 and 1
  double unifRand() {
  	return mtGen_() /( static_cast<double>(mtGen_.max()) + 1);
  }
  double operator()() { return unifRand(); }

  // return a vector of random numbers between 0 and 1
  std::vector<double> unifRandVector(uint32_t num) {
    std::vector<double> ret(num);
    std::generate(ret.begin(), ret.end(),[&](){
    	return unifRand();} );
    return ret;
  }

  template <typename T>
  T unifRand(T start, T stop) {
    return static_cast<T>((stop - start) * unifRand()) + start;
  }
  template <typename T>
  T unifRandSelection(const std::vector<T> & vec) {
    return vec[unifRand<uint64_t>(0,vec.size())];
  }
  template <typename T>
  std::vector<T> unifRandSelectionVec(const std::vector<T> & vec,
  		uint32_t amt, bool withReplacement) {
  	std::vector<T> ret;
  	std::vector<uint32_t> rSel;
  	if(withReplacement){
  		rSel = unifRandVector<uint32_t>(0,vec.size(), amt);
  	}else{
  		if(amt > vec.size()){
  			std::cout << "Error in unifRandSelectionVec, requesting"
  					" more than is in vec but said without replacement"
  					<< std::endl;
  			exit(1);
  		}else{
  			std::vector<uint32_t> rSelPos(vec.size());
  			iota<uint32_t>(rSel, 0);
  			shuffle(rSel, mtGen_);
  			rSel = getSubVector(rSelPos, 0, amt);
  		}
  	}
		for(const auto & pos : rSel){
			ret.emplace_back(vec[pos]);
		}
    return ret;
  }
  template <typename T>
  std::vector<T> unifRandVector(T start, T stop, int num) {
    std::vector<T> ret(num);
    std::generate_n(ret.begin(), num, [&]() { return unifRand(start, stop); });
    return ret;
  }
  template <typename T>
  std::vector<std::vector<T>> unifRandVecVec(T start, T stop,
  		uint32_t totalNum, uint32_t subNum) {
    std::vector<std::vector<T>> ret(totalNum);
  	std::generate(ret.begin(), ret.end(),[&](){ return unifRandVector<T>(start,stop,subNum);} );
    return ret;
  }
  // Reset the random number generator with the system clock.
  void seed() {
    std::random_device rd;
    mtGen_.seed(rd());
  }
  void seedNum(uint64_t givenSeed) {
  	mtGen_.seed(givenSeed);
  }
};

std::map<int32_t, int32_t, std::greater<int32_t>> generate(randomGenerator& gen,
                                               uint32_t start, uint32_t stop,
                                               uint32_t num,
                                               bool verbose = true);
}  // namespace bib
