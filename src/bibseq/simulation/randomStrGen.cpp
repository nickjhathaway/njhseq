//
//  randomStrGen.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/27/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "randomStrGen.hpp"

#include "bibseq/simulation/simulationCommon.hpp"
#include "bibseq/objects/counters/letterCounter.hpp"
namespace bibseq {
namespace simulation {

//even likelihood for all chars generators
std::multimap<double, std::string, std::less<double>> getEvenLikelihood(
    const std::vector<char> &letters) {
  std::multimap<double, std::string, std::less<double>> likelihoods;
  for (const auto &let : letters) {
    likelihoods.emplace(1.0 / letters.size(), std::string(1, let));
  }
  return likelihoods;
}
std::multimap<double, std::string, std::less<double>> createLikelihood(
    const std::vector<char> &letters, const std::vector<uint32_t> & counts){
	if(counts.size() != letters.size()){
		std::stringstream ss;
		ss << "Error in createLikelihood(const std::vector<char> &letters,"
				" const std::vector<uint32_t> & counts)" << std::endl;
		ss << "Size of counts differs from size of letters" << std::endl;
		ss << "Size of counts: " << counts.size() << std::endl;
		ss << "Counts: "; printVector(counts, ", ");
		ss << "Size of letters: " << letters.size() << std::endl;
		ss << "Letters: "; printVector(letters, ", ");
		throw std::runtime_error{ss.str()};
	}
	double countsSum = vectorSum(counts);

  std::multimap<double, std::string, std::less<double>> likelihoods;
  for (const auto &pos : iter::range(letters.size())) {
    likelihoods.emplace(counts[pos]/countsSum, std::string(1, letters[pos]));
  }
  return likelihoods;
}

std::string evenRandStr(uint32_t size, const std::vector<char> &letters,
                        randomGenerator &gen) {
  return randStrMap(size, getEvenLikelihood(letters), gen);
}
VecStr evenRandStrs(uint32_t size, const std::vector<char> &letters,
                    randomGenerator &gen, uint32_t strNum) {
  return randStrsMap(size, getEvenLikelihood(letters), gen, strNum);
}
VecStr evenRandStrsRandLen(uint32_t minLen, uint32_t maxLen,
                           const std::vector<char> &letters,
                           randomGenerator &gen, uint32_t strNum) {
  return randStrsRandLenMap(minLen, maxLen, getEvenLikelihood(letters), gen,
                         strNum);
}

std::string randStr(std::vector<letterCounter> counts, randomGenerator &gen){
	return randStr<letterCounter, char>(counts, gen);
}

}  // sim
}  // bib
