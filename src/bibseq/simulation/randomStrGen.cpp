//
//  randomStrGen.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/27/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "randomStrGen.hpp"

#include "bibseq/simulation/randomGenerator.hpp"
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
		std::cout << "Error in createLikelihood(const std::vector<char> &letters,"
				" const std::vector<uint32_t> & counts)" << std::endl;
		std::cout << "Size of counts differs from size of letters" << std::endl;
		std::cout << "Size of counts: " << counts.size() << std::endl;
		std::cout << "Counts: "; printVector(counts, ", ");
		std::cout << "Size of letters: " << letters.size() << std::endl;
		std::cout << "Letters: "; printVector(letters, ", ");
		exit(1);
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
