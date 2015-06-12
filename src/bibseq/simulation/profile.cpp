//
//  profile.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/27/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "profile.hpp"

namespace bibseq {
namespace simulation {

void profile::increaseCount(const std::string &ref,
                            const std::string &compare) {
  for (const auto &i : iter::range(ref.length())) {
    increaseCountAmount(ref[i], compare[i], 1);
  }
}

void profile::increaseCountAmount(const std::string &ref,
                                  const std::string &compare, uint32_t amount) {
  for (const auto &i : iter::range(ref.length())) {
    increaseCountAmount(ref[i], compare[i], amount);
  }
}

void profile::increaseCountAmount(char firstBase, char secondBase,
                                  uint32_t amount) {
  counts_[firstBase][secondBase] += amount;
}

void profile::resetTotalFrac(){
	for(const auto & row : iter::range(totals_.size())){
		totals_[row] = 0;
		fractions_[row] = 0;
	}
}

void profile::setProfile(bool ignoreGaps) {
  uint32_t totalSum = 0;
  resetTotalFrac();
  for (auto i : iter::range(counts_.size())) {
    uint32_t currentSum = 0;
    for (auto j : iter::range(counts_[i].size())) {
      if (ignoreGaps && (i == '-' || j == '-')) {
        continue;
      }
      currentSum += counts_[i][j];
      totals_[j] += counts_[i][j] - 1;
    }
    currentSum -= counts_[i].size();
    if (ignoreGaps) {
      ++currentSum;
    }
    totalSum += currentSum;
    if (currentSum != 0) {
      for (auto j : iter::range(counts_[i].size())) {
        profile_[i][j] = (double)counts_[i][j] / currentSum;
      }
    }
  }
  for (auto i : iter::range(totals_.size())) {
    fractions_[i] = (double)totals_[i] / totalSum;
  }
}

void profile::quickPrintProfile(std::ostream &out) {
  std::vector<char> dnaLetters = {'A', 'C', 'G', 'T', '-'};
  for (auto first : dnaLetters) {
    for (auto second : dnaLetters) {
      out << profile_[first][second] << " ";
    }
    out << std::endl;
  }
}

void profile::quickPrintCounts(std::ostream &out) {
  std::vector<char> dnaLetters = {'A', 'C', 'G', 'T', '-'};
  out << "\t";
  printVector(dnaLetters, "\t", out);
  for (auto first : dnaLetters) {
  	std::vector<uint32_t> row;
    for (auto second : dnaLetters) {
      row.emplace_back(counts_[first][second]);
    }
    out << first << "\t";
    printVector(row, "\t", out);
  }
}
/*
char profile::mutate(char firstBase, const bib::randomGenerator &gen,
                          const std::vector<char> &mutateTo) const {
  std::multimap<double, char> probsByChar;
  char ans = '-';
  for (const auto c : mutateTo) {
    probsByChar.emplace(profile_[firstBase - 'A'][c - 'A'], c);
  }
  double sum = 0;
  double rando = gen.unifRand();
  for (const auto &prob : probsByChar) {
    sum += prob.first;
    if (sum > rando) {
      ans = prob.second;
      break;
    }
  }
  return ans;
}

void profile::mutateInPlace(char &firstBase, const randomGenerator &gen,
                                 const std::vector<char> &mutateTo) const {
  std::multimap<double, char> probsByChar;
  for (const auto c : mutateTo) {
    probsByChar.emplace(profile_[firstBase - 'A'][c - 'A'], c);
  }
  double sum = 0;
  double rando = gen.unifRand();
  for (const auto &prob : probsByChar) {
    sum += prob.first;
    if (sum > rando) {
      firstBase = prob.second;
      return;
    }
  }
}*/



void profile::printLogOddsMatrix(const std::vector<char> &letters,
                                 bool ignoreGaps, std::ostream &out) {
	/*
  std::cout << "fractions: " << std::endl;
  for (const auto &let : letters) {
    std::cout << let << ":" << fractions_[let] << std::endl;
  }
  std::cout << "totals: " << std::endl;
  for (const auto &let : letters) {
    std::cout << let << ":" << totals_[let] << std::endl;
  }*/
  //totalCounts_.setFractions();
  setProfile(ignoreGaps);
  out << "\t";
  printVector(letters, "\t", out);
  for (const auto &let : letters) {
    out << let;
    for (const auto &secLet : letters) {
      out << "\t" << std::round( 10*( std::log10((profile_[let][secLet]) /
                                (fractions_[let] * fractions_[secLet]))));
    }
    out << std::endl;
  }
}
}  // sim
}  // bib
