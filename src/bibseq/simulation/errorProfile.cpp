//
//  errorProfile.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "errorProfile.hpp"

namespace bibseq {
namespace simulation {

void errorProfile::increaseCountAmount(const std::string &ref,
                                       const std::string &compare,
                                       uint32_t amount, const char &ignore) {
  auto mis = std::make_pair(ref.begin(), compare.begin());
  while (mis.first != ref.end()) {
    mis = std::mismatch(mis.first, ref.end(), mis.second);
    if (mis.first != ref.end()) {
      if (*mis.first != ignore && *mis.second != ignore) {
        increaseCountAmount(*mis.first, *mis.second, amount);
      }
      ++mis.first;
      ++mis.second;
    }
  }
}

void errorProfile::increaseCountAmount(char firstBase, char secondBase,
                                       uint32_t amount) {
  counts_[firstBase - 'A'][secondBase - 'A'] += amount;
}

void errorProfile::setFractions() {
  // sum across letter counts to get full amount of changes and then set the
  // fractions
  for (auto i : iter::range(counts_.size())) {
    uint32_t currentSum = 0;
    for (auto j : iter::range(counts_[i].size())) {
      currentSum += counts_[i][j];
    }
    if (currentSum != 0) {
      for (auto j : iter::range(counts_[i].size())) {
        fractions_[i][j] = (double)counts_[i][j] / currentSum;
      }
    }
  }
  // now calc the prob of change to each letter in set alphabet
  for (const auto &first : alphabet_) {
    std::multimap<double, char> currentProbsByChar;
    for (const auto &second : alphabet_) {
      currentProbsByChar.emplace(fractions_[first - 'A'][second - 'A'], second);
    }
    probsByLet_[first] = currentProbsByChar;
  }
}
void errorProfile::reset() {
  for (auto i : iter::range(fractions_.size())) {
    for (auto j : iter::range(fractions_[i].size())) {
      fractions_[i][j] = 0;
      counts_[i][j] = 0;
    }
  }
  probsByLet_.clear();
}

// mutate and return a char
char errorProfile::mutate(char firstBase, bibseq::randomGenerator &gen,
                          const std::vector<char> &mutateTo) const {
  mutateInPlace(firstBase, gen, mutateTo);
  return firstBase;
}
// mutate a char to another char
void errorProfile::mutateInPlace(char &firstBase, randomGenerator &gen,
                                 const std::vector<char> &mutateTo) const {
  double sum = 0;
  double rando = gen.unifRand();
  for (const auto &prob : probsByLet_.at(firstBase)) {
    sum += prob.first;
    if (sum > rando) {
      firstBase = prob.second;
      return;
    }
  }
}
// mutate the given sequence
void errorProfile::mutateSeqInPlace(std::string &seq, randomGenerator &gen,
                                    const std::vector<char> &mutateTo,
                                    const std::vector<double> &likelihood) {
  std::vector<double> rands = gen.unifRandVector(likelihood.size());
  for (const auto &i : iter::range(rands.size())) {
    if (rands[i] <= likelihood[i]) {
      mutateInPlace(seq[i], gen, mutateTo);
    }
  }
}
// mutate and return a new sequence
std::string errorProfile::mutateSeq(const std::string &seq,
                                    randomGenerator &gen,
                                    const std::vector<char> &mutateTo,
                                    const std::vector<double> &likelihood) {
  std::string ans = seq;
  mutateSeqInPlace(ans, gen, mutateTo, likelihood);
  return ans;
}
void errorProfile::mutateSeqInPlaceSameErrorRate(
    std::string &seq, randomGenerator &gen, const std::vector<char> &mutateTo,
    double errorRate) {
  std::vector<double> rands = gen.unifRandVector(seq.size());
  for (const auto &i : iter::range(rands.size())) {
    if (rands[i] <= errorRate) {
      mutateInPlace(seq[i], gen, mutateTo);
    }
  }
}
std::string errorProfile::mutateSeqSameErrorRate(
    const std::string &seq, randomGenerator &gen,
    const std::vector<char> &mutateTo, double errorRate) {
  std::string ans = seq;
  mutateSeqInPlaceSameErrorRate(ans, gen, mutateTo, errorRate);
  return ans;
}
void errorProfile::quickPrintProfile(std::ostream &out) {
  std::vector<VecStr> outTable;
  VecStr header = catenateVectors({"let"}, numVecToVecStr(alphabet_));
  for (auto first : alphabet_) {
    VecStr currentRow;
    currentRow.emplace_back(to_string(first));
    for (auto second : alphabet_) {
      currentRow.emplace_back(to_string(fractions_[first - 'A'][second - 'A']));
    }
    outTable.emplace_back(currentRow);
  }
  printTableOrganized(outTable, header, out);
}

void errorProfile::quickPrintCounts(std::ostream &out) {
  std::vector<VecStr> outTable;
  VecStr header = catenateVectors({"let"}, numVecToVecStr(alphabet_));
  for (auto first : alphabet_) {
    VecStr currentRow;
    currentRow.emplace_back(to_string(first));
    for (auto second : alphabet_) {
      currentRow.emplace_back(to_string(counts_[first - 'A'][second - 'A']));
    }
    outTable.emplace_back(currentRow);
  }
  printTableOrganized(outTable, header, out);
}
void errorProfile::quickPrintProbs(std::ostream &out) {
  setFractions();

  std::vector<VecStr> outTable;
  VecStr header = {"let", "secondLet", "prob"};

  for (const auto &let : alphabet_) {
    for (const auto &probs : probsByLet_[let]) {
      outTable.emplace_back(VecStr{to_string(let), to_string(probs.second),
                                   to_string(probs.first)});
      // out << let <<"\t" << probs.first << "\t" << probs.second << std::endl;
    }
  }
  printTableOrganized(outTable, header, out);
}

void errorProfile::setEqualProb() {
  reset();
  for (const auto &let : alphabet_) {
    for (const auto &letSecond : alphabet_) {
      if (let == letSecond) {
        continue;
      }
      counts_[let][letSecond] = 1;
    }
  }
  setFractions();
}
}  // sim
}  // bib
