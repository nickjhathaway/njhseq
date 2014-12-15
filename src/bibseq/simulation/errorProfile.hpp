#pragma once
//
//  errorProfile.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/simulation/randomGenerator.hpp"
#include "bibseq/objects/counters/letterCounter.hpp"
namespace bibseq {
namespace simulation {

class errorProfile {

 public:
  // constructor
  /*! \brief Default Construct DNA Alphabet
   *
   */
  errorProfile() : alphabet_({'A', 'C', 'G', 'T'}) { reset(); }
  /*! \brief Construct with Given Alphabet
   *
   */
  errorProfile(const std::vector<char>& alphabet) : alphabet_(alphabet) {
    reset();
  }
  // Members
  std::array<std::array<uint32_t, 26>, 26> counts_;
  std::array<std::array<double, 26>, 26> fractions_;
  std::vector<char> alphabet_;
  std::unordered_map<char, std::multimap<double, char> > probsByLet_;

  // functions
  void reset();
  void setEqualProb();
  // set the fractions
  void setFractions();
  // increase the counts of mismatches
  void increaseCountAmount(const std::string& ref, const std::string& compare,
                           uint32_t amount, const char& ignore = '-');

  void increaseCountAmount(char firstBase, char secondBase, uint32_t amount);
  // mutators based on errors
  char mutate(char firstBase, randomGenerator& gen,
              const std::vector<char>& mutateTo) const;
  void mutateInPlace(char& firstBase, randomGenerator& gen,
                     const std::vector<char>& mutateTo) const;

  void mutateSeqInPlace(std::string& seq, randomGenerator& gen,
                        const std::vector<char>& mutateTo,
                        const std::vector<double>& likelihood);
  std::string mutateSeq(const std::string& seq, randomGenerator& gen,
                        const std::vector<char>& mutateTo,
                        const std::vector<double>& likelihood);

  void mutateSeqInPlaceSameErrorRate(std::string& seq, randomGenerator& gen,
                                     const std::vector<char>& mutateTo,
                                     double errorRate);
  std::string mutateSeqSameErrorRate(const std::string& seq,
                                     randomGenerator& gen,
                                     const std::vector<char>& mutateTo,
                                     double errorRate);
  // printing
  void quickPrintProfile(std::ostream& out);
  void quickPrintCounts(std::ostream& out);
  void quickPrintProbs(std::ostream& out);
};

}  // simulation
}  // bib

#ifndef NOT_HEADER_ONLY
#include "errorProfile.cpp"
#endif
