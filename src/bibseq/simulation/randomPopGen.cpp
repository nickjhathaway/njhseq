/*
 * randomPopGen.cpp
 *
 *  Created on: Mar 23, 2014
 *      Author: nickhathaway
 */

#include "randomPopGen.hpp"

namespace bibseq {

void randomPopGen::runPcr(std::map<std::string, uint32_t>& startingReads,
                          double basalErrorRate, uint32_t rounds) {
  for (uint32_t round = 0; round < rounds; ++round) {
    bib::scopedStopWatch roundTime("round_" + to_string(round),true);
    std::cout << "starting round: " << round << std::endl;
    runOnePcr(startingReads, basalErrorRate);
  }
}
void randomPopGen::runOnePcr(std::map<std::string, uint32_t>& startingReads,
                             double basalErrorRate) {
  std::map<std::string, uint32_t> currentRound;
  for (const auto& reads : startingReads) {
    // now attempt to mutate for the current round of pcr, do twice the ammount
    // due to two strands being duplicated
  	for (uint32_t mut = 0; mut <reads.second * 2 ; ++mut){
    //for (const auto& mut : iter::range(reads.second * 2)) {
      ++currentRound[eProfile_.mutateSeqSameErrorRate(
            reads.first, gen_, eProfile_.alphabet_, basalErrorRate)];
    }
  }
  // now add the newly mutated reads to the current reads
  for (const auto& newReads : currentRound) {
    startingReads[newReads.first] += newReads.second;
  }
}

} /* namespace bib */
