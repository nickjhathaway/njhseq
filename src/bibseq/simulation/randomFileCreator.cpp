/*
 * randomFileCreator.cpp
 *
 *  Created on: Mar 23, 2014
 *      Author: nickhathaway
 */

#include "randomFileCreator.hpp"
#include "bibseq/simulation/randomStrGen.hpp"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
namespace bibseq {

// functions
void randomFileCreator::randomFastq(uint32_t len, uint32_t numOfSeqs,
                                    std::ostream& out, uint32_t offset,
                                    bool processed, uint32_t topAmount) {
  auto seqs = simulation::evenRandStrs(len, alphabet_, gen_, numOfSeqs);
  std::vector<std::vector<uint32_t>> quals;
  quals.reserve(numOfSeqs);
  for (uint32_t i = 0; i < numOfSeqs; ++i) {
    quals.emplace_back(gen_.unifRandVector(qualStart_, qualStop_ + 1, len));
  }
  for (const auto& readNum : iter::range(numOfSeqs)) {
    out << "@Seq." << leftPadNumStr(readNum, numOfSeqs);
    if (processed) {
      out << "_t" << gen_.unifRand<uint32_t>(1, topAmount);
    }
    out << std::endl;
    out << seqs[readNum] << std::endl;
    out << "+" << std::endl;
    out << seqInfo::getFastqString(quals[readNum], offset) << std::endl;
  }
}
void randomFileCreator::randomFasta(uint32_t strLen, uint32_t strNum,
		uint32_t width, const std::vector<char> & alphabetVec,
		randomGenerator & gen, bool processed, uint32_t topAmount,
		std::ostream & out){
  VecStr randomSeqs = simulation::evenRandStrs(strLen, alphabetVec,gen, strNum );
  for(const auto & seqPos : iter::range(randomSeqs.size())){
  	out << ">Seq." << leftPadNumStr<uint32_t>(seqPos, strNum);
  	if(processed){
  		out << "_t" << gen.unifRand<uint32_t>(1, topAmount + 1);
  	}
  	out << std::endl;
  	uint32_t currentPos = 0;
  	while(currentPos < randomSeqs[seqPos].size()){
  		out << randomSeqs[seqPos].substr(currentPos,width) << std::endl;
  		currentPos += width;
  	}
  }
}
} /* namespace bib */
