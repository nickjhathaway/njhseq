/*
 * motif.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: nickhathaway
 */

#include "motif.hpp"

namespace bibseq {


void motifSubUnit::setScoreArray(){
	if(inclusive_){
		// if inclusive give all other letters a score of zero and
		// any in the input as a score of 1 for match
		for(const auto & c : getUpperCaseLetters()){
			score_[c - 'A'] = 0;
		}
		for(const auto & aa : aas_){
			score_[aa - 'A'] = 1;
		}
	}else{
		// if exclusive give all other letters a score of 1 and
		// any in the input as a score of 0 for a non-match
		for(const auto & c : getUpperCaseLetters()){
			score_[c - 'A'] = 1;
		}
		for(const auto & aa : aas_){
			score_[aa - 'A'] = 0;
		}
	}
}
uint32_t motifSubUnit::scoreChar(char c) const{
	return score_[c - 'A'];
}
motifSubUnit motif::processInclusion(uint32_t start, uint32_t stop){
	std::vector<char> include;
	for(const auto & pos : iter::range(start + 1, stop)){
		include.emplace_back(motifOriginal_[pos]);
	}
	return motifSubUnit(include, true);
}
motifSubUnit motif::processExclusion(uint32_t start, uint32_t stop){
	std::vector<char> exclude;
	for(const auto & pos : iter::range(start + 1, stop)){
		exclude.emplace_back(motifOriginal_[pos]);
	}
	return motifSubUnit(exclude, false);
}
void motif::processMotif(){
  auto forwardBrackets = findOccurences(motifOriginal_, "{");
 // printVector(forwardBrackets);
  auto backBrackets = findOccurences(motifOriginal_, "}");
 // printVector(backBrackets);
  std::vector<std::pair<uint32_t, uint32_t>> exclusionPairs;
  for(const auto & pos : iter::range(forwardBrackets.size())){
  	exclusionPairs.emplace_back(std::pair<uint32_t,uint32_t>{forwardBrackets[pos],backBrackets[pos]} );
  }
  auto forwardBrace = findOccurences(motifOriginal_, "[");
 // printVector(forwardBrace);
  auto backBrace = findOccurences(motifOriginal_, "]");
 // printVector(backBrace);
  std::vector<std::pair<uint32_t, uint32_t>> inclusionPairs;
  for(const auto & pos : iter::range(forwardBrace.size())){
  	inclusionPairs.emplace_back(std::pair<uint32_t,uint32_t>{forwardBrace[pos],backBrace[pos]} );
  }
  std::vector<uint32_t> singles(motifOriginal_.size());
  std::map<uint32_t, uint32_t> offSets;
  for(const auto & pos : iter::range(motifOriginal_.size())){
  	offSets[pos] = pos;
  }
  bib::iota<uint32_t>(singles, 0);
  for(const auto & includ : inclusionPairs){
  	removeElements(singles,getRange(includ.first,includ.second));
  	for(const auto & off : iter::range<uint32_t>(includ.first + 1,motifOriginal_.size())){
  		offSets[off] = offSets[off] - (includ.second - includ.first);
  	}
  }
  for(const auto & exclud : exclusionPairs){
  	removeElements(singles,getRange(exclud.first,exclud.second));
  	for(const auto & off : iter::range<uint32_t>(exclud.first + 1,motifOriginal_.size())){
  		offSets[off] = offSets[off] - (exclud.second - exclud.first);
  	}
  }
  /*
  printVector(singles);
  for(const auto & offSet : offSets){
  	std::cout << offSet.first << " : " << offSet.second << std::endl;
  }*/
  //add singles
  for(const auto & sing : singles){
  	motifUnits_[offSets[sing]] = motifSubUnit(motifOriginal_.substr(sing,1) ,true);
  }
  //add inclusions
  for(const auto & includ : inclusionPairs){
  	motifUnits_[offSets[includ.first]] = processInclusion(includ.first,includ.second);
  }
  //add exclusions
  for(const auto & exclud : exclusionPairs){
  	motifUnits_[offSets[exclud.first]] = processExclusion(exclud.first,exclud.second);
  }
  /*
  for(const auto & unit : motifUnits_){
  	std::cout << "pos: " << unit.first << std::endl;
  	std::cout << "\tinclud: " << convertBoolToString(unit.second.inclusive_) <<" "; printVector(unit.second.aas_, " ");
  }*/

  //getRange(0,2);

}
uint32_t motif::scoreMotif(const std::string & possibleMotif){
	if(possibleMotif.size() !=  motifUnits_.size()){
		std::cout << "motif size doesn't equal size of the check" << std::endl;
	}
	uint32_t score = 0;
	for(const auto & cPos : iter::range(possibleMotif.size())){
		score += motifUnits_[cPos].scoreChar(possibleMotif[cPos]);
	}
	return score;
}

uint32_t motif::scoreMotif(const std::string::const_iterator & targetBegin,
		const std::string::const_iterator & targetEnd ){
	if((targetEnd - targetBegin) !=  motifUnits_.size()){
		std::cout << "motif size doesn't equal size of the check" << std::endl;
	}
	uint32_t score = 0;
	auto strIt = targetBegin;
	auto mapIt = motifUnits_.begin();
	for(; strIt!=targetEnd; ++strIt, ++mapIt){
		score += mapIt->second.scoreChar(*strIt);
	}
	return score;
}


bool motif::passMotifParameter(const std::string & possibleMotif, uint32_t scoreCutOff){
	return scoreMotif(possibleMotif) >= scoreCutOff;
}

std::vector<uint32_t> motif::findPositions(const std::string & wholeProtein, uint32_t scoreCutOff){
	uint32_t motifSize = motifUnits_.size();
	uint32_t pos = 0;
	std::vector<uint32_t> positions;
	while(pos + motifSize <= wholeProtein.size()){
		if(passMotifParameter(wholeProtein.substr(pos, motifSize), scoreCutOff)){
			positions.emplace_back(pos);
		}
		++pos;
	}
	return positions;
}

std::vector<uint32_t> motif::findPositionsFull(const std::string & wholeProtein,
		uint32_t allowableErrors){
	uint32_t sum = 0;
	auto predTest = [&sum, &allowableErrors] (const char & a, decltype(*motifUnits_.begin()) mot)
		                   {
		sum += 1 - mot.second.scoreChar(a);
		return sum <= allowableErrors;
		                   };
	uint32_t pos = 0;
	std::vector<uint32_t> positions;
	uint32_t motifSize = motifUnits_.size();
	while(pos + motifUnits_.size() <= wholeProtein.size()){
		sum = 0;
		if(std::equal(wholeProtein.begin() + pos, wholeProtein.begin() + motifSize + pos,
				motifUnits_.begin(), predTest)){
			positions.emplace_back(pos);
		}
		++pos;
	}
	return positions;
}
std::vector<uint32_t> motif::findPositionsFull(const std::string & wholeProtein,
		uint32_t allowableErrors, uint32_t start, uint32_t stop) const{
	uint32_t sum = 0;
	auto predTest = [&sum, &allowableErrors] (const char & a, decltype(*motifUnits_.begin()) mot)
		                   {
		sum += 1 - mot.second.scoreChar(a);
		return sum <= allowableErrors;
		                   };
	uint32_t pos = start;
	std::vector<uint32_t> positions;
	uint32_t motifSize = motifUnits_.size();
	while(pos + motifUnits_.size() <= stop){
		sum = 0;
		if(std::equal(wholeProtein.begin() + pos, wholeProtein.begin() + motifSize + pos,
				motifUnits_.begin(), predTest)){
			positions.emplace_back(pos);
		}
		++pos;
	}
	return positions;
}

} /* namespace bib */
