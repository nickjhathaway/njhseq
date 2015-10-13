/*
 * charCounter.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: nickhathaway
 */

#include "charCounter.hpp"

namespace bibseq {

charCounterArray::charCounterArray() :
		charCounterArray(std::vector<char> { 'A', 'C', 'G', 'T', '-' }) {

}

charCounterArray::charCounterArray(const std::vector<char>& alphabet) :
		alphabet_(alphabet), originalAlphabet_(alphabet_) {
	reset();
}

charCounterArray::charCounterArray(const std::string & str) {
	increaseCountByString(str);
	resetAlphabet(false);
	setFractions();
}

charCounterArray::charCounterArray(const std::string & str,
		const std::vector<char>& alphabet): alphabet_(alphabet), originalAlphabet_(alphabet) {
	increaseCountByString(str);
	resetAlphabet(true);
	setFractions();
}

void charCounterArray::increasePortion(const std::string & str, uint64_t len, double cnt){
	increasePortion(str.begin(), str.begin() + len, cnt);
}


///charArray
void charCounterArray::getBest(char &letter) const {
  uint32_t bestCount = 0;
  for (const auto & c : alphabet_) {
    if (chars_[c] > bestCount) {
      letter = c;
      bestCount = chars_[c];
    }
  }
}
char charCounterArray::outputBestLetter() {
  uint32_t bestCount = 0;
  char bestBase = ' ';
  for (const auto & c : alphabet_ ) {
    if (chars_[c] > bestCount) {
      bestBase = c;
      bestCount = chars_[c];
    }
  }
  return bestBase;
}

void charCounterArray::getBest(char &letter, uint32_t &quality)const {
  uint32_t bestCount = 0;
  char bestBase = ' ';
  uint32_t bestQuality = 0;

  for (const auto & c : alphabet_ ) {
    if (chars_[c] > bestCount) {
      bestBase = c;
      bestCount = chars_[c];
      bestQuality = qualities_[c];
    }
  }
  letter = bestBase;
  quality = bestQuality;
}

//for finding for insertions
void charCounterArray::getBest(char & base, uint32_t &quality, uint32_t size)const {
  uint32_t totalCount = getTotalCount();
  uint32_t nonInsertingSize = size - totalCount;
  uint32_t bestCount = nonInsertingSize;
  //find best count
  for (const auto & c : alphabet_ ) {
    if (chars_[c] > bestCount) {
    	base = c;
    	quality = qualities_[c]/size;
      bestCount = chars_[c];
    }else if(chars_[c] == bestCount){
    	//if count is the same get the let with the best qualities, if same qualities
    	//whatever one was found first
    	//if the bestCount is the non inserting, opt to not putting an insert
    	if(bestCount != nonInsertingSize){
      	if(qualities_[c]/size > quality){
      		base = c;
      		quality = qualities_[c];
      	}
    	}
    }
  }
}

void charCounterArray::outPutInfo(std::ostream &out, bool ifQualities) const {
	VecStr header;
	if(ifQualities){
		header = { "letter", "count", "qualities","fraction"};
	}else{
		header = { "letter", "count","fraction"};
	}

  double totalCount = getTotalCount();
  std::vector<VecStr> outConent;
  if (!ifQualities) {
    for (const auto &let : alphabet_) {
    	outConent.emplace_back(VecStr{to_string(let), to_string(chars_[let]), to_string(chars_[let] / totalCount) });
    }
  } else {
  	for (const auto &let : alphabet_) {
  		outConent.emplace_back(VecStr{to_string(let), to_string(chars_[let]),to_string(qualities_[let]), to_string(chars_[let] / totalCount) });
    }
  }
  printTableOrganized(outConent, header, out);
}

void charCounterArray::outPutACGTInfo(std::ostream &out) const {
	std::vector<uint32_t> counts;
	for(const auto & c : alphabet_){
		counts.emplace_back(chars_[c]);
	}
	printVector(counts," ", out);
}

void charCounterArray::outPutACGTFractionInfo(std::ostream &out) {
  setFractions();
	std::vector<double> counts;
	for(const auto & c : alphabet_){
		counts.emplace_back(fractions_[c]);
	}
	printVector(counts," ", out);
}

void charCounterArray::calcGcContent() {
  double GCContent = chars_['G'] + chars_['C'];
  gcContent = GCContent / getTotalCount();
}
double charCounterArray::computeEntrophy() {
  setFractions();
  double sum = 0;
  for (const auto &c : alphabet_) {
    if (0 != fractions_[c]) {
      if (1 == fractions_[c]) {
        return 0;
      }
      sum += fractions_[c]* std::log2(fractions_[c]);
    }
  }
  return (-1 * sum);
}
char charCounterArray::getDegenativeBase() const {
  if (chars_['A'] == 0 && chars_['C'] == 0 && chars_['G'] == 0 &&
      chars_['T'] > 0) {
    return 'T';
  } else if (chars_['A'] == 0 && chars_['C'] == 0 &&
             chars_['G'] > 0 && chars_['T'] == 0) {
    return 'G';
  } else if (chars_['A'] == 0 && chars_['C'] > 0 &&
             chars_['G'] == 0 && chars_['T'] == 0) {
    return 'C';
  } else if (chars_['A'] > 0 && chars_['C'] == 0 &&
             chars_['G'] == 0 && chars_['T'] == 0) {
    return 'A';
  } else if (chars_['A'] == 0 && chars_['C'] == 0 &&
             chars_['G'] > 0 && chars_['T'] > 0) {
    return 'K';
  } else if (chars_['A'] == 0 && chars_['C'] > 0 &&
             chars_['G'] == 0 && chars_['T'] > 0) {
    return 'Y';
  } else if (chars_['A'] > 0 && chars_['C'] == 0 &&
             chars_['G'] == 0 && chars_['T'] > 0) {
    return 'W';
  } else if (chars_['A'] == 0 && chars_['C'] > 0 &&
             chars_['G'] > 0 && chars_['T'] == 0) {
    return 'S';
  } else if (chars_['A'] > 0 && chars_['C'] == 0 &&
             chars_['G'] > 0 && chars_['T'] == 0) {
    return 'R';
  } else if (chars_['A'] > 0 && chars_['C'] > 0 &&
             chars_['G'] == 0 && chars_['T'] == 0) {
    return 'M';
  } else if (chars_['A'] == 0 && chars_['C'] > 0 &&
             chars_['G'] > 0 && chars_['T'] > 0) {
    return 'B';
  } else if (chars_['A'] > 0 && chars_['C'] == 0 &&
             chars_['G'] > 0 && chars_['T'] > 0) {
    return 'D';
  } else if (chars_['A'] > 0 && chars_['C'] > 0 &&
             chars_['G'] == 0 && chars_['T'] > 0) {
    return 'H';
  } else if (chars_['A'] > 0 && chars_['C'] > 0 &&
             chars_['G'] > 0 && chars_['T'] == 0) {
    return 'V';
  } else {
    return 'N';
  }
}
int charCounterArray::getGcDifference() { return chars_['G'] - chars_['C']; }

void charCounterArray::printDescription(std::ostream &out, bool deep) const {
  out << "charCounterArray{" << std::endl;
 /* out << "letters:" << std::endl;
  out << "\tstd::map<std::string, int>" << std::endl;
  for (const auto &let : letters) {
    out << "\t" << let.first << ":" << let.second << std::endl;
  }
  out << "fractions:" << std::endl;
  out << "\tstd::map<std::string, double>" << std::endl;
  for (const auto &frac : fractions) {
    out << "\t" << frac.first << ":" << frac.second << std::endl;
  }
  out << "qualities:" << std::endl;
  out << "\tstd::map<std::string, uint32_t>" << std::endl;
  for (const auto &qual : qualities) {
    out << "\t" << qual.first << ":" << qual.second << std::endl;
  }
  out << "allQualities:" << std::endl;
  out << "\tstd::map<std::string, std::vector<uint32_t>>" << std::endl;
  for (const auto &qual : allQualities) {
    out << "\t" << qual.first << ":" << qual.second << std::endl;
  }
  out << "gcContent:" << gcContent << std::endl << "}" << std::endl;*/
}
void charCounterArray::reset(){
	for(const auto & pos : iter::range(chars_.size())){
		chars_[pos] = 0;
		fractions_[pos] = 0;
		qualities_[pos] = 0;
		allQualities_[pos].clear();
	}
}

void charCounterArray::increaseCountOfBase(const char &base){
	chars_[base] += 1;
}
void charCounterArray::increaseCountOfBase(const char &base, double cnt){
	chars_[base] += cnt;
}
void charCounterArray::increaseCountByString(const std::string &seq){
	for(const auto & c : seq){
		chars_[c] += 1;
	}
}
void charCounterArray::increaseCountByString(const std::string &seq, double cnt){
	for(const auto & c : seq){
		chars_[c] += cnt;
	}
}

void charCounterArray::increaseCountOfBaseQual(const char &base, uint32_t qual){
	chars_[base] += 1;
	qualities_[base] += qual;
	allQualities_[base].emplace_back(qual);
}
void charCounterArray::increaseCountOfBaseQual(const char &base, uint32_t qual, double cnt){
	chars_[base] += cnt;
	qualities_[base] += qual*cnt;
	addOtherVec(allQualities_[base], std::vector<uint32_t>(cnt, qual));
}
/**
 *
 * @param seq
 * @param qualities
 * @todo check for small size of seq and qualitites
 */
void charCounterArray::increaseCountByStringQual(const std::string &seq, const std::vector<uint32_t> & qualities){
	for(const auto & pos : iter::range(seq.size())){
		chars_[seq[pos]] += 1;
		qualities_[seq[pos]] += qualities[pos];
		allQualities_[seq[pos]].emplace_back(qualities[pos]);
	}
}
void charCounterArray::increaseCountByStringQual(const std::string &seq, const std::vector<uint32_t> & qualities, double cnt){
	for(const auto & pos : iter::range(seq.size())){
		chars_[seq[pos]] += cnt;
		qualities_[seq[pos]] += qualities[pos] *cnt;
		addOtherVec(allQualities_[seq[pos]], std::vector<uint32_t>(cnt, qualities[pos]));
	}
}

void charCounterArray::setFractions(){
	setFractions(alphabet_);
	/*
	uint32_t total = 0;
	for(const auto & c : alphabet_){
		total += chars_[c];
	}
	if(total != 0){
		for(const auto & c : alphabet_ ){
			fractions_[c] = chars_[c]/static_cast<double>(total);
		}
	}*/
}
double charCounterArray::getFracDifference(const charCounterArray & otherCounter, const std::vector<char> & alph)const{
	double sum = 0;
	for(const auto & let : alph){
		sum += std::abs(fractions_[let] - otherCounter.fractions_[let]);
	}
	return sum;
}
void charCounterArray::setFractions(const std::vector<char>& alphabet){
	uint32_t total = 0;
	for(const auto & c : alphabet){
		total += chars_[c];
	}
	if(total != 0){
		for(const auto & c : alphabet ){
			fractions_[c] = chars_[c]/static_cast<double>(total);
		}
	}
}
uint32_t charCounterArray::getTotalCount() const {
	uint32_t total = 0;
	for(const auto & c : alphabet_){
		total += chars_[c];
	}
	return total;
}

std::multimap<double, char, std::less<double>>
charCounterArray::createLikelihoodMaps(bool setFractionFirst) {
  if (setFractionFirst) {
    setFractions();
  }
  std::multimap<double, char, std::less<double>> likelihoods;
  for (const auto &c : alphabet_) {
    likelihoods.emplace(fractions_[c], c);
  }
  return likelihoods;
}

void charCounterArray::resetAlphabet(bool keepOld){
	std::vector<char> present;
	for(auto i : iter::range(chars_.size())){
		if(chars_[i] > 0){
			present.emplace_back(i);
		}
	}
	if(!keepOld){
		alphabet_.clear();
	}
	for(const auto & c : present){
		if(!in(c, alphabet_)){
			alphabet_.emplace_back(c);
		}
	}
	sort(alphabet_);
}

void charCounterArray::addOtherCounts(const charCounterArray & otherCounter, bool setFractionsAfter){
	for(const auto & pos : iter::range(otherCounter.chars_.size())){
		chars_[pos]+=otherCounter.chars_[pos];
	}
	if(setFractionsAfter){
		setFractions();
	}
}

} /* namespace bib */
