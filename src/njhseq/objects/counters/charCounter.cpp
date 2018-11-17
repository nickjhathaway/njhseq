/*
 * charCounter.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: nickhathaway
 */
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "charCounter.hpp"

namespace njhseq {

Json::Value charCounter::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson("njhseq::charCounter");
	auto & chars = ret["chars_"];
	for (auto let : alphabet_) {
		chars[std::string(1, let)] = chars_[let];
	}
	auto & fracs = ret["fractions_"];
	for (auto let : alphabet_) {
		fracs[std::string(1, let)] = fractions_[let];
	}
	auto & quals = ret["qualities_"];
	for (auto let : alphabet_) {
		quals[std::string(1, let)] = qualities_[let];
	}

	auto & allQuals = ret["allQualities_"];
	for (auto let : alphabet_) {
		allQuals[std::string(1, let)] = njh::json::toJson(allQualities_[let]);
	}
	ret["allowNewChars_"] = njh::json::toJson(allowNewChars_);
	ret["alphabet_"] = njh::json::toJson(alphabet_);
	ret["originalAlphabet_"] = njh::json::toJson(originalAlphabet_);
	double GCContent = chars_['G'] + chars_['C'] + chars_['g'] + chars_['c'];
	ret["gcContent"] = njh::json::toJson(GCContent / getTotalCount());
	ret["totalCount"] = njh::json::toJson(getTotalCount());

	return ret;
}

charCounter::charCounter() :
		charCounter(std::vector<char> { 'A', 'C', 'G', 'T', '-' }) {

}

charCounter::charCounter(const std::vector<char>& alphabet) :
		alphabet_(alphabet), originalAlphabet_(alphabet_) {
	reset();
}

charCounter::charCounter(const std::string & str) {
	increaseCountByString(str);
	resetAlphabet(false);
	setFractions();
}

charCounter::charCounter(const std::string & str,
		const std::vector<char>& alphabet) :
		alphabet_(alphabet), originalAlphabet_(alphabet) {
	increaseCountByString(str);
	resetAlphabet(true);
	setFractions();
}

void charCounter::increasePortion(const std::string & str, uint64_t len,
		double cnt) {
	increasePortion(str.begin(), str.begin() + len, cnt);
}

///charArray
void charCounter::getBest(char &letter) const {
	uint32_t bestCount = 0;
	uint64_t bestQualSum = 0;
	for (const auto & c : alphabet_) {
		if (chars_[c] > bestCount) {
			letter = c;
			bestCount = chars_[c];
			bestQualSum = qualities_[c];
		}else if (chars_[c] == bestCount && qualities_[c] > bestQualSum){
			letter = c;
			bestCount = chars_[c];
			bestQualSum = qualities_[c];
		}
	}
}
char charCounter::outputBestLetter()const {
	char bestBase = ' ';
	getBest(bestBase);
	return bestBase;
}

void charCounter::getBest(char &letter, uint32_t &quality) const {
	getBest(letter);
	auto totalCount = getTotalCount();
	if (totalCount < 5) {
		//std::cout << "bestLetter: " << letter << std::endl;
		//std::cout << "qualities_[letter]: " << qualities_[letter] << std::endl;
		//std::cout << "letCount  : " << chars_[letter] << std::endl;
		//std::cout << "totalCount: " << totalCount << std::endl;
		//not casting to double below for division because it is just being cast back anyways
		quality = qualities_[letter] / chars_[letter];
		if(quality > totalCount - chars_[letter]){
			quality-=(totalCount - chars_[letter]);
		}
	} else {
		//quality = qualities_[letter] / totalCount;
		quality = qualities_[letter] /  chars_[letter];
	}
}

//for finding insertions
void charCounter::getBest(char & base, uint32_t &quality, uint32_t size) const {
	uint32_t totalCount = getTotalCount();
	//due to rounding sometimes the total count is greater than size
	//if this happens, we get overflow and size - totalCount becomes a gigantic number
	uint32_t nonInsertingSize = totalCount > size ? 0 : size - totalCount;
	if(totalCount > nonInsertingSize){
		//find best base, if it's count is better than the non inserting reads then output that
		char bestBase = outputBestLetter();
		if(chars_[bestBase] > nonInsertingSize){
			if(size < 5){
				quality = qualities_[bestBase] / chars_[bestBase];
				//due to rounding sometimes the total count is greater than size
				//if this happens, we get overflow and size - chars_[bestBase] becomes a gigantic number
				if(chars_[bestBase] < size && quality > size - chars_[bestBase]){
					quality-=(size - chars_[bestBase]);
				}
			}else{
				quality = qualities_[bestBase] / size;
			}
			base = bestBase;
		}
	}
}

void charCounter::outPutInfo(std::ostream &out, bool ifQualities) const {
	VecStr header;
	if (ifQualities) {
		header = {"letter", "count", "qualities","fraction"};
	} else {
		header = {"letter", "count","fraction"};
	}

	double totalCount = getTotalCount();
	std::vector<VecStr> outConent;
	if (!ifQualities) {
		for (const auto &let : alphabet_) {
			outConent.emplace_back(
					toVecStr(let, chars_[let], chars_[let] / totalCount));
		}
	} else {
		for (const auto &let : alphabet_) {
			outConent.emplace_back(
					toVecStr(let, chars_[let], qualities_[let],
							chars_[let] / totalCount));
		}
	}
	printTableOrganized(outConent, header, out);
}

void charCounter::outPutACGTInfo(std::ostream &out) const {
	std::vector<uint32_t> counts;
	for (const auto & c : alphabet_) {
		counts.emplace_back(chars_[c]);
	}
	printVector(counts, " ", out);
}

void charCounter::outPutACGTFractionInfo(std::ostream &out) {
	setFractions();
	std::vector<double> counts;
	for (const auto & c : alphabet_) {
		counts.emplace_back(fractions_[c]);
	}
	printVector(counts, " ", out);
}

void charCounter::calcGcContent() {
	double GCContent = chars_['G'] + chars_['C'] + chars_['g'] + chars_['c'];
	gcContent_ = GCContent / getTotalCount();
}

double charCounter::computeEntrophy() {
	setFractions();
	double sum = 0;
	for (const auto &c : alphabet_) {
		if (0 != fractions_[c]) {
			if (1 == fractions_[c]) {
				return 0;
			}
			sum += fractions_[c] * std::log2(fractions_[c]);
		}
	}
	return (-1 * sum);
}

char charCounter::getDegenativeBase() const {
	if (chars_['A'] == 0 && chars_['C'] == 0 && chars_['G'] == 0
			&& chars_['T'] > 0) {
		return 'T';
	} else if (chars_['A'] == 0 && chars_['C'] == 0 && chars_['G'] > 0
			&& chars_['T'] == 0) {
		return 'G';
	} else if (chars_['A'] == 0 && chars_['C'] > 0 && chars_['G'] == 0
			&& chars_['T'] == 0) {
		return 'C';
	} else if (chars_['A'] > 0 && chars_['C'] == 0 && chars_['G'] == 0
			&& chars_['T'] == 0) {
		return 'A';
	} else if (chars_['A'] == 0 && chars_['C'] == 0 && chars_['G'] > 0
			&& chars_['T'] > 0) {
		return 'K';
	} else if (chars_['A'] == 0 && chars_['C'] > 0 && chars_['G'] == 0
			&& chars_['T'] > 0) {
		return 'Y';
	} else if (chars_['A'] > 0 && chars_['C'] == 0 && chars_['G'] == 0
			&& chars_['T'] > 0) {
		return 'W';
	} else if (chars_['A'] == 0 && chars_['C'] > 0 && chars_['G'] > 0
			&& chars_['T'] == 0) {
		return 'S';
	} else if (chars_['A'] > 0 && chars_['C'] == 0 && chars_['G'] > 0
			&& chars_['T'] == 0) {
		return 'R';
	} else if (chars_['A'] > 0 && chars_['C'] > 0 && chars_['G'] == 0
			&& chars_['T'] == 0) {
		return 'M';
	} else if (chars_['A'] == 0 && chars_['C'] > 0 && chars_['G'] > 0
			&& chars_['T'] > 0) {
		return 'B';
	} else if (chars_['A'] > 0 && chars_['C'] == 0 && chars_['G'] > 0
			&& chars_['T'] > 0) {
		return 'D';
	} else if (chars_['A'] > 0 && chars_['C'] > 0 && chars_['G'] == 0
			&& chars_['T'] > 0) {
		return 'H';
	} else if (chars_['A'] > 0 && chars_['C'] > 0 && chars_['G'] > 0
			&& chars_['T'] == 0) {
		return 'V';
	} else {
		return 'N';
	}
}
int charCounter::getGcDifference() {
	return chars_['G'] - chars_['C'];
}

void charCounter::reset() {
	for (const auto & pos : iter::range(chars_.size())) {
		chars_[pos] = 0;
		fractions_[pos] = 0;
		qualities_[pos] = 0;
		allQualities_[pos].clear();
	}
}

void charCounter::increaseCountOfBase(const char &base) {
	chars_[base] += 1;
}
void charCounter::increaseCountOfBase(const char &base, double cnt) {
	chars_[base] += cnt;
}
void charCounter::increaseCountByString(const std::string &seq) {
	for (const auto & c : seq) {
		chars_[c] += 1;
	}
}
void charCounter::increaseCountByString(const std::string &seq, double cnt) {
	for (const auto & c : seq) {
		chars_[c] += cnt;
	}
}

void charCounter::increaseCountOfBaseQual(const char &base, uint32_t qual) {
	chars_[base] += 1;
	qualities_[base] += qual;
	allQualities_[base].emplace_back(qual);
}
void charCounter::increaseCountOfBaseQual(const char &base, uint32_t qual,
		double cnt) {
	chars_[base] += cnt;
	qualities_[base] += qual * cnt;
	addOtherVec(allQualities_[base], std::vector<uint32_t>(cnt, qual));
}

void charCounter::increaseCountByStringQual(const std::string &seq,
		const std::vector<uint32_t> & qualities) {
	for (const auto & pos : iter::range(seq.size())) {
		chars_[seq[pos]] += 1;
		qualities_[seq[pos]] += qualities[pos];
		allQualities_[seq[pos]].emplace_back(qualities[pos]);
	}
}
void charCounter::increaseCountByStringQual(const std::string &seq,
		const std::vector<uint32_t> & qualities, double cnt) {
	for (const auto & pos : iter::range(seq.size())) {
		chars_[seq[pos]] += cnt;
		qualities_[seq[pos]] += qualities[pos] * cnt;
		addOtherVec(allQualities_[seq[pos]],
				std::vector<uint32_t>(cnt, qualities[pos]));
	}
}

void charCounter::setFractions() {
	setFractions(alphabet_);
}

double charCounter::getFracDifference(const charCounter & otherCounter,
		const std::vector<char> & alph) const {
	double sum = 0;
	for (const auto & let : alph) {
		sum += std::abs(fractions_[let] - otherCounter.fractions_[let]);
	}
	return sum;
}

void charCounter::setFractions(const std::vector<char>& alphabet) {
	uint32_t total = 0;
	for (const auto & c : alphabet) {
		total += chars_[c];
	}
	if (total != 0) {
		for (const auto & c : alphabet) {
			fractions_[c] = chars_[c] / static_cast<double>(total);
		}
	} else {
		for (const auto & c : alphabet) {
			fractions_[c] = 0;
		}
	}
}

uint32_t charCounter::getTotalCount() const {
	uint32_t total = 0;
	for (const auto & c : alphabet_) {
		total += chars_[c];
	}
	return total;
}

std::multimap<double, char, std::less<double>> charCounter::createLikelihoodMaps(
		bool setFractionFirst) {
	if (setFractionFirst) {
		setFractions();
	}
	std::multimap<double, char, std::less<double>> likelihoods;
	for (const auto &c : alphabet_) {
		likelihoods.emplace(fractions_[c], c);
	}
	return likelihoods;
}

void charCounter::resetAlphabet(bool keepOld) {
	std::vector<char> present;
	for (auto i : iter::range(chars_.size())) {
		if (chars_[i] > 0) {
			present.emplace_back(i);
		}
	}
	if (!keepOld) {
		alphabet_.clear();
	}
	for (const auto & c : present) {
		if (!njh::in(c, alphabet_)) {
			alphabet_.emplace_back(c);
		}
	}
	njh::sort(alphabet_);
}

void charCounter::addOtherCounts(const charCounter & otherCounter,
		bool setFractionsAfter) {
	for (const auto & pos : iter::range(otherCounter.chars_.size())) {
		chars_[pos] += otherCounter.chars_[pos];
	}
	if (setFractionsAfter) {
		setFractions();
	}
}

void charCounter::increaseRates(substituteMatrix & mat, char refBase) const{
	for(char let : alphabet_){
		mat.mat_[refBase][let] += chars_[let];
	}
}

} /* namespace njhseq */
