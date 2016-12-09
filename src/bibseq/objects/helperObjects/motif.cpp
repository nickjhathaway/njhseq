/*
 * motif.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: nickhathaway
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "motif.hpp"

namespace bibseq {

motif::motifSubUnit::motifSubUnit() :
		aas_(getUpperCaseLetters()), inclusive_(true) {
	setScoreArray();
}
motif::motifSubUnit::motifSubUnit(const std::vector<char> & aas, bool inclusive) :
		aas_(aas), inclusive_(inclusive) {
	setScoreArray();
}
motif::motifSubUnit::motifSubUnit(const std::string & motifSub, bool inclusive) :
		inclusive_(inclusive) {
	for (const auto & c : motifSub) {
		aas_.emplace_back(c);
	}
	setScoreArray();
}

void motif::motifSubUnit::setScoreArray() {
	if (inclusive_) {
		// if inclusive give all other letters a score of zero and
		// any in the input as a score of 1 for match
		for (const auto & c : getUpperCaseLetters()) {
			score_[c - 'A'] = 0;
		}
		for (const auto & aa : aas_) {
			score_[aa - 'A'] = 1;
		}
	} else {
		// if exclusive give all other letters a score of 1 and
		// any in the input as a score of 0 for a non-match
		for (const auto & c : getUpperCaseLetters()) {
			score_[c - 'A'] = 1;
		}
		for (const auto & aa : aas_) {
			score_[aa - 'A'] = 0;
		}
	}
}
uint32_t motif::motifSubUnit::scoreChar(char c) const {
	return score_[c - 'A'];
}
motif::motifSubUnit motif::processInclusion(uint32_t start, uint32_t stop) {
	std::vector<char> include;
	for (const auto & pos : iter::range(start + 1, stop)) {
		include.emplace_back(motifOriginal_[pos]);
	}
	return motifSubUnit(include, true);
}
motif::motifSubUnit motif::processExclusion(uint32_t start, uint32_t stop) {
	std::vector<char> exclude;
	for (const auto & pos : iter::range(start + 1, stop)) {
		exclude.emplace_back(motifOriginal_[pos]);
	}
	return motifSubUnit(exclude, false);
}
void motif::processMotif() {
	auto forwardBrackets = findOccurences(motifOriginal_, "{");
	// printVector(forwardBrackets);
	auto backBrackets = findOccurences(motifOriginal_, "}");
	// printVector(backBrackets);
	std::vector<std::pair<uint32_t, uint32_t>> exclusionPairs;
	for (const auto & pos : iter::range(forwardBrackets.size())) {
		exclusionPairs.emplace_back(std::pair<uint32_t, uint32_t> {
				forwardBrackets[pos], backBrackets[pos] });
	}
	auto forwardBrace = findOccurences(motifOriginal_, "[");
	// printVector(forwardBrace);
	auto backBrace = findOccurences(motifOriginal_, "]");
	// printVector(backBrace);
	std::vector<std::pair<uint32_t, uint32_t>> inclusionPairs;
	for (const auto & pos : iter::range(forwardBrace.size())) {
		inclusionPairs.emplace_back(std::pair<uint32_t, uint32_t> {
				forwardBrace[pos], backBrace[pos] });
	}
	std::vector<uint32_t> singles(motifOriginal_.size(), 0);
	std::map<uint32_t, uint32_t> offSets;
	for (const auto & pos : iter::range(motifOriginal_.size())) {
		offSets[pos] = pos;
	}
	bib::iota<uint32_t>(singles, 0);
	for (const auto & includ : inclusionPairs) {
		removeElements(singles, getRange(includ.first, includ.second));
		for (const auto & off : iter::range<uint32_t>(includ.first + 1,
				motifOriginal_.size())) {
			offSets[off] = offSets[off] - (includ.second - includ.first);
		}
	}
	for (const auto & exclud : exclusionPairs) {
		removeElements(singles, getRange(exclud.first, exclud.second));
		for (const auto & off : iter::range<uint32_t>(exclud.first + 1,
				motifOriginal_.size())) {
			offSets[off] = offSets[off] - (exclud.second - exclud.first);
		}
	}
	/*
	 printVector(singles);
	 for(const auto & offSet : offSets){
	 std::cout << offSet.first << " : " << offSet.second << std::endl;
	 }*/
	//add singles
	for (const auto & sing : singles) {
		motifUnits_[offSets[sing]] = motifSubUnit(motifOriginal_.substr(sing, 1),
				true);
	}
	//add inclusions
	for (const auto & includ : inclusionPairs) {
		motifUnits_[offSets[includ.first]] = processInclusion(includ.first,
				includ.second);
	}
	//add exclusions
	for (const auto & exclud : exclusionPairs) {
		motifUnits_[offSets[exclud.first]] = processExclusion(exclud.first,
				exclud.second);
	}
	/*
	 for(const auto & unit : motifUnits_){
	 std::cout << "pos: " << unit.first << std::endl;
	 std::cout << "\tinclud: " << convertBoolToString(unit.second.inclusive_) <<" "; printVector(unit.second.aas_, " ");
	 }*/

	//getRange(0,2);
}

motif::motif(const std::string & inMotif) :
		motifOriginal_(inMotif) {
	processMotif();
}

uint32_t motif::scoreMotif(const std::string & possibleMotif) const {
	if (possibleMotif.size() != motifUnits_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": motif size doesn't equal size of the check"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	uint32_t score = 0;
	for (const auto & cPos : iter::range(possibleMotif.size())) {
		score += motifUnits_.at(cPos).scoreChar(possibleMotif[cPos]);
	}
	return score;
}

uint32_t motif::scoreMotif(const std::string::const_iterator & targetBegin,
		const std::string::const_iterator & targetEnd) const {
	if(targetEnd < targetBegin){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": targetBegin must come before targetEnd"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}else if (static_cast<size_t>(targetEnd - targetBegin) != motifUnits_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": motif size doesn't equal size of the check"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	uint32_t score = 0;
	auto strIt = targetBegin;
	auto mapIt = motifUnits_.begin();
	for (; strIt != targetEnd; ++strIt, ++mapIt) {
		score += mapIt->second.scoreChar(*strIt);
	}
	return score;
}

bool motif::passMotifParameter(const std::string & possibleMotif,
		uint32_t scoreCutOff) const {
	return scoreMotif(possibleMotif) >= scoreCutOff;
}

std::vector<size_t> motif::findPositions(const std::string & wholeProtein,
		uint32_t scoreCutOff) const {
	uint32_t motifSize = motifUnits_.size();
	size_t pos = 0;
	std::vector<size_t> positions;
	while (pos + motifSize <= wholeProtein.size()) {
		if (passMotifParameter(wholeProtein.substr(pos, motifSize), scoreCutOff)) {
			positions.emplace_back(pos);
		}
		++pos;
	}
	return positions;
}

std::vector<size_t> motif::findPositionsFull(const std::string & wholeProtein,
		uint32_t allowableErrors) const {
	return findPositionsFull(wholeProtein, allowableErrors, 0,
			wholeProtein.size());
}

std::vector<size_t> motif::findPositionsFull(const std::string & wholeProtein,
		uint32_t allowableErrors, size_t start, size_t stop) const {
	uint32_t sum = 0;
	auto predTest =
			[&sum, &allowableErrors] (const char & a, decltype(*motifUnits_.begin()) mot)
			{
				sum += 1 - mot.second.scoreChar(a);
				return sum <= allowableErrors;
			};
	size_t pos = start;
	std::vector<size_t> positions;
	uint32_t motifSize = motifUnits_.size();
	while (pos + motifUnits_.size() <= stop && pos + motifSize <= wholeProtein.size()) {
		sum = 0;
		if (std::equal(wholeProtein.begin() + pos,
				wholeProtein.begin() + motifSize + pos, motifUnits_.begin(),
				predTest)) {
			positions.emplace_back(pos);
		}
		++pos;
	}
	return positions;
}

size_t motif::size() const {
	return motifUnits_.size();
}

} /* namespace bibseq */
