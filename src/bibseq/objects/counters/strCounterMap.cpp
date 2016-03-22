/*
 * strCounterMap.cpp
 *
 *  Created on: Jun 1, 2014
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
#include "strCounterMap.hpp"

namespace bibseq {

void strCounter::addOtherCounts(const strCounter & otherCounter) {
	for (const auto & strCount : otherCounter.counts_) {
		counts_[strCount.first] += strCount.second;
	}
}

///str map
std::multimap<double, std::string, std::less<double>> strCounter::createLikelihoodMaps(
		bool setFractionFirst) {
	if (setFractionFirst) {
		setFractions();
	}
	std::multimap<double, std::string, std::less<double>> likelihoods;
	for (const auto &str : fractions_) {
		likelihoods.emplace(str.second, str.first);
	}
	return likelihoods;
}
void strCounter::increaseCountByString(const std::string &seq) {
	increaseCountByString(seq, 1);
}
void strCounter::increaseCountByVecStr(const VecStr &seqs) {
	for (const auto & str : seqs) {
		increaseCountByString(str, 1);
	}
}
void strCounter::increaseCountByVecStr(const VecStr &seqs,
		const std::vector<double> & counts) {
	if (counts.size() == 1) {
		for (const auto & str : seqs) {
			increaseCountByString(str, counts.front());
		}
	} else if (counts.size() == seqs.size()) {
		for (const auto & strPos : iter::range(seqs.size())) {
			increaseCountByString(seqs[strPos], counts[strPos]);
		}
	} else {
		throw std::runtime_error {
				"VecStr and counts should be the same length or counts should be only one number" };
	}
}

void strCounter::increaseCountByString(const std::string &seq, double cnt) {
	counts_[seq] += cnt;
}

uint32_t strCounter::getTotalCount() const {
	uint32_t ret = 0;
	for (const auto & let : counts_) {
		ret += let.second;
	}
	return ret;
}

void strCounter::setFractions() {
	auto total = getTotalCount();
	if (total != 0) {
		for (const auto & let : counts_) {
			fractions_[let.first] = let.second / static_cast<double>(total);
		}
	}
}
void strCounter::printAllInfo(std::ostream & out) const {
	printCountInfo(std::cout);
	printFractionInfo(std::cout);
}
void strCounter::printCountInfo(std::ostream & out) const {
	std::cout << "strs" << std::endl;
	for (const auto & s : counts_) {
		std::cout << s.first << " : " << s.second << std::endl;
	}
}
void strCounter::printFractionInfo(std::ostream & out) const {
	std::cout << "fractions" << std::endl;
	for (const auto & s : fractions_) {
		std::cout << s.first << " : " << s.second << std::endl;
	}
}

} /* namespace bib */
