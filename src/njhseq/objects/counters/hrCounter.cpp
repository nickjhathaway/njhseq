/*
 * hrCounter.cpp
 *
 *  Created on: May 23, 2014
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
#include "hrCounter.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"

namespace njhseq {

///str map
std::multimap<double, std::string, std::less<double>> hrCounter::createLikelihoodMaps(
      bool setFractionFirst){
  if (setFractionFirst) {
    setFractions();
  }
  std::multimap<double, std::string, std::less<double>> likelihoods;
  for (const auto &str : fractions_) {
  	for (const auto & run : str.second){
  		likelihoods.emplace(run.second, std::string(run.first, str.first));
  	}
  }
  return likelihoods;
}
void hrCounter::increaseCount(const std::string &seq){
	condensedSeq cSeq(seq);
	increaseCount(cSeq);
}

void hrCounter::increaseCount(const std::string &seq, double cnt){
	condensedSeq cSeq(seq);
	//cSeq.cnt_ = std::round(cnt);
	cSeq.cnt_ = cnt;
	increaseCount(cSeq);
}
void hrCounter::increaseCount(const condensedSeq & cSeq){
	for(auto i : iter::range(cSeq.seq_.size())){
		hCounts_[cSeq.seq_[i]][cSeq.counts_[i]] += cSeq.cnt_;
	}
}
void hrCounter::reset(){
	hCounts_.clear();
	fractions_.clear();
	hQuals_.clear();
}

void hrCounter::setFractions(){
	uint32_t total = 0;
	fractions_.clear();
	for(const auto & let : hCounts_){
		for(const auto & run : let.second){
			total += run.second;
		}
	}
	if(total != 0){
		for(const auto & let : hCounts_){
			for(const auto & run : let.second){
				fractions_[let.first][run.first] = run.second /static_cast<double> (total);
			}
		}
	}
}

void hrCounter::printAllInfo(std::ostream & out){
	printCountInfo(std::cout);
	printFractionInfo(std::cout);
}

void hrCounter::printCountInfo(std::ostream & out){
	table outTab(hCounts_, VecStr{"char", "hpRunSize", "count"});
	outTab.sortTable("hpRunSize", true);
	outTab.sortTable("char", true);
	outTab.outPutContentOrganized(out);
}

void hrCounter::printFractionInfo(std::ostream & out){
	table outTab(fractions_, VecStr{"char", "hpRunSize", "fraction"});
	outTab.sortTable("hpRunSize", true);
	outTab.sortTable("char", true);
	outTab.outPutContentOrganized(out);
}

condensedSeq::condensedSeq(const std::string & seq){
  uint32_t currentCount = 1;
  uint32_t i = 1;

  for (; i < seq.length(); i++) {
    if (seq[i] == seq[i - 1]) {
      ++currentCount;
    } else {
      seq_.push_back(seq[i - 1]);
      counts_.push_back(currentCount);
      currentCount = 1;
    }
  }

  seq_.push_back(seq[i - 1]);
  counts_.push_back(currentCount);
}

} /* namespace njh */
