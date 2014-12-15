//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#pragma once
/*

 * strCounterMap.hpp
 *
 *  Created on: Jun 1, 2014
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"

namespace bibseq {

class strCounterMap {

public:
	//constructor
	//members
  std::unordered_map<std::string, uint32_t> counts_;
  std::unordered_map<std::string, double> fractions_;
  //std::unordered_map<std::string, std::vector<uint32_t>> qualities_;
  //std::unordered_map<std::string, std::vector<std::vector<uint32_t>>> allQualities_;



  //functions
  void increaseCountByString(const std::string &seq);
  virtual void increaseCountByString(const std::string &seq, double cnt);
  void increaseCountByVecStr(const VecStr &seqs);
  void increaseCountByVecStr(const VecStr &seqs, const std::vector<double> & counts);
  std::multimap<double, std::string, std::less<double>> createLikelihoodMaps(
        bool setFractionFirst);
  void setFractions();
  void printAllInfo(std::ostream & out);
  void printCountInfo(std::ostream & out);
  void printFractionInfo(std::ostream & out);
};

} /* namespace bib */

#ifndef NOT_HEADER_ONLY
#include "strCounterMap.cpp"
#endif
