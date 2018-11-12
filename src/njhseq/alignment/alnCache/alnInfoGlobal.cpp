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
/*
 * alnInfoGlobal.cpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */




#include "alnInfoGlobal.hpp"
#include "njhseq/utils/vectorUtils.hpp"

namespace njhseq {


// empty construcor
alnInfoGlobal::alnInfoGlobal() :score_(0), addFromFile_(false) {}
// contructor for global alingment
alnInfoGlobal::alnInfoGlobal(const std::vector<gapInfo>& gInfos, double score, bool addFromFile)
    : gapInfos_(gInfos),score_(score), addFromFile_(addFromFile) {}

void alnInfoGlobal::writeInfoSingleLine(std::ostream &indexFile, uint64_t seq1,
                                        uint64_t seq2) const {
  indexFile << seq1 << " " << seq2 << " " << score_;
  for (const auto &g : gapInfos_) {
    indexFile << " ";
    g.writeInfoNoEnd(indexFile);
  }
  indexFile << std::endl;
}

alnInfoGlobal alnInfoGlobal::readInfo(std::stringstream& ss,
   		std::vector<std::string>& info, std::vector<std::string>& inputGapInfo,
   		std::vector<gapInfo> & gInfos, std::string & out){

	uint32_t index = 2;
  while (!ss.eof()) {
    ss >> out;
    if(index >= info.size()){
    	addOtherVec(info, VecStr(index));
    }
    info[index] = out;
    ++index;
  }

  if (index > 3) {
  	std::vector<uint32_t> positions(index - 3);
  	std::iota(positions.begin(), positions.end(), 3);
    for (auto i : positions) {
      inputGapInfo = tokenizeString(info[i], ",");
      gInfos.emplace_back(gapInfo(std::stoi(inputGapInfo[0]),
                                  std::stoi(inputGapInfo[1]),
                                  std::stoi(inputGapInfo[2])));
    }
  }
  return alnInfoGlobal(gInfos, std::stod(info[2]), true);
}

Json::Value alnInfoGlobal::toJson() const {
	Json::Value ret;
	ret["class"] = "njhseq::alnInfoGlobal";
	ret["gapInfos_"] = njh::json::toJson(gapInfos_);
	ret["score_"] = njh::json::toJson(score_);
	ret["addFromFile_"] = njh::json::toJson(addFromFile_);
	return ret;
}


}  // namespace njhseq
