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
 * alnInfoLocal.cpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */



#include "alnInfoLocal.hpp"
#include "njhseq/utils/vectorUtils.hpp"
namespace njhseq {

void alnInfoLocal::writeInfoSingleLine(std::ostream &indexFile, uint64_t seq1,
                                       uint64_t seq2) const {
  indexFile << seq1 << " " << seq2 << " " << score_ << " " << localAStart_ << " " << localASize_
            << " " << localBStart_ << " " << localBSize_;
  for (const auto &g : gapInfos_) {
    indexFile << " ";
    g.writeInfoNoEnd(indexFile);
  }
  indexFile << std::endl;
}

alnInfoLocal alnInfoLocal::readInfo(std::stringstream& ss,
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
  if (index > 7) {
  	std::vector<uint32_t> positions(index - 7);
  	std::iota(positions.begin(), positions.end(), 7);
    for (auto i : positions) {
      inputGapInfo = tokenizeString(info[i], ",");
      gInfos.emplace_back(gapInfo(std::stoi(inputGapInfo[0]),
                                  std::stoi(inputGapInfo[1]),
                                  std::stoi(inputGapInfo[2])));
    }
  }
  return alnInfoLocal(gInfos, std::stoi(info[3]), std::stoi(info[4]),
                            std::stoi(info[5]), std::stod(info[6]),
                            std::stoi(info[2]), true);
}


Json::Value alnInfoLocal::toJson() const {
	Json::Value ret;
	ret["class"] = "njhseq::alnInfoLocal";
	ret["gapInfos_"] = njh::json::toJson(gapInfos_);
	ret["localAStart_"] = njh::json::toJson(localAStart_);
	ret["localASize_"] = njh::json::toJson(localASize_);
	ret["localBStart_"] = njh::json::toJson(localBStart_);
	ret["localBSize_"] = njh::json::toJson(localBSize_);
	ret["score_"] = njh::json::toJson(score_);
	ret["addFromFile_"] = njh::json::toJson(addFromFile_);
	return ret;
}

alnInfoLocal::alnInfoLocal()
    : localAStart_(0),
      localASize_(0),
      localBStart_(0),
      localBSize_(0),
      score_(0),
      addFromFile_(false)
       {}
// constructor for local alignment
alnInfoLocal::alnInfoLocal(const std::vector<gapInfo>& gInfos, uint32_t localAStart,
             uint32_t localASize, uint32_t localBStart, uint32_t localBSize,
             double score, bool addFromFile)
    : gapInfos_(gInfos),
      localAStart_(localAStart),
      localASize_(localASize),
      localBStart_(localBStart),
      localBSize_(localBSize),
      score_(score),
      addFromFile_(addFromFile) {}

}  // namespace njhseq
