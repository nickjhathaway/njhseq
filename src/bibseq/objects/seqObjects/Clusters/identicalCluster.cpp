#include "identicalCluster.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers.h"
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
namespace njhseq {

identicalCluster::identicalCluster(const seqInfo& firstRead) : baseCluster(firstRead) {

}


void identicalCluster::addRead(const readObject& identicalRead) {
  reads_.emplace_back(std::make_shared<readObject>(identicalRead));
  seqBase_.cnt_ += identicalRead.seqBase_.cnt_;
}
////////setting of the representive quality and seq
void identicalCluster::setSeq() {
  readVecSorter::sort(reads_);
  seqBase_.seq_ = reads_.front()->seqBase_.seq_;
  seqBase_.name_ = reads_.front()->seqBase_.name_;
  updateName();
}
void identicalCluster::setRep(const std::string& repQual){
	if (repQual == "worst") {
		setWorstQualRep();
	} else if (repQual == "median") {
		setMedianQualRep();
	} else if (repQual == "average") {
		setAverageQualRep();
	} else if (repQual == "bestSeq") {
		setBestSeqRep();
	} else if (repQual == "bestQual") {
		setBestQualRep();
	} else {
		std::stringstream ss;
		ss << "Unrecognized qualRep: " << repQual << std::endl;
		ss << "Needs to be median, average, bestSeq, bestQual, or worst"
							<< std::endl;
		throw std::runtime_error{ss.str()};
	}
}
void identicalCluster::setBestQualRep() {
  setSeq();
  seqBase_.qual_.clear();
  for (auto i : iter::range(seqBase_.seq_.length())) {
    uint32_t currentQual = 0;
    for (const auto& read : reads_) {
      if (read->seqBase_.qual_[i] > currentQual) {
        currentQual = read->seqBase_.qual_[i];
      }
    }
    seqBase_.qual_.push_back(currentQual);
  }
}
void identicalCluster::setWorstQualRep() {
  setSeq();
  seqBase_.qual_.clear();
  for (auto i : iter::range(seqBase_.seq_.length())) {
    uint32_t currentQual = UINT32_MAX;
    for (const auto& read : reads_) {
      if (read->seqBase_.qual_[i] < currentQual) {
        currentQual = read->seqBase_.qual_[i];
      }
    }
    seqBase_.qual_.push_back(currentQual);
  }
}

void identicalCluster::setBestSeqRep() {
  setSeq();
  seqBase_.qual_.clear();
  seqBase_.qual_ = reads_.front()->seqBase_.qual_;
}

void identicalCluster::setAverageQualRep() {
  setSeq();
  seqBase_.qual_.clear();
  for (auto i : iter::range(seqBase_.seq_.length())) {
    int qualSum = 0;
    for (const auto& read : reads_) {
      qualSum += read->seqBase_.qual_[i];
    }
    seqBase_.qual_.push_back((int)(qualSum / seqBase_.cnt_));
  }
}
void identicalCluster::setMedianQualRep() {
  setSeq();
  seqBase_.qual_.clear();
  for (auto i : iter::range(seqBase_.seq_.length())) {
    std::vector<uint32_t> qualities;
    for (auto read : reads_) {
      qualities.push_back(read->seqBase_.qual_[i]);
    }
    uint32_t qualToInsert = vectorMedianRef(qualities);
    seqBase_.qual_.push_back(qualToInsert);
  }
}
}  // namespace njh
