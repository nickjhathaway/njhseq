#include "cluster.hpp"
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
#include "bibseq/IO/SeqIO/SeqOutput.hpp"

namespace bibseq {

cluster::cluster() : baseCluster() { rejected_ = false; }

cluster::cluster(const seqInfo &firstRead) : baseCluster(firstRead) {
  rejected_ = false;
  setLetterCount();
}

void cluster::removeRead(const std::string & stubName){
	uint32_t readPos = std::numeric_limits<uint32_t>::max();
	for(const auto & pos : iter::range(reads_.size())){
		if(reads_[pos]->getStubName(true) == stubName){
			readPos = pos;
			break;
		}
	}
	if(readPos!=std::numeric_limits<uint32_t>::max()){
		seqBase_.cnt_ -= reads_[readPos]->seqBase_.cnt_;
		reads_.erase(reads_.begin() + readPos);
	  if (seqBase_.cnt_ / 2 > firstReadCount_) {
	    needToCalculateConsensus_ = true;
	  }
	}
}
void cluster::removeReads(const std::vector<readObject> & vec){
	for(const auto & read : vec){
		removeRead(read.getStubName(true));
	}
	updateName();
}


void cluster::addRead(const cluster& cr) {
  reads_.insert(reads_.end(), cr.reads_.begin(), cr.reads_.end());
  seqBase_.cnt_ += cr.seqBase_.cnt_;
  seqBase_.frac_ = seqBase_.frac_ * reads_.size();
  seqBase_.frac_ = (seqBase_.frac_ + cr.seqBase_.frac_) / reads_.size();
  // needToCalculateConsensus = true;
  if (seqBase_.cnt_ / 2.0 >= firstReadCount_) {
    needToCalculateConsensus_ = true;
  }
}

int cluster::getBiggestReadSize() {
  int biggestSize = 0;
  for (const auto& read : reads_) {
    if (read->seqBase_.cnt_ > biggestSize) {
      biggestSize = read->seqBase_.cnt_;
    }
  }
  return biggestSize;
}
double cluster::getAverageSizeDifference() {
  double averageSizeDifference = 0;
  for (const auto& read : reads_) {
    averageSizeDifference +=
        uAbsdiff(read->seqBase_.seq_.length(), seqBase_.seq_.length()) *
        read->seqBase_.cnt_;
  }
  averageSizeDifference = averageSizeDifference / seqBase_.cnt_;
  return averageSizeDifference;
}

size_t cluster::getLargestSizeDifference() {
  size_t largestSizeDifference = 0;
  for (const auto& read : reads_) {
    if (uAbsdiff(read->seqBase_.seq_.length(),seqBase_.seq_.length()) >
        largestSizeDifference) {
      largestSizeDifference =
      		uAbsdiff(read->seqBase_.seq_.length(), seqBase_.seq_.length());
    }
  }
  return largestSizeDifference;
}


void cluster::outputInfoComp(const std::string& workingDir) const {
  std::ofstream info(bib::files::make_path(workingDir, seqBase_.name_).string());
  if (!info) {
    std::cout << "Error in opening" << workingDir << seqBase_.name_
              << std::endl;
  }
  std::map<int, int> clusterCounter;
  for (const auto& rIter : reads_) {
    ++clusterCounter[rIter->seqBase_.cnt_];
  }
  info << "clusterSize\tFrequency" << std::endl;
  for (const auto& kv : clusterCounter) {
    info << kv.first << "\t" << kv.second << std::endl;
  }
  return;
}




}  // namespace bib
