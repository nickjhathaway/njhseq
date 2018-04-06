#pragma once
/*
 * clusterUtils.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: nick
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bibseq/objects/seqObjects/Clusters/baseCluster.hpp"
#include "bibseq/readVectorManipulation/readVectorOperations/massSetters.hpp"
#include "bibseq/readVectorManipulation/readVectorOperations/massGetters.hpp"

namespace bibseq {
namespace clusterVec {

template<typename T>
void allWriteClustersInDir(const std::vector<T> &reads,
                                     const std::string &workingDir,
																		 const SeqIOOptions & opts) {
  bib::for_each(reads,
           [&opts,&workingDir](const T &read) { read.writeClustersInDir(workingDir, opts); });
}



template<typename T>
void allCalculateConsensus(std::vector<T> &reads, aligner &alignerObj,
                                  bool setToConsensus) {
	bib::for_each(reads, [&alignerObj,&setToConsensus](T &read) {
    read.calculateConsensus(alignerObj, setToConsensus);
  });
}

template<typename T>
void allSetFractionClusters(std::vector<T> &reads) {
  size_t sizeOfRead = readVec::getTotalReadCount(reads);
  int count = 0;
  for (auto &read : reads) {
    read.setFractionByCount(sizeOfRead);
    ++count;
    for (auto &subSub : read.reads_) {
      subSub->setFractionByCount(sizeOfRead);
    }
  }
}

template<typename T>
std::string returnSuperFromSub(const std::vector<T> &reads,
                                      const std::string &subName) {
  std::string notfound = "";
  for (const auto &read : reads) {
    for (const auto &sread : read.reads_) {
      if (subName == sread->seqBase_.name_) {
        return read.seqBase_.name_;
      }
    }
  }
  return notfound;
}


}  // namespace clusterVec
}  // namespace bibseq



