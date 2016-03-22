#pragma once
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
//
//  clusterCollapser.h
//  sequenceTools
//
//  Created by Nick Hathaway on 11/16/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils/utils.hpp"
#include "bibseq/objects/seqObjects/Clusters/identicalCluster.hpp"
#include "bibseq/objects/seqObjects/Clusters/cluster.hpp"
#include "bibseq/alignment.h"




namespace bibseq {

class clusterCollapser {

 public:
  clusterCollapser() {}
  template<typename T>
  static std::vector<identicalCluster> collapseToUniqueReads(
        const std::vector<T> &input, const std::string &repQual);

  template<typename T>
  static std::vector<identicalCluster> collapseIdenticalReads(
      const std::vector<T> &reads, const std::string &repQual){
    std::vector<identicalCluster> ret;
    uint32_t count = 0;
    for (const auto &read : reads) {
    	if(!getSeqBase(read).on_){
    		continue;
    	}
      ++count;
      if (count == 1) {
      	ret.emplace_back(getSeqBase(read));
        continue;
      }
      bool foundMatch = false;
      for (auto &clusterIter : ret) {
        if (getSeqBase(read).seq_ == clusterIter.seqBase_.seq_) {
          clusterIter.addRead(getSeqBase(read));
          foundMatch = true;
          break;
        }
      }
      if (!foundMatch) {
      	ret.emplace_back(getSeqBase(read));
      }
    }
    identicalCluster::setIdneticalClusterQual(ret, repQual);
    readVec::allSetLetterCount(ret);
    readVec::allUpdateName(ret);
    return ret;
  }

  // cluster down on tandems
  static void collapseTandems(std::vector<cluster> &processedReads,
                              aligner &alignerObj, int runCutOff, int kLength,
                              bool kMersByPosition, double freqCutoff,
                              bool local, bool weighHomopolyer);
};
template<typename T>
std::vector<identicalCluster> clusterCollapser::collapseToUniqueReads(
        const std::vector<T> &reads, const std::string &repQual){
	std::vector<identicalCluster> output;
	std::unordered_map<std::string, size_t> clusters;
	for(const auto & read : reads){
		auto search = clusters.find(getSeqBase(read).seq_);
		if(search == clusters.end()){
			clusters[getSeqBase(read).seq_] = output.size();
			output.emplace_back(getSeqBase(read));
		}else{
			output[search->second].addRead(getSeqBase(read));
		}
	}
	identicalCluster::setIdneticalClusterQual(output, repQual);
	return output;
}



}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "clusterCollapser.cpp"
#endif
