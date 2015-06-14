#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//  identicalCluster.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/1/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/seqObjects/baseCluster.hpp"

namespace bibseq {

class identicalCluster : public baseCluster {
 public:
  // constructors

  identicalCluster() : baseCluster() {}

  identicalCluster(const readObject& firstRead) ;
  template <typename T>
  identicalCluster(const std::vector<T>& reads, const std::string & qualRep) : baseCluster(readObject(reads.front().seqBase_)) {
  	for(const auto & readPos : iter::range<uint32_t>(1, len(reads))){
  		addRead(readObject(reads[readPos].seqBase_));
  	}
  	setRep(qualRep);
  	updateName();
  }

  virtual void addRead(const readObject& identicalRead);
  void setRep(const std::string& repQual);
  // set the quality and seq to represent the cluster
  void setSeq();
  void setBestQualRep();
  void setWorstQualRep();
  void setBestSeqRep();
  void setAverageQualRep();
  void setMedianQualRep();
  //
  template <typename T>
  static void setIdneticalClusterQual(std::vector<T>& reads,
                                      const std::string& repQual) {
    if (repQual == "worst") {
      std::for_each(reads.begin(), reads.end(),
                    [](T& read) { read.setWorstQualRep(); });
    } else if (repQual == "median") {
      std::for_each(reads.begin(), reads.end(),
                    [](T& read) { read.setMedianQualRep(); });
    } else if (repQual == "average") {
      std::for_each(reads.begin(), reads.end(),
                    [](T& read) { read.setAverageQualRep(); });
    } else if (repQual == "bestSeq") {
      std::for_each(reads.begin(), reads.end(),
                    [](T& read) { read.setBestSeqRep(); });
    } else if (repQual == "bestQual") {
      std::for_each(reads.begin(), reads.end(),
                    [](T& read) { read.setBestQualRep(); });
    } else {
    	std::stringstream ss;
      ss << "Unrecognized qualRep: " << repQual << std::endl;
      ss << "Needs to be median, average, bestSeq, bestQual, or worst"
                << std::endl;
      throw std::runtime_error{ss.str()};
    }
  }
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "identicalCluster.cpp"
#endif
