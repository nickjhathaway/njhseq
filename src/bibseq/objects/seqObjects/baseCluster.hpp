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
//  baseCluster.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/19/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/alignment.h"
#include "bibseq/IO/readObjectIO.hpp"
#include "bibseq/objects/helperObjects/probabilityProfile.hpp"
#include "bibseq/simulation.h"

namespace bibseq {

class baseCluster : public readObject {

 public:
  // constructors
  // default constructor, just set everything to zero
  baseCluster() : readObject() {
    seqBase_.cnt_ = 0;
    seqBase_.frac_ = 0;
    firstReadCount = 0;
    needToCalculateConsensus = true;
    updateName();
  }
  // constructor for adding the frist read to the cluster
  baseCluster(const readObject& firstRead);
  std::string firstReadName;
  double firstReadCount;
  std::vector<readObject> reads_;
  void addRead(const readObject& newRead);
  std::map<std::string, comparison> previousErrorChecks;
  /// consensus
  bool needToCalculateConsensus;
  // vectors to hold the alignments to the longest sequence in order to create
  // the consensus
  VecStr longestAlignments;
  VecStr longestAlingmentsRef;
  std::vector<std::vector<uint32_t>> longestAlignmentsQualities;
  // calculated consensus to hold the consensus calculated
  std::string calculatedConsensus;
  std::vector<uint32_t> calculatedConsensusQuality;
  seqInfo calcConsensusInfo_;
  void calculateConsensusOldLet(aligner& alignerObj, bool setToConsensus);
  void calculateConsensus(aligner& alignerObj, bool setToConsensus);
  void calculateConsensusTo(const seqInfo & seqBase, aligner& alignerObj, bool setToConsensus);
  void calculateConsensusOld(aligner& alignerObj, bool setToConsensus);
  void calculateConsensusToOld(const seqInfo & seqBase, aligner& alignerObj, bool setToConsensus);
  void calculateConsensusToCurrent(aligner& alignerObj, bool setToConsensus);
  // consensus comparison
  // get info about the reads in the reads vectors
  VecStr getReadNames() const;
  // writeout
  void writeOutLongestAlignments(const std::string& workingDirectory);
  std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
      calculateAlignmentsToConsensus(aligner& alingerObj);
  void writeOutAlignments(const std::string& directoryName, aligner& alignObj);
  void writeOutClusters(const std::string& directoryName,const readObjectIOOptions & ioOptions) const;
  // get infos
  double getAverageReadLength() const;

  template <typename P>
  void updateErrorProfile(P& prof, aligner& alignObj, bool local) const {
    for (const auto& read : reads_) {
      alignObj.alignVec(*this, read, local);
      prof.increaseCountAmount(alignObj.alignObjectA_.seqBase_.seq_,
                               alignObj.alignObjectB_.seqBase_.seq_,
                               read.seqBase_.cnt_);
    }
  }

  readObject createRead() const;
  virtual void printDescription(std::ostream& out, bool deep = false) const;
  // converter for all clusters
  template <class CLUSTER, class READ>
  static std::vector<CLUSTER> convertVectorToClusterVector(
      const std::vector<READ>& reads) {
    std::vector<CLUSTER> ans;
    for (const auto& read : reads) {
      ans.emplace_back(CLUSTER(read));
    }
    return ans;
  }
};
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "baseCluster.cpp"
#endif
