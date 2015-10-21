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
#include "bibseq/helpers/alignmentProfiler.hpp"
#include "bibseq/objects/helperObjects/probabilityProfile.hpp"
#include "bibseq/simulation.h"
#include "bibseq/objects/collapseObjects/collapserOpts.hpp"

namespace bibseq {

class baseCluster : public readObject {

 public:
  // constructors
	/**@brief Empty constructor
	 *
	 */
  baseCluster();
  /**@brief Constructor for adding the first read to the cluster
   *
   * @param firstRead The seed of this cluster
   */
  baseCluster(const readObject& firstRead);

  std::string firstReadName_;
  double firstReadCount_;
  std::vector<readObject> reads_;
  std::map<std::string, comparison> previousErrorChecks_;
  bool needToCalculateConsensus_;
  seqInfo calcConsensusInfo_;

  void addRead(const readObject& newRead);

  void calculateConsensus(aligner& alignerObj, bool setToConsensus);
  void calculateConsensusTo(const seqInfo & seqBase, aligner& alignerObj, bool setToConsensus);
  void calculateConsensusToCurrent(aligner& alignerObj, bool setToConsensus);
  // consensus comparison
  // get info about the reads in the reads vectors
  VecStr getReadNames() const;
  std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
      calculateAlignmentsToConsensus(aligner& alingerObj);
  void writeOutAlignments(const std::string& directoryName, aligner& alignObj);
  void writeOutClusters(const std::string& directoryName,const readObjectIOOptions & ioOptions) const;
  void alignmentProfile(const std::string& workingDirectory,
                        aligner& alignerObj, bool local, bool kmerChecking,
                        int kLength, bool kmersByPosition,
                        bool weighHomopolyers) const;
  // get infos
  double getAverageReadLength() const;

  template <typename P>
  void updateErrorProfile(P& prof, aligner& alignObj, bool local) const {
    /*if (reads.size() == 1) {
     // if reads size is 1, consensus should be that read and shouldn't be any
     // different
     return;
     }*/
    for (const auto& read : reads_) {
      /*if (read.seqBase_.seq_ == seqBase_.seq_) {
       // if the seqs are equal no point in looking for mismatches
       continue;
       }*/
      alignObj.alignVec(*this, read, local);
      prof.increaseCountAmount(alignObj.alignObjectA_.seqBase_.seq_,
                               alignObj.alignObjectB_.seqBase_.seq_,
                               read.seqBase_.cnt_);
    }
  }
  bool compare(baseCluster & read, aligner & alignerObj,
  		const comparison & errorThreshold,
  		const collapserOpts & collapserOptsObj);

  bool compareId(baseCluster & read, aligner & alignerObj,
  		const comparison & errorThreshold,
  		const collapserOpts & collapserOptsObj);

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
