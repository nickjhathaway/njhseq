#pragma once
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

namespace bibseq {

class baseCluster : public readObject {

 public:
  // constructors
  // default constructor, just set everything to zero
  baseCluster() : readObject() {
    seqBase_.cnt_ = 0;
    seqBase_.frac_ = 0;
    cumulativeFraction = 0;
    normalizedFraction = 0;
    needToCalculateConsensus = true;
    updateName();
  }
  // constructor for adding the frist read to the cluster
  baseCluster(const readObject& firstRead) : readObject(firstRead.seqBase_) {
    firstReadName = firstRead.seqBase_.name_;
    firstReadCount = firstRead.seqBase_.cnt_;
    cumulativeFraction = firstRead.cumulativeFraction;
    normalizedFraction = firstRead.normalizedFraction;
    reads_.push_back(firstRead);
    // addRead(firstRead);
    needToCalculateConsensus = true;
    remove = false;
    updateName();
  }
  std::string firstReadName;
  double firstReadCount;
  std::vector<readObject> reads_;
  void addRead(const readObject& newRead);
  std::map<std::string, errorProfile> previousErrorChecks;
  /// consensus
  bool needToCalculateConsensus;
  // vectors to hold the alignemtns to the longest sequence in order to create
  // the consensus
  VecStr longestAlignments;
  VecStr longestAlingmentsRef;
  std::vector<std::vector<uint32_t>> longestAlignmentsQualities;
  // calculated consensus to hold the consensus calculated
  std::string calculatedConsensus;
  std::vector<uint32_t> calculatedConsensusQuality;
  void calculateConsensusNew(aligner& alignerObj, bool setToConsensus);
  void calculateConsensus(aligner& alignerObj, bool setToConsensus = false);
  // consensus comparison
  // get info about the reads in the reads vectors
  VecStr getReadNames() const;
  // writeout
  void writeOutLongestAlignments(const std::string& workingDirectory);
  std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
      calculateAlignmentsToConsensus(aligner& alingerObj);
  void writeOutAlignments(const std::string& directoryName, aligner& alignObj);
  void writeOutClusters(const std::string& directoryName) const;
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
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "baseCluster.cpp"
#endif
