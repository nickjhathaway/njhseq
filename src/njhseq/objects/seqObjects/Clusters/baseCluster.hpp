#pragma once
//
//  baseCluster.hpp
//
//  Created by Nicholas Hathaway on 9/19/13.
//
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

#include "njhseq/objects/seqObjects/readObject.hpp"
#include "njhseq/alignment.h"
#include "njhseq/objects/helperObjects/probabilityProfile.hpp"
#include "njhseq/objects/collapseObjects/opts.h"


namespace njhseq {

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
  baseCluster(const seqInfo& firstRead);

  std::string firstReadName_;
  double firstReadCount_;
  std::vector<std::shared_ptr<readObject>> reads_;
  std::map<std::string, comparison> previousErrorChecks_;
  bool needToCalculateConsensus_;
  seqInfo calcConsensusInfo_;

  void addRead(const baseCluster& newRead);

	struct calculateConsensusPars {
		calculateConsensusPars(bool setToCon) :
				setToConsensus(setToCon) {

		}
		bool setToConsensus { true };
		bool convergeConsensus { false };
		uint32_t convergeAttempts = 5;
	};
  void calculateConsensus(aligner& alignerObj,                            calculateConsensusPars conPars);
  bool calculateConsensusTo(const seqInfo seqBase, aligner& alignerObj, calculateConsensusPars conPars);
  void calculateConsensusToCurrent(aligner& alignerObj,                   calculateConsensusPars conPars);
  // consensus comparison
  // get info about the reads in the reads vectors
  VecStr getReadNames() const;
  std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
      calculateAlignmentsToConsensus(aligner& alingerObj);
  void writeOutAlignments(const std::string& directoryName, aligner& alignObj);

  void writeClusters(const SeqIOOptions & ioOptions) const;
  void writeClustersInDir(const std::string &workingDir, const SeqIOOptions & ioOptions) const;

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
      alignObj.alignCache(*this, read, local);
      prof.increaseCountAmount(alignObj.alignObjectA_.seqBase_.seq_,
                               alignObj.alignObjectB_.seqBase_.seq_,
                               read->seqBase_.cnt_);
    }
  }
  bool compare(baseCluster & read, aligner & alignerObj,
  		const IterPar & runParams,
  		const CollapserOpts & collapserOptsObj);

  comparison getComparison(baseCluster & read, aligner & alignerObj, bool checkKmers) const;


  readObject createRead() const;

	bool isClusterCompletelyChimeric();
	bool isClusterAtLeastHalfChimeric();
	bool isClusterAtLeastChimericCutOff(double cutOff);


  void removeRead(const std::string &stubName);
  void removeRead(uint32_t readPos);
  void removeReads(std::vector<uint32_t> readPositions);

  // converter for all clusters
  template <class CLUSTER, class READ>
  static std::vector<CLUSTER> convertVectorToClusterVector(
      const std::vector<READ>& reads) {
    std::vector<CLUSTER> ans;
    for (const auto& read : reads) {
      ans.emplace_back(CLUSTER(getSeqBase(read)));
    }
    return ans;
  }
};
}  // namespace njhseq


