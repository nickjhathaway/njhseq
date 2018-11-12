#pragma once
/*
 * PairedCluster.hpp
 *
 *  Created on: Jul 27, 2015
 *      Author: nick
 */
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

#include "njhseq/objects/seqObjects/Paired/PairedRead.hpp"
#include "njhseq/alignment/aligner.h"
#include "njhseq/objects/collapseObjects/opts.h"

namespace njhseq {

class PairedCluster : public PairedRead {
public:

	PairedCluster(const PairedRead& firstRead) :
			PairedRead(firstRead.seqBase_, firstRead.mateSeqBase_, false) {
		firstReadName = firstRead.seqBase_.name_;
		firstReadCount = firstRead.seqBase_.cnt_;
		reads_.push_back(std::make_shared<PairedRead>(firstRead));
		setLetterCount();
		needToCalculateConsensus = false;
		remove = false;
		updateName();
/*
		if(seqBase_.name_ == "lib_00.278-2/1_t1"){
			seqBase_.outPutFastq(std::cout);
			mateSeqBase_.outPutFastq(std::cout);
			for(const auto & read : reads_){
				read.seqBase_.outPutFastq(std::cout);
				read.mateSeqBase_.outPutFastq(std::cout);
			}
			exit(1);
		}*/

	}


	//seqInfo mateSeqBase_;

	std::vector<std::shared_ptr<PairedRead>> reads_;

	//uint32_t firstAlnFlag_ = 0;
	//uint32_t secondAlnFlag_ = 0;


	void addRead(const PairedCluster & otherRead);
  std::string firstReadName;
  double firstReadCount;
  std::map<std::string, std::pair<comparison,comparison>> previousErrorChecks;
  bool needToCalculateConsensus;

	void calculateConsensus(aligner & alignerObj, bool setToConsensus);

	bool compare(PairedCluster & otherRead,
			aligner & alignerObj,
			const IterPar & errorThreshold,
			const CollapserOpts & collapserOptsObj);

	bool compareRead(PairedCluster & otherRead,
			aligner & alignerObj,
			const IterPar & errorThreshold,
			const CollapserOpts & collapserOptsObj);

	bool compareMate(PairedCluster & otherRead,
			aligner & alignerObj,
			const IterPar & errorThreshold,
			const CollapserOpts & collapserOptsObj);



	void writeOutClusters(const std::string& directoryName,const SeqIOOptions & ioOptions) const;
	void writeOutClustersWithConsensus(std::ostream & firstOut,std::ostream & secondOut )const;

	virtual ~PairedCluster(){}
	using size_type = std::string::size_type;
};

template<>
inline seqInfo::size_type len(const PairedCluster & read){
	return len(read.seqBase_) + len(read.mateSeqBase_);
}


void collapseIdenticalPairedClusters(std::vector<PairedRead> & clusters, std::string qualRep);

void clusterPairedClusters(std::vector<PairedCluster> & clusters,
		aligner & alignerObj,
		aligner & alignerObjMate,
		comparison & errorThreshold,
		const CollapserOpts & collapserOptsObj);




std::vector<PairedCluster> runClustering(std::vector<PairedCluster> currentClusters,
		const CollapserOpts & collapserOptsObj, aligner & alignerObj, aligner & alignerObjMate,
		comparison & errorThreshold);


} /* namespace njhseq */


