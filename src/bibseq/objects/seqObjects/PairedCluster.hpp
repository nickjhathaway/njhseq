#pragma once
/*
 * PairedCluster.hpp
 *
 *  Created on: Jul 27, 2015
 *      Author: nick
 */

#include "bibseq/objects/seqObjects/PairedRead.hpp"
#include "bibseq/alignment/aligner.hpp"
#include "bibseq/objects/collapseObjects/collapserOpts.hpp"

namespace bibseq {

class PairedCluster : public PairedRead {
public:

	PairedCluster(const PairedRead& firstRead) :
			PairedRead(firstRead.seqBase_, firstRead.mateSeqBase_, false) {
		firstReadName = firstRead.seqBase_.name_;
		firstReadCount = firstRead.seqBase_.cnt_;
		reads_.push_back(firstRead);
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

	std::vector<PairedRead> reads_;

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
			const comparison & errorThreshold,
			const collapserOpts & collapserOptsObj);

	bool compareRead(PairedCluster & otherRead,
			aligner & alignerObj,
			const comparison & errorThreshold,
			const collapserOpts & collapserOptsObj);

	bool compareMate(PairedCluster & otherRead,
			aligner & alignerObj,
			const comparison & errorThreshold,
			const collapserOpts & collapserOptsObj);

	bool compareId(PairedCluster & otherRead,
			aligner & alignerObj,
			const comparison & errorThreshold,
			const collapserOpts & collapserOptsObj);

	bool compareReadId(PairedCluster & otherRead,
			aligner & alignerObj,
			const comparison & errorThreshold,
			const collapserOpts & collapserOptsObj);

	bool compareMateId(PairedCluster & otherRead,
			aligner & alignerObj,
			const comparison & errorThreshold,
			const collapserOpts & collapserOptsObj);

	void writeOutClusters(const std::string& directoryName,const readObjectIOOptions & ioOptions) const;
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
		const collapserOpts & collapserOptsObj);




std::vector<PairedCluster> runClustering(std::vector<PairedCluster> currentClusters,
		const collapserOpts & collapserOptsObj, aligner & alignerObj, aligner & alignerObjMate,
		comparison & errorThreshold);


} /* namespace bibseq */


