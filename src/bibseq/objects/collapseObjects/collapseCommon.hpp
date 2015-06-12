#pragma once
//
//  collapseCommon.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/sampleCluster.hpp"
#include "bibseq/helpers/clusterCollapser.hpp"
#include "bibseq/objects/collapseObjects/collapser.hpp"
#include "bibseq/readVectorManipulation.h"
#include "bibseq/seqToolsUtils.h"
#include "bibseq/helpers/profiler.hpp"

namespace bibseq {
namespace collapse {

class clusterSetInfo {
public:
	clusterSetInfo(){}
	template<typename T>
	clusterSetInfo(const std::vector<T> & reads){
		setInfo(reads);
	}
	std::map<std::string, sampInfo> infos_;
	double totalReadCount_ = 0;
	uint32_t numberOfClusters_ = 0;
	std::vector<uint32_t> mois_;

	template<typename T>
	void setInfo(const std::vector<T> & reads){
		clear();
		for (const auto & read : reads){
			updateInfo(read);
		}
		//for cluster set info set the runReadCnt equal to the readCnt count
		resetRunReadCnt();
	}

	template<typename T>
	void updateInfo(const std::vector<T> & reads){
		for (const auto & read : reads){
			updateInfo(read);
		}
		//for cluster set info set the runReadCnt equal to the readCnt count
		resetRunReadCnt();
	}
	template<typename T>
	void updateInfo(const T & read){
		auto search = infos_.find(read.getOwnSampName());
		if(search == infos_.end()){
			infos_.emplace(read.getOwnSampName(), sampInfo(read));
		}else{
			search->second.update(read);
		}
		totalReadCount_ += read.seqBase_.cnt_;
		++numberOfClusters_;
	}


	void clear(){
		std::map<std::string, sampInfo>().swap(infos_);
		totalReadCount_ = 0;
		numberOfClusters_ = 0;
	}

	void resetRunReadCnt(){
		for(auto & i : infos_){
			i.second.runReadCnt_ = i.second.readCnt_;
		}
	}

	void updateMoi(uint32_t moi){
		mois_.emplace_back(moi);
	}

	double meanMoi()const{
		return vectorMean(mois_);
	}
	double medianMoi()const{
		return vectorMedian(mois_);
	}
	uint32_t maxMoi()const {
		return vectorMaximum(mois_);
	}

	uint32_t minMoi()const {
		return vectorMinimum(mois_);
	}

	static std::string getMoiInfoHeader(const std::string & delim){
		return vectorToString(toVecStr("p_meanMoi", "p_medianMoi", "p_minMoi", "p_maxMoi"), delim);
	}

	std::string getMoiInfo(const std::string & delim)const {
		return vectorToString(toVecStr(meanMoi(), medianMoi(), minMoi(), maxMoi()), delim);
	}
};




class clusterSet {
public:
	// constructors
	clusterSet() {
	}
	clusterSet(const std::vector<sampleCluster>& clusters) :
			clusters_(clusters), info_(clusters) {
	}
	// members
	std::vector<sampleCluster> clusters_;
	std::unordered_map<std::string, uint32_t> subClustersPositions_;

	clusterSetInfo info_;

	// functions
	void setSubClustersPositions();
	void setSetInfo();

	template<typename REF>
	void checkAgainstExpected(const std::vector<REF>& refSeqs,
			aligner& alignerObj, bool local, bool weighHomopolyers) {
		bool eventBased = true;
		for (auto& clus : clusters_) {

			std::string bestRef = profiler::getBestRef(refSeqs, clus, alignerObj,
					local, weighHomopolyers, eventBased, true, ",");
			clus.expectsString = bestRef;
			//auto bestRefs = profiler::compareToRefSingle(
			//  refSeqs, clus, alignerObj, local, false, weighHomopolyers);
			//clus.expectsString = bestRefs.front();
		}
	}
	// writing
	void writeClusters(std::string filename, std::string format, bool overWrite,
                     bool exitOnFailure) const;
	void writeClusters(std::string filename, const readObjectIOOptions & ioOptions) const;

	table getReplicateInfo()const;

	double getRMSE()const;
};

}  // namspace collapse
}  // namspace bib

#ifndef NOT_HEADER_ONLY
#include "collapseCommon.cpp"
#endif
