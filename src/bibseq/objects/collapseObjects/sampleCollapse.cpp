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
//  sampleCollapse.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/31/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "sampleCollapse.hpp"

namespace bibseq {
namespace collapse {



sampleCollapse::sampleCollapse(const std::vector<std::vector<bibseq::cluster>> &inputClusters,
               const std::string &sampName, uint32_t sizeCutOff): sampName_(sampName) {
  uint32_t singletCount = 0;
  std::map<std::string, double> allSampCounts;
  for (auto &reads : inputClusters) {
  	//clusterVec::allSetFractionClusters(reads);
    for(const auto & read : reads){
    	allSampCounts[read.getOwnSampName()] += read.seqBase_.cnt_;
    }
  }
  std::map<std::string, sampInfo> allSampInfos;
  for(const auto & sCount : allSampCounts){
  	allSampInfos[sCount.first] = sampInfo(sCount.first, sCount.second);
  }

	std::vector<sampleCluster> sampleClusterReads;

  for (const auto &reads : inputClusters) {
  	for(const auto & read : reads){
  		sampleClusterReads.emplace_back(read,allSampInfos);
  	}
  }

  addOtherVec(input_.clusters_, sampleClusterReads);

  for(const auto & read : sampleClusterReads){
  	if(read.seqBase_.cnt_ > sizeCutOff){
  		collapsed_.clusters_.emplace_back(read);
  	}else{
  		excluded_.clusters_.emplace_back(read);
  		++singletCount;
  	}
  }
  updateInitialInfos();
  updateCollapsedInfos();
}

void sampleCollapse::updateExclusionInfos() {
  excluded_.setSetInfo();
  sampleCluster::updateAllClusters(excluded_.clusters_, input_.info_.infos_);
}
// update the initial infos
void sampleCollapse::updateInitialInfos() {
  input_.setSetInfo();
}
// update the collapsed infos
void sampleCollapse::updateCollapsedInfos() {
  collapsed_.setSetInfo();

  sampleCluster::updateAllClusters(collapsed_.clusters_, collapsed_.info_.infos_);
}

void sampleCollapse::updateAfterExclustion() {
  updateExclusionInfos();
  updateCollapsedInfos();
}


void sampleCollapse::cluster(collapser &collapserObj,
                             std::map<int, std::vector<double>> iteratorMap,
                             const std::string &sortBy, aligner &alignerObj) {
  // std::cout <<"clus 1 " << std::endl;
  collapsed_.clusters_ = collapserObj.collapseCluster(
      collapsed_.clusters_, iteratorMap, sortBy, alignerObj);
  // std::cout <<"clus 2 " << std::endl;
}

void sampleCollapse::clusterOnPerId(collapser &collapserObj,
             std::map<int, std::vector<double>> iteratorMap,
             const std::string &sortBy, aligner &alignerObj){
  // std::cout <<"clus 1 " << std::endl;
  collapsed_.clusters_ = collapserObj.collapseClusterOnPerId(
      collapsed_.clusters_, iteratorMap, sortBy, alignerObj);
  // std::cout <<"clus 2 " << std::endl;
}

// excludes
void sampleCollapse::excludeChimeras(bool update) {
  for (auto &clus : collapsed_.clusters_) {
    if (clusterVec::isClusterAtLeastHalfChimeric(clus.reads_)) {
      clus.seqBase_.markAsChimeric();
      clus.remove = true;
    }
  }
  uint32_t chimeraNum = 0;
  //collapsed_.clusters_ = readVecSplitter::splitVectorOnRemoveAdd(
  //    collapsed_.clusters_, excluded_.clusters_, chimeraNum, "none", false);
  collapsed_.clusters_ = readVecSplitter::splitVectorWithNameContainingAdd(
        collapsed_.clusters_,"CHI_", excluded_.clusters_, chimeraNum, false);
  if (update) {
    updateAfterExclustion();
  }
}

// excludes
void sampleCollapse::excludeChimeras(bool update, double fracCutOff) {
  for (auto &clus : collapsed_.clusters_) {
    if (clusterVec::isClusterAtLeastChimericCutOff(clus.reads_, fracCutOff)) {
      clus.seqBase_.markAsChimeric();
      clus.remove = true;
    }
  }
  uint32_t chimeraNum = 0;
  collapsed_.clusters_ = readVecSplitter::splitVectorOnRemoveAdd(
      collapsed_.clusters_, excluded_.clusters_, chimeraNum, "none", false);
  //collapsed_.clusters_ = readVecSplitter::splitVectorWithNameContainingAdd(
  //      collapsed_.clusters_,"CHI_", excluded_.clusters_, chimeraNum, false);
  if (update) {
    updateAfterExclustion();
  }
}


void sampleCollapse::excludeFraction(double fractionCutOff, bool update) {
  uint32_t fractionCutOffNum = 0;
  collapsed_.clusters_ = readVecSplitter::splitVectorOnReadFractionAdd(
      collapsed_.clusters_, fractionCutOff, excluded_.clusters_,
      fractionCutOffNum, false);
  if (update) {
    updateAfterExclustion();
  }
}

void sampleCollapse::excludeBySampNum(uint64_t sampsRequired, bool update) {
  uint32_t splitCount = 0;
  for (auto &clus : collapsed_.clusters_) {
    if (clus.numberOfRuns() < sampsRequired) {
      clus.remove = true;
    }
  }
  collapsed_.clusters_ = readVecSplitter::splitVectorOnRemoveAdd(
      collapsed_.clusters_, excluded_.clusters_, splitCount, "none", false);
  if (update) {
    updateAfterExclustion();
  }
}
void sampleCollapse::renameClusters(const std::string &sortBy) {
  readVecSorter::sortReadVector(collapsed_.clusters_, sortBy);
  renameReadNames(collapsed_.clusters_, sampName_, true, false);
  for (auto &clus : collapsed_.clusters_) {
    if (clusterVec::isClusterAtLeastHalfChimeric(clus.reads_)) {
      clus.seqBase_.markAsChimeric();
      //clus.remove = true;
    }
  }
  readVec::allUpdateName(collapsed_.clusters_);
}

// output
std::string sampleCollapse::getSimpleSampInfoHeader(const std::string & delim){
	std::stringstream ss;
	ss << "s_Name"
				<< delim << "s_ReadCntTotUsed"
				<< delim << "s_InputClusterCnt"
				<< delim << "s_FinalClusterCnt"
				<< delim << "c_clusterID";
	return ss.str();
}
std::string sampleCollapse::getSimpleSampInfo(uint32_t clusterId, const std::string &delim) const {
  // sampName\treadsUsed(%total)\thapsUsed
	std::stringstream ss;
	ss << sampName_
				<< delim << getPercentageString(collapsed_.info_.totalReadCount_,
                                                  input_.info_.totalReadCount_)
				<< delim << collapsed_.info_.numberOfClusters_
				<< delim << collapsed_.clusters_.size()
				<< delim << clusterId;
	return ss.str();
}



VecStr sampleCollapse::getAllInfoVec(bool checkingExpected,
                                     const std::string &delim) const {
  auto mapOfInfos = getAllInfoMap(checkingExpected, delim);
  VecStr ans = getVectorOfMapValues(mapOfInfos);
  return ans;
}

std::map<std::string, std::string, std::greater<std::string>> sampleCollapse::getAllInfoMap(
		bool checkingExpected, const std::string &delim,
		uint32_t maxRepCount) const {
	std::map<std::string, std::string, std::greater<std::string>> ans;
	uint32_t pos = 0;
	for (const auto &clus : collapsed_.clusters_) {
		std::stringstream currentInfo;
		currentInfo << getSimpleSampInfo(pos, delim);
		currentInfo << delim << clus.getSampleInfo(delim);
		currentInfo << delim
				<< clus.getRepsInfo(input_.info_.infos_, excluded_.info_.infos_,
						collapsed_.info_.infos_, maxRepCount, checkingExpected, delim);
		ans.insert( { clus.getStubName(true), currentInfo.str() });
		++pos;
	}
	return ans;
}

std::vector<sampleCluster> sampleCollapse::createOutput(
    bool renameFirst, const std::string &sortBy) {
  if (renameFirst) {
    renameClusters(sortBy);
  }
  std::vector<sampleCluster> output;
  double totalReadCnt_ = 0;
  for (const auto &out : collapsed_.clusters_) {
  	totalReadCnt_ += out.seqBase_.cnt_;
  }
  std::map<std::string, sampInfo> outSampInfos{{sampName_, sampInfo(sampName_, totalReadCnt_)}};
  for (const auto &out : collapsed_.clusters_) {
    output.emplace_back(sampleCluster(out.createRead(), outSampInfos));
  }
  readVec::allUpdateName(output);
  // need to change subreads names so look up can happen correctly
  for (auto &out : output) {
    for (auto &subRead : out.reads_) {
      subRead.seqBase_.name_ = out.seqBase_.name_;
    }
  }
  return output;
}

// writing
void sampleCollapse::writeExcluded(const std::string &outDirectory,
		const readObjectIOOptions & ioOptions) const {
  excluded_.writeClusters(outDirectory + sampName_, ioOptions);
}

void sampleCollapse::writeExcludedOriginalClusters(const std::string &outDirectory,
		const readObjectIOOptions & ioOptions) const {
  for (const auto &clus : excluded_.clusters_) {
    clus.writeOutClusters(outDirectory, ioOptions);
  }
}

void sampleCollapse::writeInitial(const std::string &outDirectory,
		const readObjectIOOptions & ioOptions) const {
  input_.writeClusters(outDirectory + sampName_, ioOptions);
}

void sampleCollapse::writeFinal(const std::string &outDirectory,
		const readObjectIOOptions & ioOptions) const {
  collapsed_.writeClusters(outDirectory + sampName_, ioOptions);
}

void sampleCollapse::writeFinalOrignalClusters(const std::string &outDirectory, const readObjectIOOptions & ioOptions) const {
  for (const auto &clus : collapsed_.clusters_) {
    clus.writeOutClusters(outDirectory, ioOptions);
  }
}



}  // napsace collapse
}  // namespace bib
