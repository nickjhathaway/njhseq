
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
               const std::string &sampName, uint32_t freqCutOff): sampName_(sampName) {
  uint32_t lowFreqCount = 0;
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

  for(const auto & seq : sampleClusterReads){
  	if(seq.seqBase_.cnt_ > freqCutOff){
  		collapsed_.clusters_.emplace_back(seq);
  	}else{
  		excluded_.clusters_.emplace_back(seq);
  		++lowFreqCount;
  	}
  }
  //161104-miseq-D10-GHA-4
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


void sampleCollapse::cluster(const collapser &collapserObj,
		CollapseIterations iteratorMap, const std::string &sortBy,
		aligner &alignerObj) {
	// std::cout <<"clus 1 " << std::endl;
	collapsed_.clusters_ = collapserObj.runClustering(collapsed_.clusters_,
			iteratorMap, alignerObj);
	// std::cout <<"clus 2 " << std::endl;
}



// excludes
void sampleCollapse::excludeChimeras(bool update) {
  for (auto &clus : collapsed_.clusters_) {
    if (clus.isClusterAtLeastHalfChimeric()) {
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
    if (clus.isClusterAtLeastChimericCutOff(fracCutOff)) {
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

void sampleCollapse::excludeBySampNum(uint32_t sampsRequired, bool update) {
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
    if (clus.isClusterAtLeastHalfChimeric()) {
      clus.seqBase_.markAsChimeric();
      //clus.remove = true;
    }
  }
  readVec::allUpdateName(collapsed_.clusters_);
}

// output
std::string sampleCollapse::getSimpleSampInfoHeader(const std::string & delim) {
	return bib::conToStr(getSimpleSampInfoHeaderVec(), delim);
}

std::string sampleCollapse::getSimpleSampInfo(const std::string &delim) const {
	// sampName\treadsUsed(%total)\thapsUsed
	return bib::conToStr(getSimpleSampInfoVec(), delim);
}

VecStr sampleCollapse::getSimpleSampInfoVec() const {
	return toVecStr(sampName_,
			getPercentageString(collapsed_.info_.totalReadCount_,
					input_.info_.totalReadCount_), collapsed_.info_.numberOfClusters_,
			collapsed_.clusters_.size());
}

VecStr sampleCollapse::getSimpleSampInfoHeaderVec() {
	return VecStr { "s_Name", "s_ReadCntTotUsed", "s_InputClusterCnt",
			"s_FinalClusterCnt" };
	//
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
		currentInfo << getSimpleSampInfo(delim);
		currentInfo << delim << pos;
		currentInfo << delim << clus.getClusterInfo(delim);
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
      subRead->seqBase_.name_ = out.seqBase_.name_;
    }
  }
  return output;
}

// writing
void sampleCollapse::writeExcluded(const std::string &outDirectory,
		const SeqIOOptions & ioOptions) const {
  excluded_.writeClusters(outDirectory + sampName_, ioOptions);
}

void sampleCollapse::writeExcludedOriginalClusters(const std::string &outDirectory,
		const SeqIOOptions & ioOptions) const {
  for (const auto &clus : excluded_.clusters_) {
    clus.writeClustersInDir(outDirectory, ioOptions);
  }
}

void sampleCollapse::writeInitial(const std::string &outDirectory,
		const SeqIOOptions & ioOptions) const {
  input_.writeClusters(outDirectory + sampName_, ioOptions);
}

void sampleCollapse::writeFinal(const std::string &outDirectory,
		const SeqIOOptions & ioOptions) const {
  collapsed_.writeClusters(outDirectory + sampName_, ioOptions);
}

void sampleCollapse::writeFinalOrignalClusters(const std::string &outDirectory,
		const SeqIOOptions & ioOptions) const {
	for (const auto &clus : collapsed_.clusters_) {
		clus.writeClustersInDir(outDirectory, ioOptions);
	}
}



}  // napsace collapse
}  // namespace bib
