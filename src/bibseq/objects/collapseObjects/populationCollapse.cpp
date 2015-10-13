//
//  populationCollapse.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "populationCollapse.hpp"
namespace bibseq {
namespace collapse {


populationCollapse::populationCollapse(){

}
populationCollapse::populationCollapse(const std::vector<sampleCluster> &inputClusters,
                   const std::string &populationName)
    : input_(clusterSet(inputClusters)), populationName_(populationName) {
	if(containsSubString(populationName_, ".")){
		throw std::runtime_error{"Error in populationCollapse::populationCollapse, populationName can't contain '.', " + populationName_};
	}
	for(auto & i : input_.clusters_){
		i.setSampInfos(input_.info_.infos_);
	}
}

void populationCollapse::popCluster(
    collapser &collapserObj, std::map<int, std::vector<double>> iteratorMap,
    const std::string &sortBy, aligner &alignerObj) {

  collapsed_.clusters_ = collapserObj.collapseCluster(
      input_.clusters_, iteratorMap, sortBy, alignerObj);
  renameClusters();//before with the update right afterwards, the base name will no longer represent the underlying samples
  updateCollapsedInfos();
}

void populationCollapse::popClusterOnId(collapser &collapserObj,
                  std::map<int, std::vector<double>> iteratorMap,
                  const std::string &sortBy, aligner &alignerObj){
  collapsed_.clusters_ = collapserObj.collapseClusterOnPerId(
      input_.clusters_, iteratorMap, sortBy, alignerObj);
  renameClusters(); //before with the update right afterwards, the base name will no longer represent the underlying samples
  updateCollapsedInfos();
}

// update the initial infos
void populationCollapse::updateInitialInfos() {
  input_.setSetInfo();
}

// update the collapsed infos
void populationCollapse::updateCollapsedInfos() {
  collapsed_.setSetInfo();
  collapsed_.setSubClustersPositions();
}


void populationCollapse::writeFinal(const std::string &outDirectory,
                                    const std::string &outFormat,
                                    bool overWrite,
                                    bool exitOnFailureToWrite) const {
  collapsed_.writeClusters(outDirectory + populationName_, outFormat, overWrite,
                           exitOnFailureToWrite);
}

uint32_t populationCollapse::numOfSamps()const{
	return collapsed_.info_.infos_.size();
}

void populationCollapse::renameClusters() {
	//need to consider different ways of sort reads before renaming the clustering
  bib::sort(collapsed_.clusters_,[](const sampleCluster & clus1, const sampleCluster & clus2){
  	if(clus1.numberOfRuns() == clus2.numberOfRuns()){
  		return clus1.getCumulativeFrac() >clus2.getCumulativeFrac();
  	}else{
  		return clus1.numberOfRuns() >clus2.numberOfRuns();
  	}
  });

  renameReadNames(collapsed_.clusters_, populationName_, true, false,false);
  for (auto &clus : collapsed_.clusters_) {
    if (clusterVec::isClusterAtLeastHalfChimeric(clus.reads_)) {
      clus.seqBase_.markAsChimeric();
    }
  }
  readVec::allUpdateName(collapsed_.clusters_);
}

void populationCollapse::writeFinalInitial(const std::string &outDirectory, const readObjectIOOptions & ioOptions)
    const {
  for (const auto &clus : collapsed_.clusters_) {
    clus.writeOutClusters(outDirectory, ioOptions);
  }
}

void populationCollapse::renameToOtherPopNames(
    const std::vector<readObject> &previousPop, aligner &alignerObjIn) {
  // set up scoreing so gaps don't cost as much and end gaps don't cost anything
  gapScoringParameters refGapScore(3.0, 1, 0.0, 0.0, 0.0, 0.0);

  uint64_t maxLen = 0;
  readVec::getMaxLength(previousPop, maxLen);
  readVec::getMaxLength(collapsed_.clusters_, maxLen);
  aligner alignerObj(maxLen, refGapScore, substituteMatrix(2,-2));
  comparison noErrors;

  VecStr alreadyTakenNames;
  std::vector<uint32_t> alreadyTakenIds;
  std::map<std::string, std::vector<uint32_t>> previousIdsByExp;
  for(const auto & read : previousPop){
  	auto firstPer = read.seqBase_.name_.find(".");
  	auto firstUnder = read.seqBase_.name_.find("_", firstPer);
  	alreadyTakenNames.emplace_back(read.seqBase_.name_.substr(0, firstUnder));
  	auto prevId = bib::lexical_cast<uint32_t>(read.seqBase_.name_.substr(firstPer + 1, firstUnder -1 - firstPer));
  	auto prevExpName = read.seqBase_.name_.substr(0, firstPer);
  	alreadyTakenNames.emplace_back(prevExpName);
  	alreadyTakenIds.emplace_back(prevId);
  	previousIdsByExp[prevExpName].emplace_back(prevId);
  }
  uint32_t id = 0;
  if(bib::in(populationName_, previousIdsByExp)){
  	id = (*std::max_element(previousIdsByExp[populationName_].begin(), previousIdsByExp[populationName_].end())) + 1;
  }
  for (auto &clus : collapsed_.clusters_) {
  	bool chimeric = clus.seqBase_.name_.find("CHI") != std::string::npos;
    double bestScore = 0;
    uint32_t bestRefPos = std::numeric_limits<uint32_t>::max();
    for (const auto &refPos : iter::range(previousPop.size())) {
      alignerObj.alignVec(previousPop[refPos], clus, false);
      alignerObj.profilePrimerAlignment(previousPop[refPos], clus,false);
      if(noErrors.passErrorProfile(alignerObj.comp_) &&  alignerObj.parts_.score_ > bestScore) {
        bestScore = alignerObj.parts_.score_;
        bestRefPos = refPos;
      }
    }
    if(bestRefPos != std::numeric_limits<uint32_t>::max()){
      clus.seqBase_.name_ = previousPop[bestRefPos].seqBase_.name_;
      clus.updateName();
    }else{
      clus.seqBase_.name_ = populationName_ + "." + estd::to_string(id);
      ++id;
      clus.updateName();
    }
    if(chimeric){
    	clus.seqBase_.markAsChimeric();
    }
  }
  //fix id name
}


void populationCollapse::updateInfoWithSampCollapses(const std::map<std::string, sampleCollapse> & sampCollapses){
	for(const auto & samp : sampCollapses){
		updateInfoWithSampCollapse(samp.second);
	}
}
void populationCollapse::updateInfoWithSampCollapse(const sampleCollapse & sampCollapse){
	//only update with non-zero collapses
	if(sampCollapse.collapsed_.clusters_.size() >0){
		collapsed_.info_.updateMoi(sampCollapse.collapsed_.clusters_.size());
	}
}
}  // namespace collapse
}  // namespace bibseq
