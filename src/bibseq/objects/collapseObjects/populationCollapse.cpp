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
  for (auto &clus : collapsed_.clusters_) {
    double bestScore = 0;
    uint32_t bestRefPos = 0;
    for (const auto &refPos : iter::range(previousPop.size())) {
      alignerObj.alignVec(previousPop[refPos], clus, false);
      if (alignerObj.parts_.score_ > bestScore) {
        bestScore = alignerObj.parts_.score_;
        bestRefPos = refPos;
      }
    }
    clus.seqBase_.name_ = previousPop[bestRefPos].seqBase_.name_;
    clus.updateName();
  }
  // set the scoring back to what it was
  //alignerObj.setGapScoring(currentGapScores);
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
