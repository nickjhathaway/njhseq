//
//  populationCollapse.cpp
//
//  Created by Nicholas Hathaway on 1/2/14.
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

//AvgFracFoundAt

#include "populationCollapse.hpp"
namespace njhseq {
namespace collapse {


populationCollapse::populationCollapse(const std::string &populationName) :
		populationName_(populationName) {
	if (njh::containsSubString(populationName_, ".")) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ":Error populationName_ can't contain '.', "
				<< populationName_ << "\n";
		throw std::runtime_error { ss.str() };
	}
}

populationCollapse::populationCollapse(
		const std::vector<sampleCluster> &inputClusters,
		const std::string &populationName) :
		input_(clusterSet(inputClusters)),
		populationName_(populationName) {
	if (njh::containsSubString(populationName_, ".")) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error populationName_ can't contain '.', "
				<< populationName_ << "\n";
		throw std::runtime_error { ss.str() };
	}
	//no longer doing this cause for large sample numbers it would require way too much ram, e.g. 8000 samples would do 8000*8000 infos which adds up
//	for (auto & i : input_.clusters_) {
//		i.setSampInfosTotals(input_.info_.infos_);
//	}
	for (auto & i : input_.clusters_) {
		i.totalRepNumberInAnalysis_ = input_.info_.infos_.size();
	}

	for(auto & clus : input_.clusters_){
		clus.updateSampInfosFracs();
	}
//	std::cout <<njh::bashCT::boldRed("Sleeping......") << std::endl;;
//	using namespace std::chrono_literals;
//	std::this_thread::sleep_for(100000s);
}

void populationCollapse::addInput(
		const std::vector<sampleCluster> &inputClusters) {
	input_ = clusterSet(inputClusters);

	for (auto & i : input_.clusters_) {
		i.setSampInfosTotals(input_.info_.infos_);
	}
	for(auto & clus : input_.clusters_){
		clus.updateSampInfosFracs();
	}
}

void populationCollapse::popCluster(const collapser &collapserObj,
		CollapseIterations iteratorMap, const std::string &sortBy,
		aligner &alignerObj) {

	collapsed_.clusters_ = collapserObj.runClustering(input_.clusters_, iteratorMap, alignerObj);
	renameClusters(); //before with the update right afterwards, the base name will no longer represent the underlying samples
	for(auto & clus : collapsed_.clusters_){
		clus.updateSampInfosFracs();
	}

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
		const SeqIOOptions & ioOptions) const {
	collapsed_.writeClusters(outDirectory + populationName_, ioOptions);
}

uint32_t populationCollapse::numOfSamps() const {
	return collapsed_.info_.infos_.size();
}

void populationCollapse::renameClusters() {
	//need to consider different ways of sort reads before renaming the clustering
  njh::sort(collapsed_.clusters_,[](const sampleCluster & clus1, const sampleCluster & clus2){
  	if(clus1.numberOfRuns() == clus2.numberOfRuns()){
  		return clus1.getCumulativeFrac() >clus2.getCumulativeFrac();
  	}else{
  		return clus1.numberOfRuns() >clus2.numberOfRuns();
  	}
  });

  renameReadNames(collapsed_.clusters_, populationName_, true, false,false);
  for (auto &clus : collapsed_.clusters_) {
    if (clus.isClusterAtLeastHalfChimeric()) {
      clus.seqBase_.markAsChimeric();
    }
  }
  readVec::allUpdateName(collapsed_.clusters_);
}

void populationCollapse::writeFinalInitial(const std::string &outDirectory, const SeqIOOptions & ioOptions)
    const {
  for (const auto &clus : collapsed_.clusters_) {
    clus.writeClustersInDir(outDirectory, ioOptions);
  }
}

void populationCollapse::renameToOtherPopNames(const std::vector<readObject> &previousPop,
		aligner & alignerObj, comparison allowableErrors){

  VecStr alreadyTakenNames;
  std::vector<uint32_t> alreadyTakenIds;
  std::map<std::string, std::vector<uint32_t>> previousIdsByExp;
  for(const auto & popSeq : previousPop){
  	auto firstPer = popSeq.seqBase_.name_.find(".");
  	auto firstUnder = popSeq.seqBase_.name_.find("_", firstPer);
  	if(std::string::npos != firstUnder){
    	alreadyTakenNames.emplace_back(popSeq.seqBase_.name_.substr(0, firstUnder));
    	auto prevId = estd::stou(popSeq.seqBase_.name_.substr(firstPer + 1, firstUnder - 1 - firstPer));
    	auto prevExpName = popSeq.seqBase_.name_.substr(0, firstPer);
    	alreadyTakenNames.emplace_back(prevExpName);
    	alreadyTakenIds.emplace_back(prevId);
    	previousIdsByExp[prevExpName].emplace_back(prevId);
  	}else if(std::string::npos == firstUnder && std::string::npos == firstPer){
  		alreadyTakenNames.emplace_back(popSeq.seqBase_.name_);
  	}else{
    	alreadyTakenNames.emplace_back(popSeq.seqBase_.name_);
    	auto prevId = estd::stou(popSeq.seqBase_.name_.substr(firstPer + 1));
    	auto prevExpName = popSeq.seqBase_.name_.substr(0, firstPer);
    	alreadyTakenNames.emplace_back(prevExpName);
    	alreadyTakenIds.emplace_back(prevId);
    	previousIdsByExp[prevExpName].emplace_back(prevId);
  	}
  }

	uint32_t id = 0;
	if (njh::in(populationName_, previousIdsByExp)) {
		id = (*std::max_element(previousIdsByExp[populationName_].begin(),
				previousIdsByExp[populationName_].end())) + 1;
	}
	std::unordered_map<std::string, uint32_t> allIds;
	auto names = readVec::getNames(previousPop);
	for (const auto & name : names) {
		allIds[name] = 0;
	}
  for (auto &clus : collapsed_.clusters_) {
  	bool chimeric = clus.seqBase_.name_.find("CHI") != std::string::npos;
    double bestScore = 0;
    uint32_t bestRefPos = std::numeric_limits<uint32_t>::max();
    for (const auto refPos : iter::range(previousPop.size())) {
      alignerObj.alignCache(previousPop[refPos], clus, false);
      alignerObj.profilePrimerAlignment(previousPop[refPos], clus);
      if(allowableErrors.passErrorProfile(alignerObj.comp_) && alignerObj.parts_.score_ > bestScore) {
        bestScore = alignerObj.parts_.score_;
        bestRefPos = refPos;
      }
    }
    if(bestRefPos != std::numeric_limits<uint32_t>::max()){
      clus.seqBase_.name_ = previousPop[bestRefPos].seqBase_.name_ + "." + leftPadNumStr<size_t>(allIds[previousPop[bestRefPos].seqBase_.name_], collapsed_.clusters_.size());
      ++allIds[previousPop[bestRefPos].seqBase_.name_];
      clus.updateName();
    }else{
      clus.seqBase_.name_ = populationName_ + "." + leftPadNumStr<size_t>(id, collapsed_.clusters_.size());
      ++id;
      clus.updateName();
    }
    if(chimeric){
    	clus.seqBase_.markAsChimeric();
    }
  }
  //fix id name
}

void populationCollapse::renameToOtherPopNames(
    const std::vector<readObject> &previousPop,
		comparison allowableErrors) {
  // set up scoring so gaps don't cost as much and end gaps don't cost anything
  gapScoringParameters refGapScore(3.0, 1, 0.0, 0.0, 0.0, 0.0);

  uint64_t maxLen = 0;
  readVec::getMaxLength(previousPop, maxLen);
  readVec::getMaxLength(collapsed_.clusters_, maxLen);
  aligner alignerObj(maxLen, refGapScore, substituteMatrix(2,-2));
  alignerObj.weighHomopolymers_ = true;
  renameToOtherPopNames(previousPop, alignerObj, allowableErrors);

}


void populationCollapse::addRefMetaToName(const std::vector<readObject> &previousPop,
		comparison allowableErrors){
  // set up scoring so gaps don't cost as much and end gaps don't cost anything
  gapScoringParameters refGapScore(3.0, 1, 0.0, 0.0, 0.0, 0.0);

  uint64_t maxLen = 0;
  readVec::getMaxLength(previousPop, maxLen);
  readVec::getMaxLength(collapsed_.clusters_, maxLen);
  aligner alignerObj(maxLen, refGapScore, substituteMatrix(2,-2));
  alignerObj.weighHomopolymers_ = true;
  addRefMetaToName(previousPop, alignerObj, allowableErrors);

}
void populationCollapse::addRefMetaToName(
		const std::vector<readObject> &previousPop, aligner & alignerObj,
		comparison allowableErrors) {
	VecStr alreadyTakenNames;
	std::vector<uint32_t> alreadyTakenIds;
	std::map<std::string, std::vector<uint32_t>> previousIdsByExp;
	for (const auto & popSeq : previousPop) {
		auto firstPer = popSeq.seqBase_.name_.find(".");
		auto firstUnder = popSeq.seqBase_.name_.find("_", firstPer);
		if (std::string::npos != firstUnder) {
			alreadyTakenNames.emplace_back(
					popSeq.seqBase_.name_.substr(0, firstUnder));
			auto prevId = estd::stou(
					popSeq.seqBase_.name_.substr(firstPer + 1,
							firstUnder - 1 - firstPer));
			auto prevExpName = popSeq.seqBase_.name_.substr(0, firstPer);
			alreadyTakenNames.emplace_back(prevExpName);
			alreadyTakenIds.emplace_back(prevId);
			previousIdsByExp[prevExpName].emplace_back(prevId);
		} else if (std::string::npos == firstUnder
				&& std::string::npos == firstPer) {
			alreadyTakenNames.emplace_back(popSeq.seqBase_.name_);
		} else {
			alreadyTakenNames.emplace_back(popSeq.seqBase_.name_);
			auto prevId = estd::stou(popSeq.seqBase_.name_.substr(firstPer + 1));
			auto prevExpName = popSeq.seqBase_.name_.substr(0, firstPer);
			alreadyTakenNames.emplace_back(prevExpName);
			alreadyTakenIds.emplace_back(prevId);
			previousIdsByExp[prevExpName].emplace_back(prevId);
		}
	}


	std::unordered_map<std::string, uint32_t> allIds;
	auto names = readVec::getNames(previousPop);
	for (const auto & name : names) {
		allIds[name] = 0;
	}
	for (auto &clus : collapsed_.clusters_) {
		double bestScore = 0;
		uint32_t bestRefPos = std::numeric_limits<uint32_t>::max();
		for (const auto refPos : iter::range(previousPop.size())) {
			alignerObj.alignCache(previousPop[refPos], clus, false);
			alignerObj.profilePrimerAlignment(previousPop[refPos], clus);
			if (allowableErrors.passErrorProfile(alignerObj.comp_)
					&& alignerObj.parts_.score_ > bestScore) {
				bestScore = alignerObj.parts_.score_;
				bestRefPos = refPos;
			}
		}
		if (bestRefPos != std::numeric_limits<uint32_t>::max()) {
			MetaDataInName meta;
			if (MetaDataInName::nameHasMetaData(clus.seqBase_.name_)) {
				meta.processNameForMeta(clus.seqBase_.name_, true);
			}
			meta.addMeta("ref",
					previousPop[bestRefPos].seqBase_.name_ + "."
							+ leftPadNumStr<size_t>(
									allIds[previousPop[bestRefPos].seqBase_.name_],
									collapsed_.clusters_.size()));
			clus.seqBase_.resetMetaInName(meta);
			++allIds[previousPop[bestRefPos].seqBase_.name_];
			clus.updateName();
		}
	}
}

void populationCollapse::updateInfoWithSampCollapses(const std::map<std::string, sampleCollapse> & sampCollapses){
	for(const auto & samp : sampCollapses){
		updateInfoWithSampCollapse(samp.second);
	}
}
void populationCollapse::updateInfoWithSampCollapse(const sampleCollapse & sampCollapse){
	//only update with non-zero collapses
	if(sampCollapse.collapsed_.clusters_.size() >0){
		collapsed_.info_.updateCoi(sampCollapse.collapsed_.clusters_.size());
	}
}


VecStr populationCollapse::getPopInfoVec() const {
	double he = getExpectedHeterozygosity();
	return toVecStr(populationName_, collapsed_.info_.totalReadCount_,
			collapsed_.info_.numberOfClusters_, collapsed_.info_.infos_.size(),
			collapsed_.clusters_.size(), he, collapsed_.info_.getCoiInfoVec());
}

double populationCollapse::getExpectedHeterozygosity() const{
	//he = [𝑛/(𝑛 − 1)][1 − ∑𝑛𝑖=1 𝑝𝑖2]
	double sumOfSquares = 0;
	for(const auto & clus : collapsed_.clusters_){
		sumOfSquares += std::pow(clus.reads_.size()/static_cast<double>(collapsed_.info_.numberOfClusters_), 2);
	}
	return (collapsed_.info_.numberOfClusters_/(static_cast<double>(collapsed_.info_.numberOfClusters_) - 1)) * (1 - sumOfSquares);
}

VecStr populationCollapse::getPopInfoHeaderVec() {
	return toVecStr("p_Name","p_TotalInputReadCnt", "p_TotalInputClusterCnt",
			"p_TotalPopulationSampCnt", "p_TotalHaplotypes","p_ExpectedHeterozygosity",
			clusterSetInfo::getCoiInfoHeaderVec());
}


}  // namespace collapse
}  // namespace njhseq
