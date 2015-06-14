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
//  infoPrinter.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/8/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "infoPrinter.hpp"
namespace bibseq {

void infoPrinter::printSampleCollapseInfo(
    const std::map<std::string, collapse::sampleCollapse>& sampCollapses,
    bool checkingExpected, const std::string fileName,
    const collapse::populationCollapse& popCollapse, bool population) {
  std::ofstream infoFile;
  openTextFile(infoFile, fileName, ".txt", true, false);
  std::string delim = "\t";
  if (population) {
    infoFile << "s_Sample"
    		<< delim << sampleCluster::getPopInfoHeader(delim)
    		<< delim;
  }
  uint64_t maxRunCount = 0;
  for (auto& sampCollapse : sampCollapses) {
    if (maxRunCount < sampCollapse.second.collapsed_.info_.infos_.size()) {
      maxRunCount = sampCollapse.second.collapsed_.info_.infos_.size();
    }
  }
  infoFile << collapse::sampleCollapse::getSimpleSampInfoHeader(delim)
  					<< delim << sampleCluster::getSampleInfoHeader(delim)
  					<< delim << sampleCluster::getRepsInfoHeader(maxRunCount,checkingExpected, delim);
  infoFile << "\n";
  for (auto& sampCollapse : sampCollapses) {
    auto currentInfos = sampCollapse.second.getAllInfoMap(checkingExpected,delim, maxRunCount);
    for (const auto& info : currentInfos) {
			if (population) {
				infoFile << sampCollapse.second.sampName_;
				infoFile << delim
						<< popCollapse.collapsed_.clusters_[popCollapse.collapsed_.subClustersPositions_.at(
								info.first)].getPopInfo(
								popCollapse.collapsed_.info_.totalReadCount_,
								popCollapse.collapsed_.info_.numberOfClusters_,
								popCollapse.collapsed_.info_.infos_.size(), delim) << delim;
			}
      infoFile << info.second << "\n";
    }
  }
}

void infoPrinter::printPopulationCollapseInfo(
		const collapse::populationCollapse& popCollapse,
		const std::string& fileName, bool checkingExpected) {

	std::ofstream popInfoOutFile;
	openTextFile(popInfoOutFile, fileName, ".txt", true, false);
	std::string delim = "\t";
	popInfoOutFile << sampleCluster::getFullPopInfoHeader(collapse::clusterSetInfo::getMoiInfoHeader(delim),delim);
	if (checkingExpected) {
		popInfoOutFile << "\tbestExpected";
	}
	popInfoOutFile << "\n";
	for (const auto& clus : popCollapse.collapsed_.clusters_) {
		popInfoOutFile
				<< clus.getFullPopInfo(popCollapse.collapsed_.info_.totalReadCount_,
						popCollapse.collapsed_.info_.numberOfClusters_,
						popCollapse.collapsed_.info_.infos_.size(),
						popCollapse.collapsed_.clusters_.size(),
						popCollapse.collapsed_.info_.getMoiInfo(delim), delim);
		if (checkingExpected) {
			popInfoOutFile << "\t" << clus.expectsString;
		}
		popInfoOutFile << "\n";
	}
}

void infoPrinter::printInfoForSamps(
      const std::map<std::string, collapse::sampleCollapse>& sampCollapses,
      bool checkingExpected, std::ostream & sampOut, std::ostream & popOut,
      const collapse::populationCollapse& popCollapse, bool population, const VecStr & forSamps,const std::string & groupName,
			const std::set<uint32_t> & otherPopPositions){
  std::string delim = "\t";
  collapse::clusterSetInfo currentPopSetInfo;
  VecStr sampClusterNames;
  std::set<uint32_t> popPositions;

  for(const auto & samp : forSamps){
  	/**@todo should check to see if the sample exist*/
  	for(const auto & subClus : sampCollapses.at(samp).collapsed_.clusters_){
  		sampClusterNames.emplace_back(subClus.seqBase_.name_);
  		popPositions.emplace(popCollapse.collapsed_.subClustersPositions_.at(subClus.getStubName(true)));
  		currentPopSetInfo.updateInfo(subClus);
  	}
  	if(sampCollapses.at(samp).collapsed_.clusters_.size() > 0){
  		currentPopSetInfo.updateMoi(sampCollapses.at(samp).collapsed_.clusters_.size());
  	}
  }
  std::set<uint32_t> uniqueHaplotypes;
  std::set_difference(popPositions.begin(), popPositions.end(), otherPopPositions.begin(), otherPopPositions.end(),
                          std::inserter(uniqueHaplotypes, uniqueHaplotypes.begin()));

  if (population) {
  	sampOut << "s_Sample" << delim << "g_GroupName"
    		<< delim << sampleCluster::getPopInfoHeader(delim)
    		<< delim;
  }
  uint64_t maxRunCount = 0;
  for (auto& sampCollapse : sampCollapses) {
    if (maxRunCount < sampCollapse.second.collapsed_.info_.infos_.size()) {
      maxRunCount = sampCollapse.second.collapsed_.info_.infos_.size();
    }
  }
  sampOut << collapse::sampleCollapse::getSimpleSampInfoHeader(delim)
  					<< delim << sampleCluster::getSampleInfoHeader(delim)
  					<< delim << sampleCluster::getRepsInfoHeader(maxRunCount,checkingExpected, delim);
  sampOut << "\n";

  for (auto& samp : forSamps) {
  	const auto & sampCollapse = sampCollapses.at(samp);
    auto currentInfos = sampCollapse.getAllInfoMap(checkingExpected,delim, maxRunCount);
    for (const auto& info : currentInfos) {
			if (population) {
				sampOut << sampCollapse.sampName_;
				sampOut << delim << groupName;
				sampOut << delim
						<< popCollapse.collapsed_.clusters_[popCollapse.collapsed_.subClustersPositions_.at(
								info.first)].getPopInfo(
										currentPopSetInfo.totalReadCount_,
										currentPopSetInfo.numberOfClusters_,
										currentPopSetInfo.infos_.size(),sampClusterNames, delim) << delim;
			}
			sampOut << info.second << "\n";
    }
  }

	popOut
			<< sampleCluster::getFullPopInfoHeader(
					"p_TotalUniqueHaplotypes" + delim + collapse::clusterSetInfo::getMoiInfoHeader(delim), delim);
	if (checkingExpected) {
		popOut << delim << "bestExpected";
	}
	popOut << delim << "g_GroupName";
	popOut << delim << "g_hapsFoundOnlyInThisGroup";
	popOut << "\n";
	VecStr hapNames;
	for(const auto & uh : uniqueHaplotypes){
		hapNames.emplace_back(popCollapse.collapsed_.clusters_[uh].getStubName(false));
	}
	for (const auto & pos : popPositions) {
		//for (const auto& clus : popCollapse.collapsed_.clusters_) {
		const auto& clus = popCollapse.collapsed_.clusters_[pos];
		popOut
				<< clus.getFullPopInfo(currentPopSetInfo.totalReadCount_,
						currentPopSetInfo.numberOfClusters_,
						currentPopSetInfo.infos_.size(),
						popPositions.size(), estd::to_string(hapNames.size()) + delim + currentPopSetInfo.getMoiInfo(delim), sampClusterNames, delim);
		if (checkingExpected) {
			popOut << "\t" << clus.expectsString;
		}
		popOut << delim << groupName;
		popOut << delim << vectorToString(hapNames, ",");
		popOut << "\n";
	}
}

}  // bibseq
