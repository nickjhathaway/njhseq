//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//  alnInfoHolder.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/13/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "alnInfoHolder.hpp"

namespace bibseq {







void alnInfoMasterHolder::clearHolders(){
	localHolder_.clear();
	globalHolder_.clear();
}



alnInfoMasterHolder::alnInfoMasterHolder() {}

alnInfoMasterHolder::alnInfoMasterHolder(const gapScoringParameters & gapPars,
		const substituteMatrix & scoringArray) :
		localHolder_(
				std::unordered_map<std::string, alnInfoHolderBase<alnInfoLocal>> { {
						gapPars.getIdentifer(), alnInfoHolderBase<alnInfoLocal>(gapPars,
								scoringArray, "LOCAL") } }),
		globalHolder_(
				std::unordered_map<std::string, alnInfoHolderBase<alnInfoGlobal>> { {
						gapPars.getIdentifer(), alnInfoHolderBase<alnInfoGlobal>(gapPars,
								scoringArray, "GLOBAL") } }) {
}

void alnInfoMasterHolder::addHolder(const gapScoringParameters & gapPars,
  	  const substituteMatrix & scoringArray){
	localHolder_.emplace(gapPars.getIdentifer(), alnInfoHolderBase<alnInfoLocal>(gapPars,
								scoringArray, "LOCAL"));
	globalHolder_.emplace(gapPars.getIdentifer(), alnInfoHolderBase<alnInfoGlobal>(gapPars,
								scoringArray, "GLOBAL"));

}

void alnInfoMasterHolder::read(const std::string &masterDirName, bool verbose){
	auto allDirectories = getFiles(masterDirName, "", "directory", false, false);
  for (const auto &dir : allDirectories) {
  	if(verbose){
  		std::cout << "Reading from " << dir.first << std::endl;
  		std::cout << "Reading in local" << std::endl;
  	}
  	alnInfoHolderBase<alnInfoLocal> currentLocalHolder(dir.first, "LOCAL");
    if(verbose){
    	std::cout << "Reading in global" << std::endl;
    }
    alnInfoHolderBase<alnInfoGlobal> currentGlobalHolder(dir.first, "GLOBAL");
    localHolder_[currentLocalHolder.gapPars_.getIdentifer()] =
        currentLocalHolder;
    globalHolder_[currentGlobalHolder.gapPars_.getIdentifer()] =
        currentGlobalHolder;
  }
}


alnInfoMasterHolder::alnInfoMasterHolder(const std::string &masterDirName,
		const gapScoringParameters & gapPars, const substituteMatrix & scoringArray,
		bool verbose) {
	if (verbose) {
		std::cout << "Reading in Previous Alignments" << std::endl;
	}
	auto allDirectories = getFiles(masterDirName, "", "directory", false, false);
	if (allDirectories.empty()) {
		auto lHolder = alnInfoHolderBase<alnInfoLocal>(gapPars, scoringArray,
				"LOCAL");
		auto gHolder = alnInfoHolderBase<alnInfoGlobal>(gapPars, scoringArray,
				"GLOBAL");
		localHolder_[gapPars.getIdentifer()] = lHolder;
		globalHolder_[gapPars.getIdentifer()] = gHolder;
	} else {
		read(masterDirName, verbose);
	}
}


void alnInfoMasterHolder::write(const std::string &masterDirName, bool verbose ) {
  // write out the alignment infos
	if(!bib::files::bfs::exists(masterDirName)){
		bib::files::makeDir(bib::files::MkdirPar(masterDirName));
	}
  for (const auto &holder : localHolder_) {
  	if (verbose) {
  		std::cout << "Writing: " << holder.first << std::endl;
  	}
    holder.second.writeOutInfos(masterDirName, "alnInfo_" + holder.second.gapPars_.getIdentifer());
  }
  for (const auto &holder : globalHolder_) {
  	if (verbose) {
  		std::cout << "Writing: " << holder.first << std::endl;
  	}
    holder.second.writeOutInfos(masterDirName, "alnInfo_" + holder.second.gapPars_.getIdentifer());
  }
}

void alnInfoMasterHolder::mergeOtherHolder(const alnInfoMasterHolder & otherHolder){
	for(const auto & gH : otherHolder.globalHolder_){
		auto check = globalHolder_.find(gH.first);
		if(check == globalHolder_.end()){
			globalHolder_[gH.first] = gH.second;
		}else{
			check->second.mergeOtherHolder(gH.second);
		}
	}
	for(const auto & lH : otherHolder.localHolder_){
		auto check = localHolder_.find(lH.first);
		if(check == localHolder_.end()){
			localHolder_[lH.first] = lH.second;
		}else{
			check->second.mergeOtherHolder(lH.second);
		}
	}
}
}//bib
