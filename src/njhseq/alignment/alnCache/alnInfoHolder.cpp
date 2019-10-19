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
//
//  alnInfoHolder.hpp
//
//  Created by Nicholas Hathaway on 1/13/14.
//

#include "alnInfoHolder.hpp"

namespace njhseq {







void alnInfoMasterHolder::clearHolders(){
	localHolder_.clear();
	globalHolder_.clear();
}



alnInfoMasterHolder::alnInfoMasterHolder() {}

//alnInfoMasterHolder::alnInfoMasterHolder(const gapScoringParameters & gapPars,
//		const substituteMatrix & scoringArray) :
//		localHolder_(
//				std::unordered_map<std::string, alnInfoHolderBase<alnInfoLocal>> { {
//						gapPars.getIdentifer(), alnInfoHolderBase<alnInfoLocal>(gapPars,
//								scoringArray, "LOCAL") } }),
//		globalHolder_(
//				std::unordered_map<std::string, alnInfoHolderBase<alnInfoGlobal>> { {
//						gapPars.getIdentifer(), alnInfoHolderBase<alnInfoGlobal>(gapPars,
//								scoringArray, "GLOBAL") } }) {
//}

alnInfoMasterHolder::alnInfoMasterHolder(const gapScoringParameters & gapPars,
		const substituteMatrix & scoringArray) {
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	addHolder(gapPars, scoringArray);
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
}

void alnInfoMasterHolder::addHolder(const gapScoringParameters & gapPars,
  	  const substituteMatrix & scoringArray){
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	const auto & scoreCopy = scoringArray;

	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::array<std::array<int32_t, 127>, 127> matCopy = scoringArray.mat_;

	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::array<std::array<int32_t, 127>, 127> matCopyByHand;
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(const auto i : iter::range(scoringArray.mat_.size())){
		for(const auto & j : iter::range(scoringArray.mat_[i].size())){
			matCopyByHand[i][j] = scoringArray.mat_[i][j];
		}
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::shared_ptr<substituteMatrix> mat = std::make_shared<substituteMatrix>();
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::shared_ptr<substituteMatrix> mat2 = std::make_shared<substituteMatrix>(scoringArray);
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//substituteMatrix scores(matCopy);
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	auto subHolder = alnInfoHolderBase<alnInfoLocal>(gapPars,
//			scoringArray, "LOCAL");
	const substituteMatrix & matRef = *mat2;
	std::cout << matRef.mat_[89][89] << std::endl;
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	localHolder_.emplace(gapPars.getIdentifer(), alnInfoHolderBase<alnInfoLocal>(gapPars,
//								*mat2, "LOCAL"));
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	globalHolder_.emplace(gapPars.getIdentifer(), alnInfoHolderBase<alnInfoGlobal>(gapPars,
//								scoringArray, "GLOBAL"));
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
}

void alnInfoMasterHolder::addHolder(const gapScoringParameters & gapPars,
  	  const substituteMatrix & scoringArray){
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	const auto & scoreCopy = scoringArray;

	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::array<std::array<int32_t, 127>, 127> matCopy = scoringArray.mat_;

	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::array<std::array<int32_t, 127>, 127> matCopyByHand;
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(const auto i : iter::range(scoringArray.mat_.size())){
		for(const auto & j : iter::range(scoringArray.mat_[i].size())){
			matCopyByHand[i][j] = scoringArray.mat_[i][j];
		}
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::shared_ptr<substituteMatrix> mat = std::make_shared<substituteMatrix>();
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::shared_ptr<substituteMatrix> mat2 = std::make_shared<substituteMatrix>(scoringArray);
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//substituteMatrix scores(matCopy);
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	auto subHolder = alnInfoHolderBase<alnInfoLocal>(gapPars,
//			scoringArray, "LOCAL");
	const substituteMatrix & matRef = *mat2;
	std::cout << matRef.mat_[89][89] << std::endl;
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	localHolder_.emplace(gapPars.getIdentifer(), alnInfoHolderBase<alnInfoLocal>(gapPars,
//								*mat2, "LOCAL"));
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	globalHolder_.emplace(gapPars.getIdentifer(), alnInfoHolderBase<alnInfoGlobal>(gapPars,
//								scoringArray, "GLOBAL"));
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
}

void alnInfoMasterHolder::read(const std::string &masterDirName, bool verbose) {
	auto allDirectories = getFiles(masterDirName, "", "directory", false, false);
	for (const auto &dir : allDirectories) {
		if (verbose) {
			std::cout << "Reading from " << dir.first << std::endl;
			std::cout << "Reading in local" << std::endl;
		}
		alnInfoHolderBase<alnInfoLocal> currentLocalHolder(dir.first, "LOCAL");
		if (verbose) {
			std::cout << "Read in  " << currentLocalHolder.infos_.size() << " LOCAL alignments" << std::endl;
		}
		if (verbose) {
			std::cout << "Reading in global" << std::endl;
		}
		alnInfoHolderBase<alnInfoGlobal> currentGlobalHolder(dir.first, "GLOBAL");
		if (verbose) {
			std::cout << "Read in  " << currentGlobalHolder.infos_.size() << " GLOBAL alignments" << std::endl;
		}
		localHolder_[currentLocalHolder.gapPars_.getIdentifer()] =
				currentLocalHolder;
		globalHolder_[currentGlobalHolder.gapPars_.getIdentifer()] =
				currentGlobalHolder;
		if (verbose) {
			std::cout << "Read in check: " << globalHolder_[currentGlobalHolder.gapPars_.getIdentifer()].infos_.size() << " GLOBAL alignments for : " << currentGlobalHolder.gapPars_.getIdentifer() << std::endl;
		}
	}
	if(verbose){
		uint32_t alignmentsReadIn = 0;
		for(const auto & alnHold : globalHolder_){
			alignmentsReadIn += alnHold.second.infos_.size();
		}
		for(const auto & alnHold : localHolder_){
			alignmentsReadIn += alnHold.second.infos_.size();
		}
		std::cout << "Global Alignments Holders: " << globalHolder_.size() << std::endl;
		std::cout << "Local Alignments Holders: " << localHolder_.size() << std::endl;
		std::cout << "Read in: " << alignmentsReadIn << " alignments" << std::endl;
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
	if(!njh::files::bfs::exists(masterDirName)){
		njh::files::makeDir(njh::files::MkdirPar(masterDirName));
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
}//njh
