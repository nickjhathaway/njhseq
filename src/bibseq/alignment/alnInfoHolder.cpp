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

//  alnInfoHolder.cpp
//
//  Created by Nicholas Hathaway on 1/13/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "alnInfoHolder.hpp"

namespace bibseq {





void alnInfoLocal::writeInfoSingleLine(std::ostream &indexFile, uint64_t seq1,
                                       uint64_t seq2) const {
  indexFile << seq1 << " " << seq2 << " " << score_ << " " << localAStart_ << " " << localASize_
            << " " << localBStart_ << " " << localBSize_;
  for (const auto &g : gapInfos_) {
    indexFile << " ";
    g.writeInfoNoEnd(indexFile);
  }
  indexFile << std::endl;
}

alnInfoLocal alnInfoLocal::readInfo(std::stringstream& ss,
		std::vector<std::string>& info, std::vector<std::string>& inputGapInfo,
		std::vector<gapInfo> & gInfos, std::string & out){

	uint32_t index = 2;
  while (!ss.eof()) {
    ss >> out;
    if(index >= info.size()){
    	addOtherVec(info, VecStr(index));
    }
    info[index] = out;
    ++index;
  }
  if (index > 7) {
  	std::vector<uint32_t> positions(index - 7);
  	std::iota(positions.begin(), positions.end(), 7);
    for (auto i : positions) {
      inputGapInfo = tokenizeString(info[i], ",");
      gInfos.emplace_back(gapInfo(std::stoi(inputGapInfo[0]),
                                  std::stoi(inputGapInfo[1]),
                                  std::stoi(inputGapInfo[2])));
    }
  }
  return alnInfoLocal(gInfos, std::stoi(info[3]), std::stoi(info[4]),
                            std::stoi(info[5]), std::stod(info[6]),
                            std::stoi(info[2]), true);
}

void alnInfoGlobal::writeInfoSingleLine(std::ostream &indexFile, uint64_t seq1,
                                        uint64_t seq2) const {
  indexFile << seq1 << " " << seq2 << " " << score_;
  for (const auto &g : gapInfos_) {
    indexFile << " ";
    g.writeInfoNoEnd(indexFile);
  }
  indexFile << std::endl;
}

alnInfoGlobal alnInfoGlobal::readInfo(std::stringstream& ss,
   		std::vector<std::string>& info, std::vector<std::string>& inputGapInfo,
   		std::vector<gapInfo> & gInfos, std::string & out){

	uint32_t index = 2;
  while (!ss.eof()) {
    ss >> out;
    if(index >= info.size()){
    	addOtherVec(info, VecStr(index));
    }
    info[index] = out;
    ++index;
  }

  if (index > 3) {
  	std::vector<uint32_t> positions(index - 3);
  	std::iota(positions.begin(), positions.end(), 3);
    for (auto i : positions) {
      inputGapInfo = tokenizeString(info[i], ",");
      gInfos.emplace_back(gapInfo(std::stoi(inputGapInfo[0]),
                                  std::stoi(inputGapInfo[1]),
                                  std::stoi(inputGapInfo[2])));
    }
  }
  return alnInfoGlobal(gInfos, std::stod(info[2]), true);
}



alnInfoMasterHolder::alnInfoMasterHolder() {}

alnInfoMasterHolder::alnInfoMasterHolder(const gapScoringParameters & gapPars,
		const substituteMatrix & scoringArray) :
		localHolder_(
				std::unordered_map<std::string, alnInfoHolderBase<alnInfoLocal>> { {
						gapPars.getIdentifer(), alnInfoHolderBase<alnInfoLocal>(gapPars,
								scoringArray, "LOCAL") } }), globalHolder_(
				std::unordered_map<std::string, alnInfoHolderBase<alnInfoGlobal>> { {
						gapPars.getIdentifer(), alnInfoHolderBase<alnInfoGlobal>(gapPars,
								scoringArray, "GLOBAL") } }) {
}
;


alnInfoMasterHolder::alnInfoMasterHolder(const std::string &masterDirName, const gapScoringParameters & gapPars,
	  const substituteMatrix & scoringArray, bool verbose ) {
	alignment::alnCacheDirSearchLock.lock();
	auto fullPath = bib::files::normalize(masterDirName).string();
	auto search = alignment::alnCacheDirLocks.find(fullPath);
	if(search == alignment::alnCacheDirLocks.end()){
		alignment::alnCacheDirLocks.emplace(std::piecewise_construct, std::make_tuple(fullPath), std::tuple<>());
	}
	auto realSearch = alignment::alnCacheDirLocks.find(fullPath);
	std::lock_guard<std::mutex> lock(realSearch->second);
	alignment::alnCacheDirSearchLock.unlock();
  if(verbose) {
  	std::cout << "Reading in Previous Alignments" << std::endl;
  }
  auto allDirectories = getFiles(masterDirName, "", "directory", false, false);
  //std::cout << bib::bashCT::boldGreen("here?") << std::endl;
  for (const auto &dir : allDirectories) {
  	//std::cout << bib::bashCT::boldGreen("here1?") << std::endl;
  	if(verbose){
  		std::cout << "Reading from " << dir.first << std::endl;
  		std::cout << "Reading in local" << std::endl;
  	}
  	alnInfoHolderBase<alnInfoLocal> currentLocalHolder(dir.first, "LOCAL");
    if(verbose){
    	std::cout << "Reading in global" << std::endl;
    }
    alnInfoHolderBase<alnInfoGlobal> currentGlobalHolder(dir.first, "GLOBAL");

    //std::cout << "local: " << currentLocalHolder.gapPars_.getIdentifer()
    //          << std::endl;
    //std::cout << "global: " << currentGlobalHolder.gapPars_.getIdentifer()
    //          << std::endl;

    localHolder_[currentLocalHolder.gapPars_.getIdentifer()] =
        currentLocalHolder;
    globalHolder_[currentGlobalHolder.gapPars_.getIdentifer()] =
        currentGlobalHolder;
  }
  //std::cout << bib::bashCT::boldGreen("here2?") << std::endl;
  if(allDirectories.empty()){
  	//std::cout << bib::bashCT::boldGreen("here3?") << std::endl;
  	//std::cout << gapPars.getIdentifer() << std::endl;
  	//std::cout << bib::bashCT::boldGreen("here3.1?") << std::endl;
  	auto lHolder = alnInfoHolderBase<alnInfoLocal>(gapPars,scoringArray, "LOCAL");
  	//std::cout << bib::bashCT::boldGreen("here3.2?") << std::endl;
  	auto gHolder = alnInfoHolderBase<alnInfoGlobal>(gapPars,scoringArray, "GLOBAL");
  	//std::cout << bib::bashCT::boldGreen("here3.3?") << std::endl;
  	localHolder_[gapPars.getIdentifer()] = lHolder;
  	//localHolder_[gapPars.getIdentifer()] = alnInfoHolderBase<alnInfoLocal>(gapPars,scoringArray, "LOCAL");
  	//std::cout << bib::bashCT::boldGreen("here3.4?") << std::endl;
  	globalHolder_[gapPars.getIdentifer()] = gHolder;
  	//globalHolder_[gapPars.getIdentifer()] = alnInfoHolderBase<alnInfoGlobal>(gapPars,scoringArray, "GLOBAL");
  	//std::cout << bib::bashCT::boldGreen("here3.5?") << std::endl;
  	/*localHolder_ =
  					std::unordered_map<std::string, alnInfoHolderBase<alnInfoLocal>> { {
  							gapPars.getIdentifer(), alnInfoHolderBase<alnInfoLocal>(gapPars,scoringArray, "LOCAL") } };
  	globalHolder_ =
  					std::unordered_map<std::string, alnInfoHolderBase<alnInfoGlobal>> { {
  							gapPars.getIdentifer(), alnInfoHolderBase<alnInfoGlobal>(gapPars,scoringArray, "GLOBAL") } };*/
  	//std::cout << bib::bashCT::boldGreen("here4?") << std::endl;
  }
  //std::cout << bib::bashCT::boldGreen("here5?") << std::endl;
}


void alnInfoMasterHolder::write(const std::string &masterDirName, bool verbose ) {

	auto fullPath = bib::files::normalize(masterDirName).string();
	alignment::alnCacheDirSearchLock.lock();
	auto search = alignment::alnCacheDirLocks.find(fullPath);
	if(search == alignment::alnCacheDirLocks.end()){
		alignment::alnCacheDirLocks.emplace(std::piecewise_construct, std::make_tuple(fullPath), std::tuple<>());
	}
	auto realSearch = alignment::alnCacheDirLocks.find(fullPath);

	std::lock_guard<std::mutex> lock(realSearch->second);
	alignment::alnCacheDirSearchLock.unlock();

  int directoryStatus =
      mkdir(masterDirName.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  if (directoryStatus != 0) {
    // if directory already exits do nothing i guess
  } else {
    chmod(masterDirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  // write out the alignment infos
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
