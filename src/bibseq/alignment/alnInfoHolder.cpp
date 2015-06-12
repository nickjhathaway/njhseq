//
//  alnInfoHolder.hpp
//  sequenceTools
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
	  const substituteMatrix & scoringArray):localHolder_(std::unordered_map<std::string, alnInfoHolderBase<alnInfoLocal>> {{gapPars.getIdentifer(),alnInfoHolderBase<alnInfoLocal>(gapPars,scoringArray, "LOCAL")} }),
	  		globalHolder_(std::unordered_map<std::string, alnInfoHolderBase<alnInfoGlobal>>{{gapPars.getIdentifer(),alnInfoHolderBase<alnInfoGlobal>(gapPars,scoringArray, "GLOBAL")} }){};


alnInfoMasterHolder::alnInfoMasterHolder(const std::string &masterDirName, const gapScoringParameters & gapPars,
	  const substituteMatrix & scoringArray, bool verbose ) {
	auto fullPath = bib::files::normalize(masterDirName).string();
	auto search = alignment::alnCacheDirLocks.find(fullPath);
	if(search == alignment::alnCacheDirLocks.end()){
		alignment::alnCacheDirLocks.emplace(std::piecewise_construct, std::make_tuple(fullPath), std::tuple<>());
	}
	auto realSearch = alignment::alnCacheDirLocks.find(fullPath);
	std::lock_guard<std::mutex> lock(realSearch->second);

  if(verbose) std::cout << "Reading in previous alignments" << std::endl;
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

    //std::cout << "local: " << currentLocalHolder.gapPars_.getIdentifer()
    //          << std::endl;
    //std::cout << "global: " << currentGlobalHolder.gapPars_.getIdentifer()
    //          << std::endl;

    localHolder_[currentLocalHolder.gapPars_.getIdentifer()] =
        currentLocalHolder;
    globalHolder_[currentGlobalHolder.gapPars_.getIdentifer()] =
        currentGlobalHolder;
  }
  if(allDirectories.empty()){
  	*this = alnInfoMasterHolder(gapPars, scoringArray);
  }
}


void alnInfoMasterHolder::write(const std::string &masterDirName, bool verbose ) {
	auto fullPath = bib::files::normalize(masterDirName).string();
	auto search = alignment::alnCacheDirLocks.find(fullPath);
	if(search == alignment::alnCacheDirLocks.end()){
		alignment::alnCacheDirLocks.emplace(std::piecewise_construct, std::make_tuple(fullPath), std::tuple<>());
	}
	auto realSearch = alignment::alnCacheDirLocks.find(fullPath);
	std::lock_guard<std::mutex> lock(realSearch->second);


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
