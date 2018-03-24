#pragma once
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
/*
 * alnInfoHolderBase.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */




#include "bibseq/alignment/alnCache/alnInfoGlobal.hpp"
#include "bibseq/alignment/alnCache/alnInfoLocal.hpp"
#include "bibseq/alignment/alignerUtils.h"
#include "bibseq/IO/fileUtils.hpp"

namespace bibseq {
template<typename T>
class alnInfoHolderBase {
public:
	// Constructors
	alnInfoHolderBase() {
	}
	alnInfoHolderBase(const std::string& directoryName,
			const std::string & alnType) :
			alnType_(alnType) {
		auto allFiles = getFiles(directoryName, "", "file", false, false);
		for (const auto &file : allFiles) {
			std::cout << file.first << std::endl;
			std::cout << "INDEX_" + alnType + ".txt" << std::endl;
			if (bib::containsSubString(file.first, "README_" + alnType + ".txt")) {
				readReadMeFile(file.first);
			} else if (bib::containsSubString(file.first,
					"INDEX_" + alnType + ".txt")) {
				readIndex(file.first);
			}
		}
	}
	alnInfoHolderBase(const gapScoringParameters & gapPars,
			const substituteMatrix & scoringArray, const std::string & alnType) :
			gapPars_(gapPars), scoring_(scoringArray), alnType_(alnType) {

	}

	//hasher
  std::hash<std::string> hStr;
  //members
  gapScoringParameters gapPars_;
  substituteMatrix scoring_;
  std::string alnType_;
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, T>>
        infos_;

  /**@brief convert to json representation
   *
   * @return Json::Value object
   */
  Json::Value toJson() const {
  	Json::Value ret;
  	ret["class"] = bib::TypeName::get<alnInfoHolderBase>();
  	ret["gapPars_"] = bib::json::toJson(gapPars_);
  	ret["scoring_"] = bib::json::toJson(scoring_);
  	ret["alnType_"] = bib::json::toJson(alnType_);
  	ret["infos_"] = bib::json::toJson(infos_);
  	return ret;
  }

  // functions
	void addAlnInfo(const std::string& seq1, const std::string& seq2,
									const T& info){
		uint64_t seq1Hash = hStr(seq1);
		uint64_t seq2Hash = hStr(seq2);
		//std::cout << "adding " << seq2 << " to " << seq1 << std::endl;
		addAlnInfo(seq1Hash, seq2Hash, info);
	}

	void addAlnInfo(uint64_t seq1Hash, uint64_t seq2Hash,
									const T& info){
	  if (infos_.find(seq1Hash) == infos_.end()) {
	  	//std::cout << "adding " << seq2Hash << " to " << seq1Hash << std::endl;
	    infos_.emplace(
	        seq1Hash, std::unordered_map<uint64_t, T>{{seq2Hash, info}});
	  } else if (infos_[seq1Hash].find(seq2Hash) != infos_[seq1Hash].end()) {
	  	//@TODO Doing nothing now if there are duplicates but perhaps should do a check to see if they are the same at least
	    //std::cout << "already added " << std::endl << seq2Hash << std::endl;
	    //std::cout << "to " << std::endl << seq1Hash << std::endl;
	  } else {
	    infos_[seq1Hash].emplace(seq2Hash, info);
	  }
	}


	bool checkForAlnInfo(const std::string& seq1, const std::string& seq2){
		//std::cout << "bool checkForAlnInfo(const std::string& seq1, const std::string& seq2)" << std::endl;
	  uint64_t seq1Hash = hStr(seq1);
	  uint64_t seq2Hash = hStr(seq2);
		return checkForAlnInfo(seq1Hash, seq2Hash);
	}
	bool checkForAlnInfo(uint64_t seq1Hash, uint64_t seq2Hash){
		//std::cout << "bool checkForAlnInfo(uint64_t seq1Hash, uint64_t seq2Hash)" << std::endl;
		auto check = infos_.find(seq1Hash);
	  if (check == infos_.end()) {
	    return false;
	  } else if (check->second.find(seq2Hash) == check->second.end()) {
	    return false;
	  } else {
	    return true;
	  }
	}

	bool getAlnInfo(const std::string& seq1, const std::string& seq2, T & info){
		return getAlnInfo(hStr(seq1), hStr(seq2), info);
	}

	bool getAlnInfo(uint64_t seq1Hash, uint64_t seq2Hash, T & info){
		auto check = infos_.find(seq1Hash);
	  if (check == infos_.end()) {
	    return false;
	  } else {
	  	auto check2 = check->second.find(seq2Hash);
	  	if (check2 == check->second.end()) {
	  		return false;
	  	}else{
	  		info = check2->second;
	  		return true;
	  	}
	  }
	}

	void readReadMeFile(const std::string& filename){
	  //std::cout << "reading "+ alnType_ +"readme" << std::endl;
	  std::ifstream inFile;
	  inFile.open(filename);
	  std::vector<std::vector<int32_t>> arrayInVectors;
	  bool readingGaps = false;
	  bool readingArray = false;
	  for (std::string line; std::getline(inFile, line);) {
	    if (readingGaps) {
	      //std::cout << "gapLine: " << line << std::endl;
	      auto gapParsInfo =
	          convertStringToVector<int32_t>(line, ",");
	      //printVector(gapParsInfo);
	      gapPars_ =
	          gapScoringParameters(gapParsInfo[0], gapParsInfo[1],
	          		                    gapParsInfo[2],gapParsInfo[3],
																 gapParsInfo[4], gapParsInfo[5],
																 gapParsInfo[6], gapParsInfo[7],
																 gapParsInfo[8], gapParsInfo[9]);
	      readingGaps = false;
	    }
	    if (readingArray) {
	      arrayInVectors.emplace_back(convertStringToVector<int32_t>(line, ","));
	    }
	    if (line == "GapPars:") {
	      readingArray = false;
	      readingGaps = true;
	    } else if (line == "ScoringArray:") {
	      readingArray = true;
	      readingGaps = false;
	    }
	  }
		std::vector<uint32_t> positions(arrayInVectors.size(), 0);
		std::iota(positions.begin(), positions.end(), 0);
	  for (auto i : positions) {
	    for (auto j : positions) {
	      scoring_.mat_[i][j] = arrayInVectors[i][j];
	    }
	  }
	}
	void readIndex(const std::string& indexFilename){
	  //std::cout << "Reading " + alnType_  + " index file" << std::endl;
	  std::ifstream inFile;
	  inFile.open(indexFilename);
	  std::vector<std::string> info(800, "");
	  std::vector<std::string> inputGapInfo(3, "");
	  std::string out;
	  std::stringstream ss;
	  std::vector<gapInfo> gInfos;
	  for (std::string line; std::getline(inFile, line);) {
	    inputGapInfo.clear();
	    gInfos.clear();
	    ss.clear();
	    ss << line;
	  	uint32_t index = 0;
	    while (index < 2) {
	      ss >> out;
	      info[index] = out;
	      ++index;
	    }
	    addAlnInfo(std::stoul(info[0]), std::stoul(info[1]),
	    		T::readInfo(ss, info, inputGapInfo,gInfos, out ) );
	  }
	}
	void writeOutInfos(std::string parentDirectory,
										 std::string outDirectoryName) const{
		bfs::path outDir = bib::files::join(parentDirectory,
				bib::appendAsNeededRet(outDirectoryName, "/"));

		if(!bib::files::bfs::exists(outDir)){
			bib::files::makeDir(bib::files::MkdirPar(outDir));
		}
		// write out info for readMeFile
		bfs::path readMeFilename = bib::files::join(outDir.string(), "README_" + alnType_ + ".txt");
		if (!bib::files::bfs::exists(readMeFilename)) {
			std::ofstream outReadMeFile;
			openTextFile(outReadMeFile, readMeFilename, ".txt", false,
									 false);
			outReadMeFile << "GapPars:" << std::endl;
			gapPars_.writePars(outReadMeFile);
			outReadMeFile << "ScoringArray:" << std::endl;
			std::vector<uint32_t> positions(scoring_.mat_.size(), 0);
			std::iota(positions.begin(), positions.end(), 0);
			for (auto i : positions) {
				for (auto j : positions) {
					if (0 != j) {
						outReadMeFile << ",";
					}
					outReadMeFile << scoring_.mat_[i][j];
				}
				outReadMeFile << std::endl;
			}
		}
		std::ofstream indexFile;
		bfs::path indexFilename = bib::files::join(outDir.string(), "INDEX_" + alnType_ + ".txt");
		indexFile.open(indexFilename.string(), std::ios::app);

		// openTextFile(indexFile, , ".txt" ,false, true);
		uint64_t alnCount = 0;
		uint64_t alnAddFromFileCount = 0;
		uint64_t alnFromThisCount = 0;
		for (const auto &infoMap : infos_) {
			for (const auto &infoMapSub : infoMap.second) {
				++alnCount;
				if (!infoMapSub.second.addFromFile_) {
					infoMapSub.second.writeInfoSingleLine(indexFile, infoMap.first,
																								infoMapSub.first);
					++alnFromThisCount;
				} else {
					++alnAddFromFileCount;
				}
			}
		}
	}

	void mergeOtherHolder(const alnInfoHolderBase<T> & other ){
		bool fail = false;
		std::stringstream ss;
		if(other.alnType_ != alnType_){
			ss << "Error in mergeOtherHolder, trying to merge aln holder of a different type" << std::endl;
			ss << "current holder type: " << alnType_ << ", trying to add type: " << other.alnType_ << std::endl;
			fail = true;
		}
		if(other.gapPars_ != gapPars_){
			ss << "Error in mergeOtherHolder, trying to merge aln holder of a gap parameters" << std::endl;
			ss << "current gap type: " << gapPars_.getIdentifer() << ", trying to add type: " << other.gapPars_.getIdentifer() << std::endl;
			fail = true;
		}
		if(other.scoring_.mat_ != scoring_.mat_){
			ss << "Error in mergeOtherHolder, trying to merge aln holder with different scoring" << std::endl;
			ss << "current scoring: " << std::endl;
			scoring_.printScores(ss);
			ss << "other scoring: " << std::endl;
			other.scoring_.printScores(ss);
			fail = true;
		}
		if(fail){
			throw std::runtime_error{ss.str()};
		}
		for(const auto & s1 : other.infos_){
			for(const auto & s2 : s1.second){
				if(s2.second.addFromFile_){
					continue;
				}
				if(!checkForAlnInfo(s1.first,s2.first)){
					//does not already contain this alignment
					addAlnInfo(s1.first, s2.first, s2.second);
				}
			}
		}
	}
};


}  // namespace bibseq


