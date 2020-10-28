#pragma once
//
//  seqToolsUtils.hpp
//
//  Created by Nicholas Hathaway on 4/27/13.
//
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
#include "njhseq/objects/seqObjects/readObject.hpp"
#include "njhseq/objects/seqObjects/sffObject.hpp"

#include "njhseq/seqToolsUtils/aminoAcidInfo.hpp"
#include "njhseq/alignment.h"
#include "njhseq/helpers/seqUtil.hpp"
#include "njhseq/utils.h"
#include "njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp"
#include "njhseq/objects/helperObjects/probabilityProfile.hpp"
#include "njhseq/objects/dataContainers/graphs/ReadCompGraph.hpp"


#include "njhseq/objects/Meta/MetaDataInName.hpp"

#include <njhcpp/graphics/colorUtils.hpp>


namespace njhseq {

table getSeqPortionCounts(const SeqIOOptions & opts, size_t position, uint32_t size, bool back = false);

void processRunCutoff(uint32_t& runCutOff, const std::string& runCutOffString,
		int counter);

uint32_t processRunCutoff(const std::string& runCutOffString, uint64_t counter);

template<typename T>
Json::Value genDetailMinTreeData(const std::vector<T> & reads,
		uint32_t numThreads) {
	uint64_t maxSize = 0;
	readVec::getMaxLength(reads, maxSize);
	aligner alignerObj(maxSize, gapScoringParameters(5, 1),
			substituteMatrix(2, -2));
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	std::mutex alignerLock;
	return genDetailMinTreeData(reads, alignerObj, aligners, alignerLock,
			numThreads, comparison(), false, false, false);
}

template<typename T>
Json::Value genDetailMinTreeData(const std::vector<T> & reads,
		uint32_t numThreads,
		const comparison &allowableErrors, bool settingEventsLimits) {
	uint64_t maxSize = 0;
	readVec::getMaxLength(reads, maxSize);
	aligner alignerObj(maxSize, gapScoringParameters(5, 1),
			substituteMatrix(2, -2));
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	std::mutex alignerLock;
	return genDetailMinTreeData(reads, alignerObj, aligners, alignerLock,
			numThreads, allowableErrors, settingEventsLimits, false, false);
}

template<typename T>
Json::Value genDetailMinTreeData(const std::vector<T> & reads,
		aligner & alignerObj,
		std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
		std::mutex & alignerLock, uint32_t numThreads,
		const comparison &allowableErrors, bool settingEventsLimits,
		bool justBest, bool doTies
		){
	auto graph = genReadComparisonGraph(reads, alignerObj, aligners, alignerLock,
			numThreads);
	std::vector<std::string> popNames;
	for (const auto & n : graph.nodes_) {
		if (n->on_) {
			popNames.emplace_back(n->name_);
		}
	}
	auto nameColors = getColorsForNames(popNames);

	if(justBest){
		if (settingEventsLimits) {
			graph.turnOffEdgesWithComp(allowableErrors,
					[](const comparison & comp1, const comparison & cutOff) {
						//std::cout << comp1.toJson() << std::endl;
						return comp1.distances_.getNumOfEvents(true) >= cutOff.distances_.overLappingEvents_;
					});
		}
		graph.setJustBestConnection(doTies);
	}else{
		if(settingEventsLimits){
			graph.turnOffEdgesWithComp(allowableErrors,
					[](const comparison & comp1, const comparison & cutOff){
				//std::cout << comp1.toJson() << std::endl;
				return comp1.distances_.getNumOfEvents(true) >= cutOff.distances_.overLappingEvents_;
			});
		}else{
			comparison maxEvents = graph.setMinimumEventConnections();
		}
	}
	auto treeData = graph.toD3Json(njh::color("#000000"), nameColors);
	return treeData;
}

template<typename T>
uint32_t getMismatches(const T & read1,
				const T & read2,
				aligner alignerObj){
	alignerObj.alignCache(read1.seqBase_,read2.seqBase_, false);
	alignerObj.profilePrimerAlignment(read1.seqBase_, read2.seqBase_);
	return alignerObj.comp_.hqMismatches_;
};
Json::Value genMinTreeData(const std::vector<readObject> & reads,
		aligner & alignerObj,
		std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
		std::mutex & alignerLock, uint32_t numThreads);

Json::Value genMinTreeData(const std::vector<readObject> & reads);
Json::Value genMinTreeData(const std::vector<readObject> & reads, aligner & alignerObj);
std::string genHtmlStrForPsuedoMintree(std::string jsonFileName);
std::string genHtmlStrForPsuedoMintree(std::string jsonFileName, std::string javaScriptFile);





uint32_t countSeqs(const SeqIOOptions & opts, bool verbose);

void setUpSampleDirs(
    const std::string& sampleNamesFilename,
		const std::string& mainDirectoryName,
		bool separatedDirs,
		bool verbose = true);

std::vector<char> getAAs(bool noStopCodon);

template <typename T>
static std::vector<readObject> createCondensedObjects(std::vector<T> reads) {
  std::vector<readObject> ans;
  readVec::allSetCondensedSeq(reads);
  for (const auto& read : reads) {
    ans.emplace_back(readObject(seqInfo(read.seqBase_.name_, read.condensedSeq,
                                        read.condensedSeqQual)));
  }
  return ans;
}


template <typename T>
std::unordered_map<std::string, std::string> renameReadNames(std::vector<T>& reads, const std::string& stub,
                     bool processed, bool keepChimeraFlag = true,
                     bool keepCompFlag = true,
                     const std::string& sortBy = "none") {
	uint64_t count = 0;
  if ("none" != sortBy) {
    readVecSorter::sortReadVector(reads, sortBy);
  }
  VecStr originalNames = readVec::getNames(reads);
  uint64_t maxSize = reads.size();
  for (auto& seq : reads) {
    bool chimera = getSeqBase(seq).name_.find("CHI") != std::string::npos;
    bool comp = getSeqBase(seq).name_.find("_Comp") != std::string::npos;
    getSeqBase(seq).name_ = stub;
    if (chimera && keepChimeraFlag) {
    	getSeqBase(seq).markAsChimeric();
    }
    if (comp && keepCompFlag) {
    	getSeqBase(seq).name_ += "_Comp";
    }
    getSeqBase(seq).name_ += "." + leftPadNumStr(count, maxSize);
    if (processed) {
    	getSeqBase(seq).name_ += "_t" + estd::to_string(getSeqBase(seq).cnt_);
    }
    ++count;
  }
  std::unordered_map<std::string, std::string> ret;
  for(const auto pos : iter::range(reads.size())){
  	ret[getSeqBase(reads[pos]).name_] = originalNames[pos];
  }
  return ret;
}

template <typename T>
void renameReadNamesNewClusters(std::vector<T>& reads, const std::string& stub,
                                bool processed, bool keepChimeraFlag = true,
                                bool keepCompFlag = true,
                                const std::string& sortBy = "none") {
  renameReadNames(reads, stub, processed, keepChimeraFlag, keepCompFlag, sortBy);
  for (auto& read : reads) {
    read.reads_.front()->seqBase_.name_ = read.seqBase_.name_;
  }
}






template <typename T>
VecStr findLongestSharedSeqFromReads(const std::vector<T>& reads) {
  VecStr seqs;
  for (const auto& rIter : reads) {
    seqs.push_back(rIter.seqBase_.seq_);
  }
  return seqUtil::findLongestShardedMotif(seqs);
}


void processKrecName(readObject& read, bool post);



/////////tools for finding additional location output
std::string makeIDNameComparable(const std::string& idName);

std::string findIdNameFromFileName(const std::string& filename);

std::string processFileNameForID(const std::string& fileName);

std::string findAdditonalOutLocation(const std::string& locationFile,
                                     const std::string& fileName);


///


template <typename READ, typename REF>
std::vector<baseReadObject> alignToSeq(const std::vector<READ>& reads,
                                       const REF& reference, aligner& alignObj,
                                       bool local, bool usingQuality) {
  std::vector<baseReadObject> output;
  output.emplace_back(baseReadObject(reference));
  for (const auto& read : reads) {
    alignObj.alignCache(reference, read, local);
    output.emplace_back(baseReadObject(seqInfo(
        read.seqBase_.name_ + "_score:" + std::to_string(alignObj.parts_.score_),
        alignObj.alignObjectB_.seqBase_.seq_,
        alignObj.alignObjectB_.seqBase_.qual_)));
  }
  return output;
}
template <typename READ, typename REF>
std::vector<baseReadObject> alignToSeqVec(const std::vector<READ>& reads,
                                          const REF& reference,
                                          aligner& alignObj, bool local,
                                          bool usingQuality) {
  std::vector<baseReadObject> output;
  output.emplace_back(baseReadObject(reference));
  for (const auto& read : reads) {
    alignObj.alignCache(reference, read, local);
    output.emplace_back(baseReadObject(seqInfo(
        read.seqBase_.name_ + "_score:" + std::to_string(alignObj.parts_.score_),
        alignObj.alignObjectB_.seqBase_.seq_,
        alignObj.alignObjectB_.seqBase_.qual_)));
  }
  return output;
}
template <typename READ, typename REF>
VecStr alignToSeqStrings(const std::vector<READ>& reads, const REF& reference,
                         aligner& alignObj, bool local, bool usingQuality) {
  VecStr output;
  output.push_back(reference.seqBase_.seq_);
  for (const auto read : reads) {
    alignObj.alignCache(reference, read, local);
    output.push_back(alignObj.alignObjectB_.seqBase_.seq_);
  }
  return output;
}






template <typename T>
uint32_t seqQualSizeAgreementCheck(const std::vector<T>& reads) {
  uint32_t count = 0;
  for (const auto& read : reads) {
    if (read.seqBase_.seq_.size() != read.seqBase_.qual_.size()) {
      std::cout << "name:" << read.seqBase_.name_ << std::endl;
      std::cout << "seq:" << read.seqBase_.seq_ << std::endl;
      std::cout << "qual:" << read.seqBase_.qual_ << std::endl;
      std::cout << "sSize:" << read.seqBase_.seq_.size() << std::endl;
      std::cout << "qSize:" << read.seqBase_.qual_.size() << std::endl;
      ++count;
    }
  }
  std::cout << getPercentageString(count, reads.size()) << std::endl;
  return count;
}
template <typename T>
std::unordered_map<std::string, njh::color> getColorsForNames(
    const std::vector<T>& reads, double sat, double lum) {
  std::unordered_map<std::string, njh::color> colorsForName;
  uint32_t count = 0;
  std::vector<njh::color> colors = njh::evenHuesAll(sat, lum, len(reads));
  for (const auto& read : reads) {
    colorsForName[read.getReadId()] = colors[count];
    ++count;
  }
  return colorsForName;
}



template<typename READ>
std::vector<READ> vecStrToReadObjs(const VecStr & strs, const std::string & stubName){
	std::vector<READ> ans;
	for(const auto strPos : iter::range(strs.size())){
		ans.emplace_back(READ(seqInfo(stubName + "." + leftPadNumStr(strPos, strs.size()), strs[strPos])));
	}
	return ans;
}

template<typename T>
VecStr readObjsToVecStr(const std::vector<T> & vec){
	VecStr ans;
	for(const auto & read : vec){
		ans.emplace_back(read.seqBase_.seq_);
	}
	return ans;
}

VecStr createDegenStrs(const std::string & str);
static const std::unordered_map<char, std::vector<char>> degenBaseExapnd = {
		{ 'N', std::vector<char>{'A', 'C','G', 'T'}},
		{ 'K', std::vector<char>{'G', 'T'}},
		{ 'Y', std::vector<char>{'C', 'T'}},
		{ 'W', std::vector<char>{'A', 'T'}},
		{ 'S', std::vector<char>{'C', 'G'}},
		{ 'R', std::vector<char>{'A', 'G'}},
		{ 'M', std::vector<char>{'A', 'C'}},
		{ 'B', std::vector<char>{'C', 'G', 'T'}},
		{ 'D', std::vector<char>{'A', 'G', 'T'}},
		{ 'H', std::vector<char>{'A', 'C', 'T'}},
		{ 'V', std::vector<char>{'A', 'C', 'G'}}
};

template<typename READ, typename REF>
bool checkPossibleChiByRefs(const READ & read, const std::vector<REF> & refSeqs, table& outInfo, aligner & alignerObj,
		const comparison & chiOverlap, bool breakAtFirst,   bool weightHomopolymers ){
	bool foundAnExactMatch = false;
	bool foundAChimera = false;
	std::string firstChiName = "";
	std::string secondChiName = "";
	uint32_t inflectionPoint = UINT32_MAX;
	uint32_t inflectionPointPar1 = UINT32_MAX;
	uint32_t inflectionPointPar2 = UINT32_MAX;
	for(const auto  refPos : iter::range(len(refSeqs))){
		auto & ref = refSeqs[refPos];
		alignerObj.alignCache(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_);
		if(alignerObj.comp_.distances_.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
			return foundAChimera;
		}
	}
	for(const auto refPos : iter::range(len(refSeqs))){
		auto & ref = refSeqs[refPos];
		alignerObj.alignCache(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_);
		if(alignerObj.comp_.distances_.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
		} else if(alignerObj.comp_.distances_.mismatches_.size() > 0 && !foundAnExactMatch){
			auto savedMismatches = alignerObj.comp_.distances_.mismatches_;
			//check front
			alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, false, true,
					false, 0,
					alignerObj.getAlignPosForSeqBPos(
							savedMismatches.begin()->second.seqBasePos));
			bool passFront = chiOverlap.passErrorProfile(alignerObj.comp_);

			//check back
			alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, false, true,
					false,
					alignerObj.getAlignPosForSeqBPos(
							savedMismatches.rbegin()->second.seqBasePos + 1));
			bool passBack = chiOverlap.passErrorProfile(alignerObj.comp_);

			auto firstRefAlignA = alignerObj.alignObjectA_;
			auto firstRefAlignB = alignerObj.alignObjectB_;
			//std::cout << "pass front: " << passFront << std::endl;
			//std::cout << "pass back: " << passBack << std::endl;
			if(passFront){
				for(const auto  secondRefPos : iter::range(refPos + 1, len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}

					alignerObj.alignCache(secondRef.seqBase_, read.seqBase_, false);
					//check to see if from the mismatch on from the other one ref matches
					//to another ref
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false);
					if(alignerObj.comp_.distances_.mismatches_.size() < 1){
						continue;
					}
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false,
												alignerObj.getAlignPosForSeqBPos(savedMismatches.begin()->second.seqBasePos));

					//alignerObj.errors_.printDescription(std::cout, true);
					if(chiOverlap.passErrorProfile(alignerObj.comp_)){
						firstChiName = ref.seqBase_.name_;
						secondChiName = secondRef.seqBase_.name_;
						inflectionPoint = savedMismatches.begin()->second.seqBasePos;
						inflectionPointPar1 = savedMismatches.begin()->second.refBasePos;
						std::string refTwoPortion = alignerObj.alignObjectA_.seqBase_.seq_.substr(0, alignerObj.getAlignPosForSeqAPos(savedMismatches.begin()->second.refBasePos) + 1);
						inflectionPointPar2 = removeCharReturn(refTwoPortion, '-').size() - 1;

						foundAChimera = true;
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, estd::to_string(read.seqBase_.frac_),
							estd::to_string(inflectionPoint),
							firstChiName, "1", estd::to_string(inflectionPointPar1),
							secondChiName, "1", estd::to_string(inflectionPointPar2)});
						if(foundAChimera && breakAtFirst){
							//VecStr{"readName", "fraction", "par1", "par1Frac", "par2", "par2Frac", "inflectionPoint"}
							break;
						}
					}
				}
			}
			if(foundAChimera && breakAtFirst){
				break;
			}
			if(passBack){
				for(const auto secondRefPos : iter::range(refPos + 1,len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}
					alignerObj.alignCache(secondRef.seqBase_, read.seqBase_, false);
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false);
					if(alignerObj.comp_.distances_.mismatches_.size() < 1){
						continue;
					}
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, false, true, false,
												0, alignerObj.getAlignPosForSeqBPos(savedMismatches.rbegin()->second.seqBasePos + 1));
					if(chiOverlap.passErrorProfile(alignerObj.comp_)){
						foundAChimera = true;
						firstChiName = ref.seqBase_.name_;
						secondChiName = secondRef.seqBase_.name_;
						inflectionPoint = savedMismatches.rbegin()->second.seqBasePos;
						inflectionPointPar1 = savedMismatches.rbegin()->second.refBasePos;
						std::string refTwoPortion = alignerObj.alignObjectA_.seqBase_.seq_.substr(0, alignerObj.getAlignPosForSeqAPos(savedMismatches.rbegin()->second.refBasePos) + 1);
						inflectionPointPar2 = removeCharReturn(refTwoPortion, '-').size() - 1;
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, estd::to_string(read.seqBase_.frac_),
							estd::to_string(inflectionPoint),
							firstChiName, "1", estd::to_string(inflectionPointPar1),
							secondChiName, "1", estd::to_string(inflectionPointPar2)});
						if(foundAChimera && breakAtFirst){
							break;
						}
					}
				}
			}
		}
	}
	return foundAChimera;
}



table getSeqPosTab(const std::string & str);









struct ExtractionInfo {
	ExtractionInfo(const table & allStatsTab, const table & allprofileTab,
			const table & profileBySampTab) :
			allStatsTab_(allStatsTab), allProfileTab_(allprofileTab), profileBySampTab_(
					profileBySampTab) {

	}
	ExtractionInfo() {
	}

	table allStatsTab_;
	table allProfileTab_;
	table profileBySampTab_;
};

ExtractionInfo collectExtractionInfo(const std::string & dirName,
		const std::string & indexToDir, const std::string & sampNames);

ExtractionInfo collectExtractionInfoDirectName(const std::string & dirName,
		const std::string & indexToDir, const std::string & sampNames);




class readsWithSnps {
public:

	readsWithSnps(const std::map<uint32_t, mismatch> & mismatches,
			uint64_t firstReadPos) :
			mismatches_(mismatches), readPositions_( { firstReadPos }) {
		mismatchId_ = "";
		for(const auto & hm : mismatches_){
			mismatchId_.append(estd::to_string(hm.first) + estd::to_string(hm.second.seqBase));
		}
	}
	std::string mismatchId_;
	std::map<uint32_t, mismatch> mismatches_;
	std::vector<uint64_t> readPositions_;

	void addReadPos(uint64_t readPos){
		readPositions_.emplace_back(readPos);
	}

	using size_type = std::vector<uint64_t>::size_type;

};

template<>
inline readsWithSnps::size_type len(const readsWithSnps & reads){
	return reads.readPositions_.size();
}


template<typename T>
std::unordered_map<std::string, std::vector<T>> splitSeqsByMetaField(const std::vector<T> & seqs, const std::string & field){
	std::unordered_map<std::string, std::vector<T>> ret;
	for(const auto & seq : seqs){
		MetaDataInName meta(getSeqBase(seq).name_);
		ret[meta.getMeta(field)].emplace_back(seq);
	}
	return ret;
}

template<typename T>
std::vector<char> determineAlphabet(const std::vector<T> & seqs){
	std::set<char> retSet;
	for(const auto & seq : seqs){
		for(const char c : getSeqBase(seq).seq_){
			retSet.emplace(c);
		}
	}
	return std::vector<char>{retSet.begin(), retSet.end()};
}

}  // namespace njh


