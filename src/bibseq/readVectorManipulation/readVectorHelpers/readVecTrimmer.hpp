#pragma once
//
//  trimmer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/22/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
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

#include "bibseq/alignment.h"
#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"
#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/collapseObjects/opts/IterPar.hpp"
#include "bibseq/objects/kmer/kmerInfo.hpp"
#include "bibseq/objects/helperObjects/probabilityProfile.hpp"

namespace bibseq {

struct FullTrimReadsPars {
	struct trimSeqPars {
		bool includeSequence_;
		bool sequenceToLowerCase_;
		bool removePreviousSameBases_;
		bool local_ = true;
		uint32_t within_ = std::numeric_limits<uint32_t>::max();
	};
  // parameters
	FullTrimReadsPars();
	void initForKSharedTrim();

  uint32_t maxLength = std::numeric_limits<uint32_t>::max();
  uint32_t numberOfEndBases = std::numeric_limits<uint32_t>::max();
  uint32_t numberOfFowardBases = std::numeric_limits<uint32_t>::max();
  std::string backSeq = "";
  std::string forwardSeq = "";
  trimSeqPars tSeqPars_;
  comparison allowableErrors;
  bool keepOnlyOn = false;

  uint32_t kmerLength = 10;
  uint32_t windowLength = 25;
  uint32_t precision = 10;

  char base = 'N';
  uint32_t qual = 2;

};

class readVecTrimmer {
public:




  template <class T>
  static void trimToMaxLength(std::vector<T> &reads, size_t maxLength);
  template <class T>
  static void trimToMaxLength(T &read, size_t maxLength);
  static void trimToMaxLength(seqInfo &read, size_t maxLength);


  template <class T>
  static void trimAtFirstQualScore(std::vector<T> &reads, const uint32_t qualCutOff);
  template <class T>
  static void trimAtFirstQualScore(T &read, const uint32_t qualCutOff);
  static void trimAtFirstQualScore(seqInfo &read, const uint32_t qualCutOff);

  template <class T>
  static void trimAtFirstBase(std::vector<T> &reads, const char base);
  template <class T>
  static void trimAtFirstBase(T &read, const char base);
  static void trimAtFirstBase(seqInfo &read, const char base);

  template <class T>
  static void trimAtLastBase(std::vector<T> &reads, const char base);
  template <class T>
  static void trimAtLastBase(T &read, const char base);
  static void trimAtLastBase(seqInfo &read, const char base);

  template <class T>
  static void trimAtFirstSeq(std::vector<T> &reads, const std::string & seq);
  template <class T>
  static void trimAtFirstSeq(T &read, const std::string & seq);
  static void trimAtFirstSeq(seqInfo &read, const std::string & seq);

  template <class T>
  static void trimAtLastSeq(std::vector<T> &reads, const std::string & seq);
  template <class T>
  static void trimAtLastSeq(T &read, const std::string & seq);
  static void trimAtLastSeq(seqInfo &read, const std::string & seq);

  template <class T>
  static void trimOffEndBases(std::vector<T> &reads, size_t endBases);
  template <class T>
  static void trimOffEndBases(T &read, size_t endBases);
  static void trimOffEndBases(seqInfo &read, size_t endBases);

  template <class T>
  static void trimOffForwardBases(std::vector<T> &reads, size_t forwardBases);
  template <class T>
  static void trimOffForwardBases(T &read, size_t forwardBases);
  static void trimOffForwardBases(seqInfo &read, size_t forwardBases);

  template <class T>
  static void trimEnds(std::vector<T> &reads, size_t forwardBases, size_t endBases);
  template <class T>
  static void trimEnds(T &read, size_t forwardBases, size_t endBases);
  static void trimEnds(seqInfo &read, size_t forwardBases, size_t endBases);

  template <class T>
  static void trimEndsOfReadsToSharedSeq(std::vector<T> &reads, bool verbose);

	template<class T>
	static void trimAtSequence(std::vector<T> &reads,
			seqInfo &reversePrimer, aligner &alignObj, comparison allowableErrors,
			FullTrimReadsPars::trimSeqPars tSeqPars);
	template<class T>
	static void trimAtSequence(T &read, seqInfo &reversePrimer,
			aligner &alignObj, comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars);
	static void trimAtSequence(seqInfo &read, seqInfo &reversePrimer,
			aligner &alignObj, comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars);

	template<class T>
	static void trimBeforeSequence(std::vector<T> &reads,
			seqInfo &forwardSeq, aligner &alignObj, comparison allowableErrors,
			FullTrimReadsPars::trimSeqPars tSeqPars);
	template<class T>
	static void trimBeforeSequence(T &read, seqInfo &forwardSeq,
			aligner &alignObj, comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars);
	static void trimBeforeSequence(seqInfo &read, seqInfo &forwardSeq,
			aligner &alignObj, comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars);

	template<class T>
	static void trimBetweenSequences(std::vector<T> &reads,
			seqInfo &forwardSeq, seqInfo &backSeq, aligner &alignObj,
			comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars);
	template<class T>
	static void trimBetweenSequences(T &read, seqInfo &forwardSeq,
			seqInfo &backSeq, aligner &alignObj, comparison allowableErrors,
			FullTrimReadsPars::trimSeqPars tSeqPars);
	static void trimBetweenSequences(seqInfo &read, seqInfo &forwardSeq,
			seqInfo &backSeq, aligner &alignObj, comparison allowableErrors,
			FullTrimReadsPars::trimSeqPars tSeqPars);

	template<typename T>
	static void trimToMostCommonKmer(std::vector<T> & seqs,
			FullTrimReadsPars pars,
			aligner &alignObj);

	template<typename T>
	static void trimFromMostCommonKmer(std::vector<T> & seqs,
			FullTrimReadsPars pars,
			aligner &alignObj);

	template<typename T>
	static void trimBetweenMostCommonKmers(std::vector<T> & seqs, FullTrimReadsPars pars,
			aligner &alignObj);




};



template <class T>
void readVecTrimmer::trimToMaxLength(std::vector<T> &reads, size_t maxLength) {
  bib::for_each(reads, [&](T & read){ trimToMaxLength(read, maxLength);} );
  return;
}

template <class T>
void readVecTrimmer::trimToMaxLength(T &read, size_t maxLength) {
	trimToMaxLength(getSeqBase(read), maxLength);
}

template <class T>
void readVecTrimmer::trimAtFirstQualScore(std::vector<T> &reads, const uint32_t qualCutOff){
  bib::for_each(reads, [&](T & read){ trimAtFirstQualScore(read, qualCutOff);} );
  return;
}

template <class T>
void readVecTrimmer::trimAtFirstQualScore(T &read, const uint32_t qualCutOff){
	trimAtFirstQualScore(getSeqBase(read), qualCutOff);
}

template <class T>
void readVecTrimmer::trimAtFirstBase(std::vector<T> &reads, const char base){
  bib::for_each(reads, [&](T & read){ trimAtFirstBase(read, base);} );
  return;
}

template <class T>
void readVecTrimmer::trimAtFirstBase(T &read, const char base){
	trimAtFirstBase(getSeqBase(read), base);
}

template <class T>
void readVecTrimmer::trimAtLastBase(T &read, const char base){
	trimAtLastBase(getSeqBase(read), base);
}

template <class T>
void readVecTrimmer::trimAtLastSeq(T &read, const std::string & seq){
	trimAtLastSeq(getSeqBase(read), seq);
}

template <class T>
void readVecTrimmer::trimAtFirstSeq(T &read, const std::string & seq){
	trimAtFirstSeq(getSeqBase(read), seq);
}



template <class T>
void readVecTrimmer::trimOffEndBases(std::vector<T> &reads, size_t endBases) {
	bib::for_each(reads, [&](T & read){ trimOffEndBases(read, endBases);} );
  return;
}
template <class T>
void readVecTrimmer::trimOffEndBases(T &read, size_t endBases){
	trimOffEndBases(getSeqBase(read), endBases);
}

template <class T>
void readVecTrimmer::trimOffForwardBases(std::vector<T> &reads,
                                         size_t forwardBases) {
	bib::for_each(reads, [&](T & read){ trimOffForwardBases(read, forwardBases);} );
  return;
}
template <class T>
void readVecTrimmer::trimOffForwardBases(T &read, size_t forwardBases) {
	trimOffForwardBases(getSeqBase(read), forwardBases)	;
  return;
}

template<class T>
void readVecTrimmer::trimEnds(std::vector<T> &reads, size_t forwardBases,
		size_t endBases) {
	bib::for_each(reads, [&](T & read) {trimEnds(read, forwardBases, endBases);});
}

template <class T>
void readVecTrimmer::trimEnds(T &read, size_t forwardBases,
		size_t endBases) {
	trimEnds(getSeqBase(read), forwardBases, endBases)	;
  return;
}

template<class T>
void readVecTrimmer::trimAtSequence(std::vector<T> &reads,
		seqInfo &reversePrimer, aligner &alignObj, comparison allowableErrors,
		FullTrimReadsPars::trimSeqPars tSeqPars) {
	bib::for_each(reads,
			[&](T& read) {trimAtSequence(read, reversePrimer, alignObj, allowableErrors, tSeqPars);});
}

template<class T>
void readVecTrimmer::trimAtSequence(T &read, seqInfo &reversePrimer,
		aligner &alignObj, comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars) {
	trimAtSequence(getSeqBase(read), reversePrimer, alignObj,
			allowableErrors, tSeqPars);
}


template<class T>
void readVecTrimmer::trimBeforeSequence(std::vector<T> &reads,
		seqInfo &forwardSeq, aligner &alignObj, comparison allowableErrors,
		FullTrimReadsPars::trimSeqPars tSeqPars) {
	bib::for_each(reads,
			[&](T& read) {trimBeforeSequence(read, forwardSeq, alignObj, allowableErrors, tSeqPars);});
}

template<class T>
void readVecTrimmer::trimBeforeSequence(T &read, seqInfo &forwardSeq,
		aligner &alignObj, comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars) {
	trimBeforeSequence(getSeqBase(read), forwardSeq, alignObj,
			allowableErrors, tSeqPars);
}

template <class T>
void readVecTrimmer::trimBetweenSequences(
		std::vector<T> &reads,
					seqInfo &forwardSeq, seqInfo &backSeq, aligner &alignObj,
					comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars) {
	bib::for_each(reads, [&](T& read){ trimBetweenSequences(read, forwardSeq,backSeq, alignObj, allowableErrors) ;});
  return;
}

template <class T>
void readVecTrimmer::trimBetweenSequences(T &read,
		seqInfo &forwardSeq, seqInfo &backSeq, aligner &alignObj,
							comparison allowableErrors, FullTrimReadsPars::trimSeqPars tSeqPars){
	trimBetweenSequences(getSeqBase(read), forwardSeq, backSeq, alignObj, allowableErrors, tSeqPars);
	return;
}





template <class T>
void readVecTrimmer::trimEndsOfReadsToSharedSeq(std::vector<T> &reads,
                                                bool verbose) {
  VecStr longestSharedSeq = findLongestSharedSeqFromReads(reads);
  if (verbose) {
    std::cout << vectorToString(longestSharedSeq, ", ") << std::endl;
    std::vector<size_t> farthestLocations;
    for (const auto &rIter : reads) {
      std::vector<size_t> farthestLocation =
          findOccurences(rIter.seqBase_.seq_, longestSharedSeq[0]);
      farthestLocations.push_back(vectorMaximum(farthestLocation));
    }
    /*
    std::cout << "Max: " << vectorMaximum(farthestLocations) << std::endl;
    std::cout << "Min: " << vectorMinimum(farthestLocations) << std::endl;
    std::cout << "average: " << vectorMean(farthestLocations) << std::endl;
    std::cout << "median: " << vectorMedian(farthestLocations) << std::endl;
    */
  }
  for (auto &rIter : reads) {
    std::vector<size_t> farthestLocation =
        findOccurences(rIter.seqBase_.seq_, longestSharedSeq[0]);
    size_t maxDistance = vectorMaximum(farthestLocation);
    rIter.trimBack(maxDistance + longestSharedSeq[0].size());
  }
  return;
}


template<typename T>
void readVecTrimmer::trimToMostCommonKmer(std::vector<T> & seqs,
		FullTrimReadsPars pars,
		aligner &alignObj) {

	std::vector<kmerInfo> kInfos;
	probabilityProfile profile(pars.kmerLength);

  for(const auto & readPos : iter::range(len(seqs))){
  	auto & seq = seqs[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(seq) < pars.windowLength ){
  		getSeqBase(seq).on_ = false;
  		continue;
  	}
  	uint32_t startPos = len(seq) - pars.windowLength;
  	uint32_t stopPos = len(seq) - pars.kmerLength + 1;
  	kInfos.emplace_back(kmerInfo(getSeqBase(seq).seq_.substr(startPos), pars.kmerLength, false));
  	for(const auto & pos : iter::range(startPos, stopPos)){
  		profile.add(getSeqBase(seq).seq_.substr(pos, pars.kmerLength), false);
  	}
  }

  std::unordered_map<std::string, uint32_t> kCounts;
  for(const auto & kInfo : kInfos){
  	for(const auto & k : kInfo.kmers_){
  		++kCounts[k.first];
  	}
  }
  profile.updateProfile();
  std::string bestK = "";
  double bestProb = 0;
  uint32_t bestCount = 0;
	for (const auto & kBest : kCounts) {
		if(kBest.second> bestCount){
			bestK = kBest.first;
			bestCount = kBest.second;
			bestProb = profile.getProbabilityOfKmer(kBest.first);
		}else if(kBest.second == bestCount){
			auto prob = profile.getProbabilityOfKmer(kBest.first);
			if(prob > bestProb){
				bestK = kBest.first;
				bestProb = prob;
			}
		}
	}
	seqInfo bestSeq(bestK + ":" + estd::to_string(bestProb), bestK);
  pars.tSeqPars_.within_ = pars.windowLength + 5;
  for(auto & seq : seqs){
  	//skip if the length was too small
  	if(!getSeqBase(seq).on_){
  		continue;
  	}
		readVecTrimmer::trimAtSequence(getSeqBase(seq), bestSeq, alignObj,
				pars.allowableErrors, pars.tSeqPars_);
  }
}

template<typename T>
void readVecTrimmer::trimFromMostCommonKmer(std::vector<T> & seqs,
		FullTrimReadsPars pars,
		aligner &alignObj) {

	std::vector<kmerInfo> kInfos;
	probabilityProfile profile(pars.kmerLength);

  for(const auto & readPos : iter::range(len(seqs))){
  	auto & seq = seqs[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(seq) < pars.windowLength ){
  		getSeqBase(seq).on_ = false;
  		continue;
  	}
  	uint32_t startPos = 0;
  	uint32_t stopPos = pars.windowLength - pars.kmerLength + 1;
  	kInfos.emplace_back(kmerInfo(getSeqBase(seq).seq_.substr(0, stopPos), pars.kmerLength, false));
  	for(const auto & pos : iter::range(startPos, stopPos)){
  		profile.add(getSeqBase(seq).seq_.substr(pos, pars.kmerLength), false);
  	}
  }

  std::unordered_map<std::string, uint32_t> kCounts;
  for(const auto & kInfo : kInfos){
  	for(const auto & k : kInfo.kmers_){
  		++kCounts[k.first];
  	}
  }
  profile.updateProfile();
  std::string bestK = "";
  double bestProb = 0;
  uint32_t bestCount = 0;
	for (const auto & kBest : kCounts) {
		if(kBest.second> bestCount){
			bestK = kBest.first;
			bestCount = kBest.second;
			bestProb = profile.getProbabilityOfKmer(kBest.first);
		}else if(kBest.second == bestCount){
			auto prob = profile.getProbabilityOfKmer(kBest.first);
			if(prob > bestProb){
				bestK = kBest.first;
				bestProb = prob;
			}
		}
	}
	seqInfo bestSeq(bestK + ":" + estd::to_string(bestProb), bestK);
  pars.tSeqPars_.within_ = pars.windowLength + 5;
  for(auto & seq : seqs){
  	//skip if the length was too small
  	if(!getSeqBase(seq).on_){
  		continue;
  	}
		readVecTrimmer::trimBeforeSequence(getSeqBase(seq), bestSeq, alignObj,
				pars.allowableErrors, pars.tSeqPars_);
  }
}

template<typename T>
void readVecTrimmer::trimBetweenMostCommonKmers(std::vector<T> & seqs, FullTrimReadsPars pars,
		aligner &alignObj) {
	trimFromMostCommonKmer(seqs, pars, alignObj);
	bib::for_each(seqs, [](T & seq) {getSeqBase(seq).on_ = true;});
	trimToMostCommonKmer(seqs, pars, alignObj);
}




}  // namespace bib

