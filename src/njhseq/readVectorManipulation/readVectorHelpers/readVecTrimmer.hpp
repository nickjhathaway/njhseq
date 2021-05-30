#pragma once
//
//  trimmer.hpp
//
//  Created by Nicholas Hathaway on 10/22/13.
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

#include "njhseq/alignment.h"
#include "njhseq/utils.h"
#include "njhseq/readVectorManipulation/readVectorOperations.h"
#include "njhseq/objects/seqObjects/readObject.hpp"
#include "njhseq/objects/collapseObjects/opts/IterPar.hpp"
#include "njhseq/objects/seqObjects/seqKmers/seqWithKmerInfo.hpp"
#include "njhseq/objects/helperObjects/probabilityProfile.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers/trimming/trimPars.h"


namespace njhseq {



class readVecTrimmer {
public:



	struct trimCircularGenomeToRefPars{

		seqInfo refSeq_;
		bool doNotReOrientDirection_ = false;
		bool preferHeader_ = false;
		uint32_t kmerLength_= 7;
		uint32_t padding_ = 150;
		uint32_t extend_ = 0;
		uint32_t extendSeqCheckLenTo_ = 100;
		uint32_t extendSeqCheckLenFrom_ = 15;

		bool mark_ = false;
	};

	static std::vector<seqInfo> trimCircularGenomeToRef(seqInfo seq,
			const trimCircularGenomeToRefPars & pars,
			aligner & alignerObj);
  template <class T>
  static std::vector<seqInfo> trimCircularGenomeToRef(std::vector<T> &seqs, const trimCircularGenomeToRefPars & pars,
			aligner & alignerObj);
  template <class T>
  static std::vector<seqInfo> trimCircularGenomeToRef(T &seq, const trimCircularGenomeToRefPars & pars,
			aligner & alignerObj);



	struct TrimEdgesByLowEntropyPars {
		uint32_t windowStep = 5;
		uint32_t windowSize = 50;
		uint32_t kLen = 1;
		double entropyCutOff = 1.5;
		bool mark = false;

	};

	struct TrimEdgesByLowEntropyRes {
		TrimEdgesByLowEntropyRes(uint32_t start, uint32_t end);
		uint32_t start_ { std::numeric_limits<uint32_t>::max() };
		uint32_t end_ { std::numeric_limits<uint32_t>::max() };
	};

	static TrimEdgesByLowEntropyRes determineTrimPostionsByLowEntropy(
			const std::string & seq, const TrimEdgesByLowEntropyPars & pars);


  template <class T>
  static void trimEdgesForLowEntropy(std::vector<T> &seqs, const TrimEdgesByLowEntropyPars&  pars);
  template <class T>
  static void trimEdgesForLowEntropy(T &seq, const TrimEdgesByLowEntropyPars&  pars);
  static void trimEdgesForLowEntropy(seqInfo &seq, const TrimEdgesByLowEntropyPars&  pars);




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
  static void trimToLastQualScore(std::vector<T> &reads, const uint32_t qualCutOff);
  template <class T>
  static void trimToLastQualScore(T &read, const uint32_t qualCutOff);
  static void trimToLastQualScore(seqInfo &read, const uint32_t qualCutOff);

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
  static void trimAtLstripQualScore(std::vector<T> &reads, const uint32_t qualCutOff);
  template <class T>
  static void trimAtLstripQualScore(T &read, const uint32_t qualCutOff);
  static void trimAtLstripQualScore(seqInfo &read, const uint32_t qualCutOff);


  template <class T>
  static void trimAtRstripQualScore(std::vector<T> &reads, const uint32_t qualCutOff);
  template <class T>
  static void trimAtRstripQualScore(T &read, const uint32_t qualCutOff);
  static void trimAtRstripQualScore(seqInfo &read, const uint32_t qualCutOff);

  template <class T>
  static void trimAtLstripBase(std::vector<T> &reads, const char base);
  template <class T>
  static void trimAtLstripBase(T &read, const char base);
  static void trimAtLstripBase(seqInfo &read, const char base);

  template <class T>
  static void trimAtRstripBase(std::vector<T> &reads, const char base);
  template <class T>
  static void trimAtRstripBase(T &read, const char base);
  static void trimAtRstripBase(seqInfo &read, const char base);

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
			const seqInfo &reversePrimer,
			aligner &alignObj,
			const comparison & allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);
	template<class T>
	static void trimAtSequence(T &read,
			const seqInfo &reversePrimer,
			aligner &alignObj,
			const comparison & allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);
	static void trimAtSequence(seqInfo &read,
			const seqInfo &reversePrimer,
			aligner &alignObj,
			const comparison  & allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);

	template<class T>
	static void trimBeforeSequence(std::vector<T> &reads,
			const seqInfo &forwardSeq,
			aligner &alignObj,
			const comparison & allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);
	template<class T>
	static void trimBeforeSequence(T &read,
			const seqInfo &forwardSeq,
			aligner &alignObj,
			const comparison &allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);
	static void trimBeforeSequence(seqInfo &read,
			const seqInfo &forwardSeq,
			aligner &alignObj,
			const comparison  & allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);

	template<class T>
	static void trimBetweenSequences(std::vector<T> &reads,
			const seqInfo &forwardSeq,
			const seqInfo &backSeq,
			aligner &alignObj,
			const comparison & allowableErrors,
			const FullTrimReadsPars::trimSeqPars& tSeqPars);
	template<class T>
	static void trimBetweenSequences(T &read,
			const seqInfo &forwardSeq,
			const seqInfo &backSeq,
			aligner &alignObj,
			const comparison & allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);
	static void trimBetweenSequences(seqInfo &read,
			const seqInfo &forwardSeq,
			const seqInfo &backSeq,
			aligner &alignObj,
			const comparison & allowableErrors,
			const FullTrimReadsPars::trimSeqPars & tSeqPars);

	template<typename T>
	static void trimToMostCommonKmer(std::vector<T> & seqs,
			const FullTrimReadsPars & pars,
			aligner &alignObj);

	template<typename T>
	static void trimFromMostCommonKmer(std::vector<T> & seqs,
			const FullTrimReadsPars & pars,
			aligner &alignObj);

	template<typename T>
	static void trimBetweenMostCommonKmers(std::vector<T> & seqs,
			const FullTrimReadsPars &pars,
			aligner &alignObj);


	template<typename SEQYPTE, typename REFSEQ>
	static void trimSeqToRefByGlobalAln(SEQYPTE & seq, const std::vector<REFSEQ> & refSeqs,
			const std::vector<kmerInfo> & refSeqsKInfos,aligner &alignObj);

	template<typename SEQYPTE, typename REFSEQ>
	static void trimSeqToRefByGlobalAln(std::vector<SEQYPTE> & seq, const std::vector<REFSEQ> & refSeqs,
			const std::vector<kmerInfo> & refSeqsKInfos,aligner &alignerObj);

	template<typename SEQYPTE>
	static void trimSeqToRefByGlobalAln(SEQYPTE &seq,
			const std::vector<seqWithKmerInfo> &refSeqs, aligner &alignObj);

	template<typename SEQYPTE>
	static void trimSeqToRefByGlobalAln(std::vector<SEQYPTE> &seq,
			const std::vector<seqWithKmerInfo> &refSeqs, aligner &alignerObj);

	struct GlobalAlnTrimPars{
		uint32_t startInclusive_ = std::numeric_limits<uint32_t>::max();
		uint32_t endInclusive_ = std::numeric_limits<uint32_t>::max();
		bool needJustOneEnd_ = false;
	};

	template<typename SEQYPTE, typename REFSEQ>
	static void trimSeqToRefByGlobalAln(SEQYPTE & seq,
			const REFSEQ & refSeq,
			const GlobalAlnTrimPars & trimPars,
			aligner &alignerObj);

	template<typename SEQYPTE, typename REFSEQ>
	static void trimSeqToRefByGlobalAln(std::vector<SEQYPTE> & seqs,
			const REFSEQ & refSeq,
			const GlobalAlnTrimPars & trimPars,
			aligner &alignerObj);



	struct BreakUpRes {
		BreakUpRes(const seqInfo & seqBase, uint32_t start, uint32_t end,
				const std::string & pat);
		seqInfo seqBase_;
		uint32_t start_;
		uint32_t end_;
		std::string pat_;
	};

	static std::vector<BreakUpRes> breakUpSeqOnPat(const seqInfo & seq, const std::string & pattern);
	static std::vector<BreakUpRes> breakUpSeqOnPat(const seqInfo & seq, const njh::PatPosFinder & pFinder);

};


template <class T>
std::vector<seqInfo> readVecTrimmer::trimCircularGenomeToRef(std::vector<T> &seqs, const trimCircularGenomeToRefPars & pars,
		aligner & alignerObj){
	std::vector<seqInfo> ret;
	for(const auto & seq : seqs){
		auto res = trimCircularGenomeToRef(getSeqBase(seq), pars, alignerObj);
		addOtherVec(ret, res);
	}
	return ret;
}

template <class T>
std::vector<seqInfo> readVecTrimmer::trimCircularGenomeToRef(T &seq, const trimCircularGenomeToRefPars & pars,
		aligner & alignerObj){
	return trimCircularGenomeToRef(getSeqBase(seq));
}



template<typename SEQYPTE, typename REFSEQ>
void readVecTrimmer::trimSeqToRefByGlobalAln(std::vector<SEQYPTE> & seqs,
		const REFSEQ & refSeq,
		const GlobalAlnTrimPars & trimPars,
		aligner &alignerObj){
	for(auto & seq : seqs){
		trimSeqToRefByGlobalAln(seq, refSeq, trimPars, alignerObj);
	}
}

template<typename SEQYPTE, typename REFSEQ>
void readVecTrimmer::trimSeqToRefByGlobalAln(SEQYPTE & seq,
		const REFSEQ & refSeq,
		const GlobalAlnTrimPars & gtrimPars,
		aligner &alignerObj){
	auto startInclusive = gtrimPars.startInclusive_;
	auto endIncludsive = gtrimPars.endInclusive_;
	if(std::numeric_limits<uint32_t>::max() == startInclusive){
		startInclusive = 0;
	}
	if(std::numeric_limits<uint32_t>::max() == endIncludsive){
		endIncludsive = len(getSeqBase(refSeq)) - 1;
	}
	if (startInclusive>= getSeqBase(refSeq).seq_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << "error startInclusive, " << startInclusive
				<< ", greater than ref " << getSeqBase(refSeq).name_ << " length: "
				<< getSeqBase(refSeq).seq_.size() << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (endIncludsive >= getSeqBase(refSeq).seq_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << "error endInclusvie, " << endIncludsive
				<< ", greater than ref " << getSeqBase(refSeq).name_ << " length: "
				<< getSeqBase(refSeq).seq_.size() << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (endIncludsive < startInclusive) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << "error endInclusvie, " << endIncludsive
				<< ", less than or equal to startInclusive,  " << startInclusive
				<< "\n";
		throw std::runtime_error { ss.str() };
	}

	alignerObj.alignCacheGlobal(getSeqBase(refSeq), getSeqBase(seq));

	//ref positions
	uint32_t refAlnStartPos = alignerObj.getAlignPosForSeqAPos(startInclusive);
	uint32_t refAlnStopPos = alignerObj.getAlignPosForSeqAPos(endIncludsive);
	//aligned bases
	uint32_t firstAlignedBase = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
	uint32_t lastAlignedBase = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
	//front
	uint32_t trimFront = std::numeric_limits<uint32_t>::max();
	if (refAlnStartPos >= firstAlignedBase && refAlnStartPos < lastAlignedBase) {
		trimFront = alignerObj.getSeqPosForAlnBPos(refAlnStartPos);
	} else {
		getSeqBase(seq).on_ = false;
	}
	//back
	uint32_t trimBack = std::numeric_limits<uint32_t>::max();
	if (refAlnStopPos > firstAlignedBase && refAlnStopPos <= lastAlignedBase) {
		trimBack = alignerObj.getSeqPosForAlnBPos(refAlnStopPos);
	} else {
		getSeqBase(seq).on_ = false;
	}


	if (std::numeric_limits<uint32_t>::max() != trimFront
			&& std::numeric_limits<uint32_t>::max() != trimBack) {
		getSeqBase(seq).setClip(trimFront, trimBack);
	} else if (std::numeric_limits<uint32_t>::max() != trimFront) {
		getSeqBase(seq).trimFront(trimFront);
		if(gtrimPars.needJustOneEnd_ ){
			getSeqBase(seq).on_ =  true;
		}
	} else if (std::numeric_limits<uint32_t>::max() != trimBack) {
		getSeqBase(seq).trimBack(trimBack + 1);
		if(gtrimPars.needJustOneEnd_){
			getSeqBase(seq).on_ =  true;
		}
	}
}


template<typename SEQYPTE>
void readVecTrimmer::trimSeqToRefByGlobalAln(SEQYPTE &seq,
		const std::vector<seqWithKmerInfo> &refSeqs, aligner &alignerObj){

	uint32_t bestIndex = std::numeric_limits<uint32_t>::max();
	if (1 == refSeqs.size()) {
		bestIndex = 0;
	} else {
		double kmerCutOff = 0.8;
		kmerInfo seqKInfo(getSeqBase(seq).seq_, refSeqs.front().kInfo_.kLen_, false);
	  double bestScore = std::numeric_limits<double>::lowest();
	  std::vector<uint32_t> bestRefs;
	  while(bestRefs.empty()){
	    for (const auto refPos : iter::range(refSeqs.size())) {
	      const auto & ref = refSeqs[refPos];
	      if(refSeqs[refPos].kInfo_.compareKmers(seqKInfo).second < kmerCutOff){
	       	continue;
	      }
	      alignerObj.alignCacheGlobal(ref, getSeqBase(seq));
				double currentScore = 0;
				alignerObj.profileAlignment(ref, getSeqBase(seq), false, true, false);
				currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
				if (currentScore == bestScore) {
					bestRefs.push_back(refPos);
				}
				if (currentScore > bestScore) {
					bestRefs.clear();
					bestRefs.push_back(refPos);
					bestScore = currentScore;
				}
			}
	    if(kmerCutOff < 0 && bestRefs.empty()){
	    		std::stringstream ss;
	    		ss << __PRETTY_FUNCTION__ << ", error, couldn't find a matching ref for : " << getSeqBase(seq).name_ << "\n";
	    		throw std::runtime_error{ss.str()};
	    }
	    kmerCutOff -= .1;
	  }
	  bestIndex = bestRefs.front();
	}
	alignerObj.alignCacheGlobal(refSeqs[bestIndex], getSeqBase(seq));

	//getTrimFront
	uint32_t trimFront = std::numeric_limits<uint32_t>::max();
	if('-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
		 '-' == alignerObj.alignObjectA_.seqBase_.seq_.front() ){
		auto refAlnPos = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
		trimFront = alignerObj.getSeqPosForAlnBPos(refAlnPos);
	}else{
		if('-' == alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() ){
			getSeqBase(seq).on_ = false;
		}
	}

	//getTrimBack
	uint32_t trimBack = std::numeric_limits<uint32_t>::max();
	if('-' != alignerObj.alignObjectB_.seqBase_.seq_.back() &&
		 '-' == alignerObj.alignObjectA_.seqBase_.seq_.back() ){
		auto refAlnPos = alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
		trimBack = alignerObj.getSeqPosForAlnBPos(refAlnPos);
	}else{
		if('-' == alignerObj.alignObjectB_.seqBase_.seq_.back() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back() ){
			getSeqBase(seq).on_ = false;
		}
	}
	if (std::numeric_limits<uint32_t>::max() != trimFront
			&& std::numeric_limits<uint32_t>::max() != trimBack) {
		getSeqBase(seq).setClip(trimFront, trimBack);
	} else if (std::numeric_limits<uint32_t>::max() != trimFront) {
		getSeqBase(seq).trimFront(trimFront);
	} else if (std::numeric_limits<uint32_t>::max() != trimBack) {
		getSeqBase(seq).trimBack(trimBack + 1);
	}
}

template<typename SEQYPTE>
void readVecTrimmer::trimSeqToRefByGlobalAln(std::vector<SEQYPTE> &seqs,
		const std::vector<seqWithKmerInfo> &refSeqs, aligner &alignerObj){
	for(auto & seq : seqs	){
		trimSeqToRefByGlobalAln(seq, refSeqs, alignerObj);
	}
}



template<typename SEQYPTE, typename REFSEQ>
void readVecTrimmer::trimSeqToRefByGlobalAln(std::vector<SEQYPTE> & seqs,
		const std::vector<REFSEQ> & refSeqs,
		const std::vector<kmerInfo> & refSeqsKInfos,
		aligner &alignerObj){
	for(auto & seq : seqs	){
		trimSeqToRefByGlobalAln(seq, refSeqs, refSeqsKInfos, alignerObj);
	}
}

template<typename SEQYPTE, typename REFSEQ>
void readVecTrimmer::trimSeqToRefByGlobalAln(SEQYPTE & seq,
		const std::vector<REFSEQ> & refSeqs,
		const std::vector<kmerInfo> & refSeqsKInfos,
		aligner &alignerObj){

	uint32_t bestIndex = std::numeric_limits<uint32_t>::max();
	if (1 == refSeqs.size()) {
		bestIndex = 0;
	} else {
		double kmerCutOff = 0.8;
		kmerInfo seqKInfo(getSeqBase(seq).seq_, refSeqsKInfos.front().kLen_, false);
	  double bestScore = std::numeric_limits<double>::lowest();
	  std::vector<uint32_t> bestRefs;
	  while(bestRefs.empty()){
	    for (const auto refPos : iter::range(refSeqs.size())) {
	      const auto & ref = refSeqs[refPos];
	      if(refSeqsKInfos[refPos].compareKmers(seqKInfo).second < kmerCutOff){
	       	continue;
	      }
	      alignerObj.alignCacheGlobal(ref, getSeqBase(seq));
				double currentScore = 0;
				alignerObj.profileAlignment(ref, getSeqBase(seq), false, true, false);
				currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
				if (currentScore == bestScore) {
					bestRefs.push_back(refPos);
				}
				if (currentScore > bestScore) {
					bestRefs.clear();
					bestRefs.push_back(refPos);
					bestScore = currentScore;
				}
			}
	    if(kmerCutOff < 0 && bestRefs.empty()){
	    		std::stringstream ss;
	    		ss << __PRETTY_FUNCTION__ << ", error, couldn't find a matching ref for : " << getSeqBase(seq).name_ << "\n";
	    		throw std::runtime_error{ss.str()};
	    }
	    kmerCutOff -= .1;
	  }
	  bestIndex = bestRefs.front();
	}
	alignerObj.alignCacheGlobal(refSeqs[bestIndex], getSeqBase(seq));

	//getTrimFront
	uint32_t trimFront = std::numeric_limits<uint32_t>::max();
	if('-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
		 '-' == alignerObj.alignObjectA_.seqBase_.seq_.front() ){
		auto refAlnPos = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
		trimFront = alignerObj.getSeqPosForAlnBPos(refAlnPos);
	}else{
		if('-' == alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() ){
			getSeqBase(seq).on_ = false;
		}
	}

	//getTrimBack
	uint32_t trimBack = std::numeric_limits<uint32_t>::max();
	if('-' != alignerObj.alignObjectB_.seqBase_.seq_.back() &&
		 '-' == alignerObj.alignObjectA_.seqBase_.seq_.back() ){
		auto refAlnPos = alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
		trimBack = alignerObj.getSeqPosForAlnBPos(refAlnPos);
	}else{
		if('-' == alignerObj.alignObjectB_.seqBase_.seq_.back() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back() ){
			getSeqBase(seq).on_ = false;
		}
	}
	if (std::numeric_limits<uint32_t>::max() != trimFront
			&& std::numeric_limits<uint32_t>::max() != trimBack) {
		getSeqBase(seq).setClip(trimFront, trimBack);
	} else if (std::numeric_limits<uint32_t>::max() != trimFront) {
		getSeqBase(seq).trimFront(trimFront);
	} else if (std::numeric_limits<uint32_t>::max() != trimBack) {
		getSeqBase(seq).trimBack(trimBack + 1);
	}
}



template <class T>
void readVecTrimmer::trimEdgesForLowEntropy(std::vector<T> &seqs, const TrimEdgesByLowEntropyPars&  pars){
  njh::for_each(seqs, [&](T & seq){ trimEdgesForLowEntropy(seq, pars);} );
  return;
}

template <class T>
void readVecTrimmer::trimEdgesForLowEntropy(T &seq, const TrimEdgesByLowEntropyPars&  pars){
	trimEdgesForLowEntropy(getSeqBase(read), pars);
}



template <class T>
void readVecTrimmer::trimToMaxLength(std::vector<T> &reads, size_t maxLength) {
  njh::for_each(reads, [&](T & read){ trimToMaxLength(read, maxLength);} );
  return;
}

template <class T>
void readVecTrimmer::trimToMaxLength(T &read, size_t maxLength) {
	trimToMaxLength(getSeqBase(read), maxLength);
}

template <class T>
void readVecTrimmer::trimAtFirstQualScore(std::vector<T> &reads, const uint32_t qualCutOff){
  njh::for_each(reads, [&](T & read){ trimAtFirstQualScore(read, qualCutOff);} );
  return;
}

template <class T>
void readVecTrimmer::trimToLastQualScore(std::vector<T> &reads, const uint32_t qualCutOff){
  njh::for_each(reads, [&](T & read){ trimToLastQualScore(read, qualCutOff);} );
  return;
}



template <class T>
void readVecTrimmer::trimAtFirstQualScore(T &read, const uint32_t qualCutOff){
	trimAtFirstQualScore(getSeqBase(read), qualCutOff);
}

template <class T>
void readVecTrimmer::trimToLastQualScore(T &read, const uint32_t qualCutOff){
	trimToLastQualScore(getSeqBase(read), qualCutOff);
}


template <class T>
void readVecTrimmer::trimAtFirstBase(std::vector<T> &reads, const char base){
  njh::for_each(reads, [&](T & read){ trimAtFirstBase(read, base);} );
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
void readVecTrimmer::trimAtLstripQualScore(std::vector<T> &reads, const uint32_t qualCutOff){
  njh::for_each(reads, [&](T & read){ trimAtLstripQualScore(read, qualCutOff);} );
  return;
}

template <class T>
void readVecTrimmer::trimAtLstripQualScore(T &read, const uint32_t qualCutOff){
	trimAtLstripQualScore(getSeqBase(read), qualCutOff);
}

template <class T>
void readVecTrimmer::trimAtRstripQualScore(std::vector<T> &reads, const uint32_t qualCutOff){
  njh::for_each(reads, [&](T & read){ trimAtRstripQualScore(read, qualCutOff);} );
  return;
}

template <class T>
void readVecTrimmer::trimAtRstripQualScore(T &read, const uint32_t qualCutOff){
	trimAtRstripQualScore(getSeqBase(read), qualCutOff);
}

template <class T>
void readVecTrimmer::trimAtLstripBase(std::vector<T> &reads, const char base){
  njh::for_each(reads, [&](T & read){ trimAtLstripBase(read, base);} );
  return;
}


template <class T>
void readVecTrimmer::trimAtLstripBase(T &read, const char base){
	trimAtLstripBase(getSeqBase(read), base);
}

template <class T>
void readVecTrimmer::trimAtRstripBase(std::vector<T> &reads, const char base){
  njh::for_each(reads, [&](T & read){ trimAtRstripBase(read, base);} );
  return;
}


template <class T>
void readVecTrimmer::trimAtRstripBase(T &read, const char base){
	trimAtRstripBase(getSeqBase(read), base);
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
	njh::for_each(reads, [&](T & read){ trimOffEndBases(read, endBases);} );
  return;
}
template <class T>
void readVecTrimmer::trimOffEndBases(T &read, size_t endBases){
	trimOffEndBases(getSeqBase(read), endBases);
}

template <class T>
void readVecTrimmer::trimOffForwardBases(std::vector<T> &reads,
                                         size_t forwardBases) {
	njh::for_each(reads, [&](T & read){ trimOffForwardBases(read, forwardBases);} );
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
	njh::for_each(reads, [&](T & read) {trimEnds(read, forwardBases, endBases);});
}

template <class T>
void readVecTrimmer::trimEnds(T &read, size_t forwardBases,
		size_t endBases) {
	trimEnds(getSeqBase(read), forwardBases, endBases)	;
  return;
}

template<class T>
void readVecTrimmer::trimAtSequence(std::vector<T> &reads,
		const seqInfo &reversePrimer, aligner &alignObj, const comparison & allowableErrors,
		const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	njh::for_each(reads,
			[&](T& read) {trimAtSequence(read, reversePrimer, alignObj, allowableErrors, tSeqPars);});
}

template<class T>
void readVecTrimmer::trimAtSequence(T &read, const seqInfo &reversePrimer,
		aligner &alignObj, const comparison & allowableErrors, const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	trimAtSequence(getSeqBase(read), reversePrimer, alignObj,
			allowableErrors, tSeqPars);
}


template<class T>
void readVecTrimmer::trimBeforeSequence(std::vector<T> &reads,
		const seqInfo &forwardSeq, aligner &alignObj, const comparison & allowableErrors,
		const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	njh::for_each(reads,
			[&](T& read) {trimBeforeSequence(read, forwardSeq, alignObj, allowableErrors, tSeqPars);});
}

template<class T>
void readVecTrimmer::trimBeforeSequence(T &read, const seqInfo &forwardSeq,
		aligner &alignObj, const comparison & allowableErrors, const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	trimBeforeSequence(getSeqBase(read), forwardSeq, alignObj,
			allowableErrors, tSeqPars);
}

template <class T>
void readVecTrimmer::trimBetweenSequences(
		std::vector<T> &reads,
					const seqInfo &forwardSeq, const seqInfo &backSeq, aligner &alignObj,
					const comparison & allowableErrors, const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	njh::for_each(reads, [&](T& read){ trimBetweenSequences(read, forwardSeq,backSeq, alignObj, allowableErrors) ;});
  return;
}

template <class T>
void readVecTrimmer::trimBetweenSequences(T &read,
		const seqInfo &forwardSeq, const seqInfo &backSeq, aligner &alignObj,
							const comparison & allowableErrors, const FullTrimReadsPars::trimSeqPars & tSeqPars){
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
		const FullTrimReadsPars  & pars,
		aligner &alignObj) {

	std::vector<kmerInfo> kInfos;
	probabilityProfile profile(pars.kmerLength);

  for(const auto readPos : iter::range(len(seqs))){
  	auto & seq = seqs[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(seq) < pars.windowLength ){
  		getSeqBase(seq).on_ = false;
  		continue;
  	}
  	uint32_t startPos = len(seq) - pars.windowLength;
  	uint32_t stopPos = len(seq) - pars.kmerLength + 1;
  	kInfos.emplace_back(kmerInfo(getSeqBase(seq).seq_.substr(startPos), pars.kmerLength, false));
  	for(const auto pos : iter::range(startPos, stopPos)){
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

	auto parsCopy = pars.tSeqPars_;
	parsCopy.within_ = pars.windowLength + 5;
	for (auto & seq : seqs) {
		//skip if the length was too small
		if (!getSeqBase(seq).on_) {
			continue;
		}
		readVecTrimmer::trimAtSequence(getSeqBase(seq), bestSeq, alignObj,
				pars.allowableErrors, parsCopy);
	}
}

template<typename T>
void readVecTrimmer::trimFromMostCommonKmer(std::vector<T> & seqs,
		const FullTrimReadsPars & pars,
		aligner &alignObj) {

	std::vector<kmerInfo> kInfos;
	probabilityProfile profile(pars.kmerLength);

  for(const auto readPos : iter::range(len(seqs))){
  	auto & seq = seqs[readPos];
  	//first turn off short sequences if they are too short and then continue on
  	if(len(seq) < pars.windowLength ){
  		getSeqBase(seq).on_ = false;
  		continue;
  	}
  	uint32_t startPos = 0;
  	uint32_t stopPos = pars.windowLength - pars.kmerLength + 1;
  	kInfos.emplace_back(kmerInfo(getSeqBase(seq).seq_.substr(0, stopPos), pars.kmerLength, false));
  	for(const auto pos : iter::range(startPos, stopPos)){
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
	auto parsCopy = pars.tSeqPars_;
	parsCopy.within_ = pars.windowLength + 5;
  for(auto & seq : seqs){
  	//skip if the length was too small
  	if(!getSeqBase(seq).on_){
  		continue;
  	}
		readVecTrimmer::trimBeforeSequence(getSeqBase(seq), bestSeq, alignObj,
				pars.allowableErrors, parsCopy);
  }
}

template<typename T>
void readVecTrimmer::trimBetweenMostCommonKmers(std::vector<T> & seqs, const FullTrimReadsPars & pars,
		aligner &alignObj) {
	trimFromMostCommonKmer(seqs, pars, alignObj);
	njh::for_each(seqs, [](T & seq) {getSeqBase(seq).on_ = true;});
	trimToMostCommonKmer(seqs, pars, alignObj);
}




}  // namespace njh

