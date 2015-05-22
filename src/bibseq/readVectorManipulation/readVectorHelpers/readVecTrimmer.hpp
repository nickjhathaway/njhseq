#pragma once
//
//  trimmer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/22/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/alignment.h"
#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorOperations.h"

namespace bibseq {

class readVecTrimmer {

 public:
  template <class T>
  static void trimToMaxLength(std::vector<T> &reads, size_t maxLength);
  template <class T>
  static void trimToMaxLength(T &read, size_t maxLength);

  template <class T>
  static void trimOffEndBases(std::vector<T> &reads, size_t endBases);
  template <class T>
  static void trimOffEndBases(T &read, size_t endBases);

  template <class T>
  static void trimOffForwardBases(std::vector<T> &reads, size_t forwardBases);
  template <class T>
  static void trimOffForwardBases(T &read, size_t forwardBases);

  template <class T>
  static void trimEndsOfReadsToSharedSeq(std::vector<T> &reads, bool verbose);

  template <class T>
  static void trimAtSequence(std::vector<T> &reads, readObject &reversePrimer,
                             runningParameters runParmas, aligner &alignObj,
                             bool includeSequence, bool sequenceToLowerCase,
                             bool weighHomopolyer, bool removePreviousSameBases);
  template <class T>
  static void trimBeforeSequence(std::vector<T> &reads, readObject &forwardSeq,
                                 runningParameters runParmas, aligner &alignObj,
                                 bool includeSequence, bool sequenceToLowerCase,
                                 bool weighHomopolyer, bool removePreviousSameBases);
  template <class T>
  static void trimBetweenSequences(std::vector<T> &reads,
                                   readObject &forwardSeq, readObject &backSeq,
                                   runningParameters runParams,
                                   aligner &alignObj, bool includeSequence,
                                   bool sequenceToLowerCase,
                                   bool weighHomopolyer, bool removePreviousSameBases);


  template <class T>
  static void trimAtSequenceIdentity(std::vector<T> &reads, readObject &reversePrimer,
                             runningParameters runParmas, aligner &alignObj,
                             bool includeSequence, bool sequenceToLowerCase,
                             bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff);
  template <class T>
  static void trimAtSequenceIdentity(T &read, readObject &reversePrimer,
                             runningParameters runParmas, aligner &alignObj,
                             bool includeSequence, bool sequenceToLowerCase,
                             bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff);
  template <class T>
  static void trimBetweenSequencesIdentity(T &read,
                                   readObject &forwardSeq, readObject &backSeq,
                                   runningParameters runParams,
                                   aligner &alignObj, bool includeSequence,
                                   bool sequenceToLowerCase,
                                   bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff);

  template <class T>
  static void trimBeforeSequenceIdentity(std::vector<T> &reads, readObject &forwardSeq,
                                 runningParameters runParmas, aligner &alignObj,
                                 bool includeSequence, bool sequenceToLowerCase,
                                 bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff);
  template <class T>
  static void trimBeforeSequenceIdentity(T &reads, readObject &forwardSeq,
                                 runningParameters runParmas, aligner &alignObj,
                                 bool includeSequence, bool sequenceToLowerCase,
                                 bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff);

  template <class T>
  static void trimBetweenSequencesIdentity(std::vector<T> &reads,
                                   readObject &forwardSeq, readObject &backSeq,
                                   runningParameters runParams,
                                   aligner &alignObj, bool includeSequence,
                                   bool sequenceToLowerCase,
                                   bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff);
};

template <class T>
void readVecTrimmer::trimToMaxLength(std::vector<T> &reads, size_t maxLength) {
  for_each(reads, [&](T & read){ trimToMaxLength(read, maxLength);} );
  return;
}

template <class T>
void readVecTrimmer::trimToMaxLength(T &read, size_t maxLength) {
  if (read.seqBase_.seq_.size() > maxLength) {
    read.setClip(0, maxLength);
  }
  return;
}

template <class T>
void readVecTrimmer::trimOffEndBases(std::vector<T> &reads, size_t endBases) {
	for_each(reads, [&](T & read){ trimOffEndBases(read, endBases);} );
  return;
}
template <class T>
void readVecTrimmer::trimOffEndBases(T &read, size_t endBases){
	read.trimBack(read.seqBase_.seq_.size() - endBases);
}

template <class T>
void readVecTrimmer::trimOffForwardBases(std::vector<T> &reads,
                                         size_t forwardBases) {
	for_each(reads, [&](T & read){ trimOffForwardBases(read, forwardBases);} );
  return;
}
template <class T>
void readVecTrimmer::trimOffForwardBases(T &read, size_t forwardBases) {
	 read.trimFront(forwardBases);
  return;
}

template <class T>
void readVecTrimmer::trimAtSequence(std::vector<T> &reads,
                                    readObject &reversePrimer,
                                    runningParameters runParmas,
                                    aligner &alignObj, bool includeSequence,
                                    bool sequenceToLowerCase,
                                    bool weighHomopolyer, bool removePreviousSameBases) {
	for_each(reads, [&](T& read){ trimAtSequenceIdentity(read, reversePrimer, runParmas, alignObj, includeSequence, sequenceToLowerCase,weighHomopolyer, removePreviousSameBases, 0.5) ;});
  return;
}
template <class T>
void readVecTrimmer::trimAtSequenceIdentity(T &read,
                                    readObject &reversePrimer,
                                    runningParameters runParmas,
                                    aligner &alignObj, bool includeSequence,
                                    bool sequenceToLowerCase,
                                    bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff) {
	 	std::pair<int, int> rPos = alignObj.findReversePrimer(read, reversePrimer);
	 	alignObj.rearrangeObjs(read.seqBase_, reversePrimer.seqBase_, true);

		alignObj.profilePrimerAlignment(read, reversePrimer, weighHomopolyer);
		bool passInspection = runParmas.errors_.passErrorProfile(alignObj.errors_);
		if(alignObj.distances_.queryCoverage_ < queryCutOff){
			passInspection = false;
		}
		read.remove = !passInspection;
		/*
		alignObj.profilePrimerAlignment(read, reversePrimer, weighHomopolyer);
		bool totalyOk = true;
		if (alignObj.numberOfLargeGaps > runParmas.errors.largeBaseIndel ||
				alignObj.numberOfOneIndel > runParmas.errors.oneBaseIndel ||
				alignObj.numberOfTwoIndel > runParmas.errors.twoBaseIndel) {
			totalyOk = false;
		}
		// change this if we should be doing full reverse primer, the next part
		// states that if at least half of the primer is found then it's ok, i think
		// that's ok for now; nick 09.17.2013
		// if (reversePrimer.seqBase_.seq_.size()>alignObj.score) {

		if (alignObj.alignObjectA.seqBase_.seq_.size() > alignObj.score) {
			totalyOk = false;
		}
		// std::cout<<"totallyOk: "<<convertBoolToString(totalyOk)<<std::endl;
		if (alignObj.alignObjectA.seqBase_.seq_.size() <
		reversePrimer.seqBase_.seq_.size() / 2) {
			totalyOk = false;
		}
		read.remove = !totalyOk;
		 */
		if (includeSequence) {
			read.setClip(0, rPos.second);
			if (sequenceToLowerCase) {
				if(removePreviousSameBases){
					while (read.seqBase_.seq_[rPos.first] ==
								 read.seqBase_.seq_[rPos.first - 1]) {
						--rPos.first;
					}
				}
				changeSubStrToLowerToEnd(read.seqBase_.seq_, rPos.first);
			}
		} else {
			if(removePreviousSameBases){
				bool foundPrevious = false;
				while (read.seqBase_.seq_[rPos.first] ==
										 read.seqBase_.seq_[rPos.first - 1]) {
					foundPrevious = true;
					--rPos.first;
				}

				//set clip will keep the second postion you give it so you have to minus one more if you found any
				if(foundPrevious){
					--rPos.first;
				}
			}
			read.setClip(0, rPos.first - 1);
		}
  return;
}

template <class T>
void readVecTrimmer::trimAtSequenceIdentity(std::vector<T> &reads,
                                    readObject &reversePrimer,
                                    runningParameters runParmas,
                                    aligner &alignObj, bool includeSequence,
                                    bool sequenceToLowerCase,
                                    bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff) {
	for_each(reads, [&](T& read){ trimAtSequenceIdentity(read, reversePrimer, runParmas, alignObj, includeSequence,
			sequenceToLowerCase,weighHomopolyer, removePreviousSameBases, queryCutOff) ;});
  return;
}

template <class T>
void readVecTrimmer::trimBeforeSequenceIdentity(T &read, readObject &forwardSeq,
                               runningParameters runParmas, aligner &alignObj,
                               bool includeSequence, bool sequenceToLowerCase,
                               bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff){
  std::pair<int, int> fPos = alignObj.findReversePrimer(read, forwardSeq);
  alignObj.rearrangeObjs(read.seqBase_, forwardSeq.seqBase_, true);
  //need to reordered and profile aligner
  alignObj.profilePrimerAlignment(read, forwardSeq, weighHomopolyer);
	bool passInspection = runParmas.errors_.passErrorProfile(alignObj.errors_);
	if(alignObj.distances_.queryCoverage_ < queryCutOff){
		passInspection = false;
	}
	read.remove = !passInspection;
  if (includeSequence) {
    read.trimFront(fPos.first);
    if (sequenceToLowerCase) {
    	if(removePreviousSameBases){
        while (read.seqBase_.seq_[fPos.second] ==
               read.seqBase_.seq_[fPos.second + 1]) {
          ++fPos.second;
        }
    	}
      changeSubStrToLowerFromBegining(read.seqBase_.seq_, fPos.second);
    }
  } else {
  	if(removePreviousSameBases){
      while (read.seqBase_.seq_[fPos.second] ==
             read.seqBase_.seq_[fPos.second + 1]) {
        ++fPos.second;
      }
  	}
    read.trimFront(fPos.second + 1);
  }
}
template <class T>
void readVecTrimmer::trimBetweenSequencesIdentity(T &read,
                                 readObject &forwardSeq, readObject &backSeq,
                                 runningParameters runParams,
                                 aligner &alignObj, bool includeSequence,
                                 bool sequenceToLowerCase,
                                 bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff){
	trimAtSequenceIdentity(read, backSeq, runParams, alignObj, includeSequence, sequenceToLowerCase, weighHomopolyer, removePreviousSameBases, queryCutOff);
	trimBeforeSequenceIdentity(read, forwardSeq, runParams, alignObj, includeSequence, sequenceToLowerCase, weighHomopolyer, removePreviousSameBases, queryCutOff);
	return;
}
template <class T>
void readVecTrimmer::trimBeforeSequenceIdentity(std::vector<T> &reads, readObject &forwardSeq,
                               runningParameters runParmas, aligner &alignObj,
                               bool includeSequence, bool sequenceToLowerCase,
                               bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff){
	for_each(reads, [&](T& read){ trimBeforeSequenceIdentity(read, forwardSeq, runParmas, alignObj, includeSequence,
			sequenceToLowerCase,weighHomopolyer, removePreviousSameBases, queryCutOff) ;});
	return;
}
template <class T>
void readVecTrimmer::trimBetweenSequencesIdentity(std::vector<T> &reads,
                                 readObject &forwardSeq, readObject &backSeq,
                                 runningParameters runParams,
                                 aligner &alignObj, bool includeSequence,
                                 bool sequenceToLowerCase,
                                 bool weighHomopolyer, bool removePreviousSameBases, double queryCutOff){
	for_each(reads, [&](T& read){ trimBetweenSequencesIdentity(read, forwardSeq,backSeq, runParams, alignObj, includeSequence,
			sequenceToLowerCase,weighHomopolyer, removePreviousSameBases, queryCutOff) ;});
	return;
}


template <class T>
void readVecTrimmer::trimBeforeSequence(std::vector<T> &reads,
                                        readObject &forwardSeq,
                                        runningParameters runParmas,
                                        aligner &alignObj, bool includeSequence,
                                        bool sequenceToLowerCase,
                                        bool weighHomopolyer, bool removePreviousSameBases) {
	for_each(reads, [&](T& read){ trimBeforeSequenceIdentity(read, forwardSeq, runParmas, alignObj, includeSequence,
			sequenceToLowerCase,weighHomopolyer, removePreviousSameBases, 0.5) ;});
  return;
}

template <class T>
void readVecTrimmer::trimBetweenSequences(
    std::vector<T> &reads, readObject &forwardSeq, readObject &backSeq,
    runningParameters runParams, aligner &alignObj, bool includeSequence,
    bool sequenceToLowerCase, bool weighHomopolyer, bool removePreviousSameBases) {
	for_each(reads, [&](T& read){ trimBetweenSequencesIdentity(read, forwardSeq,backSeq, runParams, alignObj, includeSequence,
			sequenceToLowerCase, weighHomopolyer, removePreviousSameBases, 0.5) ;});
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
    std::cout << "Max: " << vectorMaximum(farthestLocations) << std::endl;
    std::cout << "Min: " << vectorMinimum(farthestLocations) << std::endl;
    std::cout << "average: " << vectorMean(farthestLocations) << std::endl;
    std::cout << "median: " << vectorMedian(farthestLocations) << std::endl;
  }
  for (auto &rIter : reads) {
    std::vector<size_t> farthestLocation =
        findOccurences(rIter.seqBase_.seq_, longestSharedSeq[0]);
    size_t maxDistance = vectorMaximum(farthestLocation);
    rIter.trimBack(maxDistance + longestSharedSeq[0].size());
  }
  return;
}
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "readVecTrimmer.cpp"
#endif