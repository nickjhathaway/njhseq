#pragma once
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
//  seqToolsUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 4/27/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/seqObjects/seqWithKmerInfo.hpp"
#include "bibseq/objects/seqObjects/sffObject.hpp"

#include "bibseq/seqToolsUtils/aminoAcidInfo.hpp"
#include "bibseq/alignment.h"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/utils.h"
#include "bibseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp"
#include "bibseq/simulation/simulationCommon.hpp"

#include "bibseq/simulation/errorProfile.hpp"
#include "bibseq/objects/helperObjects/probabilityProfile.hpp"

#include <bibcpp/graphics/colorUtils.hpp>


namespace bibseq {


template<typename T>
uint32_t getMismatches(const T & read1,
				const T & read2,
				aligner alignerObj, bool weightHomopolymers){
	alignerObj.alignVec(read1.seqBase_,read2.seqBase_, false);
	alignerObj.profilePrimerAlignment(read1.seqBase_, read2.seqBase_, weightHomopolymers);
	return alignerObj.comp_.hqMismatches_;
};


Json::Value genMinTreeData(const std::vector<readObject> & reads);
Json::Value genMinTreeData(const std::vector<readObject> & reads, aligner & alignerObj);



uint32_t countSeqs(const readObjectIOOptions & opts, bool verbose);

void setUpSampleDirs(
    const std::string& sampleNamesFilename,
		const std::string& mainDirectoryName,
		bool separatedDirs);
std::pair<uint32_t,uint32_t> getMaximumRounds(double cnt);
std::vector<char> getAAs(bool noStopCodon);
class cluster;
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
template<typename T>
simulation::mismatchProfile simProfileWithReads(const std::vector<T> & reads, aligner & alignerObj,
		randomGenerator & gen, uint32_t sampNum, bool local, bool ignoreGaps){
  std::vector<uint32_t> randomFirstPos = gen.unifRandVector<uint32_t>(0, len(reads),sampNum);
  std::vector<uint32_t> randomSecondPos = gen.unifRandVector<uint32_t>(0, len(reads),sampNum);
  simulation::mismatchProfile ret;
  for(const auto & pos : iter::range(sampNum)){
  	alignerObj.alignVec(reads[randomFirstPos[pos]],
  			reads[randomSecondPos[pos]], local);
  	ret.increaseCountAmount(alignerObj.alignObjectA_.seqBase_.seq_,
  	  				alignerObj.alignObjectB_.seqBase_.seq_, 1);
  }
  ret.setFractions();
  return ret;
}



std::string genHtmlStrForPsuedoMintree(std::string jsonFileName);




template <typename T>
void findForwardPrimerAndSort(
    std::vector<T>& inputReads,
    const std::map<std::string, std::pair<std::string, std::string>>& primers,
    std::map<std::string, std::vector<T>>& reads,
    std::vector<T>& unrecognizedPrimer, bool condensed) {
  size_t maxSize = 0;
  for (const auto& currentPrimer : primers) {
    if (currentPrimer.second.first.size() > maxSize) {
      maxSize = currentPrimer.second.first.size();
    }
  }
  for (auto& read : inputReads) {
    if (read.seqBase_.seq_.size() < maxSize) {
      unrecognizedPrimer.push_back(read);
    } else {

      bool foundMatch = false;
      for (const auto& currentPrimer : primers) {
        if (currentPrimer.first == "barcodeSize") {
          continue;
        }
        // find forward primer or if it isn't found to the unrecognized primer
        // vector you go
        std::string primerComparer;
        std::string readComparer;
        if (condensed) {
          primerComparer = condenseSeqSimple(currentPrimer.second.first);
          read.createCondensedSeq();
          readComparer = read.condensedSeq.substr(0, primerComparer.size());
        } else {
          primerComparer = currentPrimer.second.first;
          readComparer =
              read.seqBase_.seq_.substr(0, currentPrimer.second.first.size());
        }
        if (primerComparer == readComparer) {
          changeSubStrToLowerFromBegining(read.seqBase_.seq_,
                                          currentPrimer.second.first.length());
          if (reads.find(currentPrimer.first) == reads.end()) {
            reads.insert({currentPrimer.first, {read}});
          } else {
            reads[currentPrimer.first].push_back(read);
          }
          foundMatch = true;
          break;
        }
      }
      if (!foundMatch) {
        unrecognizedPrimer.push_back(read);
      }
    }
  }
}

template <typename T>
void findForwardPrimerAndSortDegenerative(
    std::vector<T>& inputReads,
    const std::map<std::string, std::pair<std::string, std::string>>& primers,
    std::map<std::string, std::vector<T>>& reads,
    std::vector<T>& unrecognizedPrimer,uint32_t variableStart,
    double queryCoverage,
    uint32_t numOfMismatches, double percentageGaps) {
  uint64_t maxSize = 0;
  for (const auto& currentPrimer : primers) {
    if (currentPrimer.second.first.size() > maxSize) {
      maxSize = currentPrimer.second.first.size();
    }
  }
  uint64_t maxReadLength = 0;
  readVec::getMaxLength(inputReads, maxReadLength);
  auto scoringMatrixMap = substituteMatrix::createDegenScoreMatrix(1,-1);

  // create alinger class object
  auto gapPars = gapScoringParameters (7, 1, 7, 1, 7, 1);
  aligner alignerObj =
      aligner(maxReadLength, gapPars, scoringMatrixMap);

  int count = 1;
  std::cout << "Finding forward primer with degenerative bases in mind"
            << std::endl;
  for (auto& read : inputReads) {
    if (count % 5000 == 0) {
      std::cout << "Curently on " << count << " of " << inputReads.size()
                << std::endl;
    }
    ++count;
    if (read.seqBase_.seq_.size() < maxSize) {
      unrecognizedPrimer.push_back(read);
    } else {
      bool foundMatch = false;
      for (const auto& currentPrimer : primers) {
        if (currentPrimer.first == "barcodeSize") {
          continue;
        }
        // find forward primer or if it isn't found to the unrecognized primer
        // vector you go
        auto primer =
            baseReadObject(seqInfo("primer", currentPrimer.second.first));
        auto readBegin =
                    baseReadObject(seqInfo("readBegin", read.seqBase_.seq_.substr(0, variableStart + maxSize)));
        /*auto forwardPosition = alignerObj.findReversePrimer(read.seqBase_.seq_, primer.seqBase_.seq_);
        alignerObj.rearrangeObjs(read.seqBase_, primer.seqBase_, true);
        alignerObj.profilePrimerAlignment(read, primer, false);*/
        auto forwardPosition = alignerObj.findReversePrimer(readBegin.seqBase_.seq_, primer.seqBase_.seq_);
				alignerObj.rearrangeObjs(readBegin, primer.seqBase_, true);
				alignerObj.profilePrimerAlignment(readBegin, primer, false);
        if (forwardPosition.first <= variableStart
            && alignerObj.comp_.distances_.queryCoverage_ >= queryCoverage
            && alignerObj.comp_.distances_.percentageGaps_ <= percentageGaps
            && alignerObj.comp_.hqMismatches_ <= numOfMismatches) {
          changeSubStrToLowerFromBegining(read.seqBase_.seq_,
                                          currentPrimer.second.first.length());
          reads[currentPrimer.first].emplace_back(read);
          foundMatch = true;
          break;
        }
      }
      if (!foundMatch) {
        unrecognizedPrimer.push_back(read);
      }
    }
  }
}

/**@b
 *
 * @param read
 * @param reversePrimer
 * @param runParmas
 * @param alignObj
 * @param reverserPrimerToLowerCase
 * @param weighHomopolyers
 * @param queryCoverageCutOff
 * @param percentIdentityCutoff
 *
 * @todo needs to be optimized
 */
template <typename T>
void checkForReversePrimer(T& read, readObject& reversePrimer,
                           runningParameters runParmas, aligner& alignObj,
                           bool reverserPrimerToLowerCase,
                           bool weighHomopolyers, double queryCoverageCutOff,
                           double percentIdentityCutoff) {
  auto rPos = alignObj.findReversePrimer(read, reversePrimer);
  alignObj.rearrangeObjs(read.seqBase_, reversePrimer.seqBase_, true);
  alignObj.profilePrimerAlignment(read, reversePrimer, weighHomopolyers);
  bool primerGood = true;
  if (alignObj.comp_.largeBaseIndel_ > runParmas.errors_.largeBaseIndel_ ||
      alignObj.comp_.oneBaseIndel_ > runParmas.errors_.oneBaseIndel_ ||
      alignObj.comp_.twoBaseIndel_ > runParmas.errors_.twoBaseIndel_) {
    primerGood = false;
  }
  if (alignObj.comp_.distances_.percentIdentity_ < percentIdentityCutoff) {
    primerGood = false;
  }
  if (alignObj.comp_.distances_.queryCoverage_ < queryCoverageCutOff) {
    primerGood = false;
  }

  read.remove = !primerGood;
  read.seqBase_.on_ = primerGood;
  if (primerGood) {
    if (reverserPrimerToLowerCase) {
      while (read.seqBase_.seq_[rPos.first] ==
             read.seqBase_.seq_[rPos.first - 1]) {
        --rPos.first;
      }
      changeSubStrToLowerToEnd(read.seqBase_.seq_, rPos.first);
    }
    read.setClip(0, rPos.second);
  }
}

template<typename T>
void checkForReversePrimer(std::vector<T>& reads, readObject& reversePrimer,
		runningParameters runParmas, aligner& alignObj,
		bool reverserPrimerToLowerCase, bool weighHomopolyers,
		double queryCoverageCutOff, double percentIdentityCutoff) {
	for (auto& read : reads) {
		checkForReversePrimer(read, reversePrimer, runParmas, alignObj,
				reverserPrimerToLowerCase, weighHomopolyers, queryCoverageCutOff,
				percentIdentityCutoff);
	}
}





template <typename T>
void sortReadsBySeqFront(const std::vector<T>& inputReads,
                         std::map<std::string, std::vector<readObject>>& reads,
                         size_t compareSize, std::vector<T>& unrecognizedPrimer,
                         size_t cutOffMinSize) {
  std::map<std::string, std::vector<readObject>> tempReads;

  for (const auto& read : inputReads) {
    if (read.seqBase_.seq_.size() < compareSize) {
      unrecognizedPrimer.push_back(read);
    } else {
      if (tempReads.find(read.seqBase_.seq_.substr(0, compareSize)) ==
          tempReads.end()) {
        tempReads.insert({read.seqBase_.seq_.substr(0, compareSize), {read}});
      } else {
        tempReads[read.seqBase_.seq_.substr(0, compareSize)].push_back(read);
      }
    }
  }
  for (const auto& read : tempReads) {
    if (read.second.size() > cutOffMinSize) {
      reads.insert(read);
    } else {
      addOtherVec(unrecognizedPrimer, read.second);
    }
  }
}

template <typename T>
void renameReadNames(std::vector<T>& reads, const std::string& stub,
                     bool processed, bool keepChimeraFlag = true,
                     bool keepCompFlag = true,
                     const std::string& sortBy = "none") {
	uint64_t count = 0;
  if (sortBy != "none") {
    readVecSorter::sortReadVector(reads, sortBy);
  }
  uint64_t maxSize = reads.size();
  for (auto& read : reads) {
    bool chimera = read.seqBase_.name_.find("CHI") != std::string::npos;
    bool comp = read.seqBase_.name_.find("_Comp") != std::string::npos;
    read.seqBase_.name_ = stub;
    if (chimera && keepChimeraFlag) {
      read.seqBase_.markAsChimeric();
    }
    if (comp && keepCompFlag) {
    	read.seqBase_.name_ += "_Comp";
    }
    read.seqBase_.name_ += "." + leftPadNumStr(count, maxSize);
    if (processed) {
      read.seqBase_.name_ += "_t" + to_string(read.seqBase_.cnt_);
    }
    ++count;
  }
}

template <typename T>
void renameReadNamesNewClusters(std::vector<T>& reads, const std::string& stub,
                                bool processed, bool keepChimeraFlag = true,
                                bool keepCompFlag = true,
                                const std::string& sortBy = "none") {
  renameReadNames(reads, stub, processed, keepChimeraFlag, keepCompFlag, sortBy);
  for (auto& read : reads) {
    read.reads_[0].seqBase_.name_ = read.seqBase_.name_;
  }
}

VecStr getRecurrenceStemNames(const VecStr& allNames);



void processRunCutoff(uint32_t& runCutOff, const std::string& runCutOffString, int counter);

uint32_t processRunCutoff(const std::string& runCutOffString, uint64_t counter);




template <typename T>
VecStr findLongestSharedSeqFromReads(const std::vector<T>& reads) {
  VecStr seqs;
  for (const auto& rIter : reads) {
    seqs.push_back(rIter.seqBase_.seq_);
  }
  return seqUtil::findLongestShardedMotif(seqs);
}

void makeMultipleSampleDirectory(const std::string& barcodeFilename,
                                 const std::string& mainDirectoryName);

void makeSampleDirectoriesWithSubDirectories(const std::string&,
                                             const std::string&);

void processKrecName(readObject& read, bool post);



/////////tools for finding additional location output
std::string makeIDNameComparable(const std::string& idName);

std::string findIdNameFromFileName(const std::string& filename);

std::string processFileNameForID(const std::string& fileName);

std::string findAdditonalOutLocation(const std::string& locationFile,
                                     const std::string& fileName);
VecStr getPossibleDNASubstringsForProtein(const std::string& seq,
                                          const std::string& protein,
                                          const std::string& seqType = "DNA");
VecStr findPossibleDNA(const std::string& seq, const std::string& protein,
                       const std::string& seqType = "DNA",
                       bool checkComplement = true);

VecStr getAllCycloProteinFragments(const std::string& protein);

std::multimap<int, std::string> getProteinFragmentSpectrum(
    const std::string& protein);

std::vector<int> getRealPossibleWeights(const std::vector<int>& spectrum);



VecStr organizeLexicallyKmers(const std::string& input, size_t colNum);
VecStr organizeLexicallyKmers(const VecStr& input, int colNum);
uint64_t smallestSizePossible(uint64_t weight);

template <typename READ, typename REF>
std::vector<baseReadObject> alignToSeq(const std::vector<READ>& reads,
                                       const REF& reference, aligner& alignObj,
                                       bool local, bool usingQuality) {
  std::vector<baseReadObject> output;
  output.emplace_back(baseReadObject(reference));
  for (const auto& read : reads) {
    alignObj.alignVec(reference, read, local);
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
    alignObj.alignVec(reference, read, local);
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
    alignObj.alignVec(reference, read, local);
    output.push_back(alignObj.alignObjectB_.seqBase_.seq_);
  }
  return output;
}

int64_t getPossibleNumberOfProteins(
    int64_t weight, std::unordered_map<int64_t, int64_t>& cache);
probabilityProfile randomMotifSearch(const VecStr& dnaStrings, int kLength,
                                     int numberOfRuns, bool gibs, int gibsRuns,
                                     randomGenerator& gen);
probabilityProfile randomlyFindBestProfile(const VecStr& dnaStrings,
                                           const std::vector<VecStr>& allKmers,
                                           int numberOfKmers,
                                           randomGenerator& gen);
probabilityProfile randomlyFindBestProfileGibs(
    const VecStr& dnaStrings, const std::vector<VecStr>& allKmers,
    int numberOfKmers, int runTimes, randomGenerator& gen);
/*
template <typename T>
struct letterCompositionSorter {
  letterCompositionSorter(const std::string& letter) : _letter(letter) {}
  std::string _letter;
  bool operator()(const T& first, const T& second) const {
    return first.counter.fractions.at(_letter) <
           second.counter.fractions.at(_letter);
  }
};
template <typename T>
void sortByLetterComposition(std::vector<T>& vec, const std::string& letter,
                             bool decending) {
  letterCompositionSorter<T> comparer(letter);
  if (decending) {
    std::sort(vec.begin(), vec.end(), comparer);
  } else {
    std::sort(vec.rbegin(), vec.rend(), comparer);
  }
}*/
template <typename T>
std::vector<std::vector<T>> getAllVecCycloFragments(const std::vector<T>& vec) {
  std::vector<std::vector<T>> ans;
  ans.emplace_back(std::vector<int>(0));
  for (auto i : iter::range<uint32_t>(1, vec.size())) {
    for (auto j : iter::range(vec.size())) {
      std::vector<T> currentFragment;
      if (j + i > vec.size()) {
        currentFragment =
            catenateVectors(getSubVector(vec, j, vec.size() - j),
                            getSubVector(vec, 0, i - (vec.size() - j)));
      } else {
        currentFragment = getSubVector(vec, j, i);
      }
      ans.push_back(currentFragment);
    }
  }
  ans.push_back(vec);
  std::sort(ans.begin(), ans.end());
  return ans;
}
template <typename T>
std::vector<std::vector<T>> getAllVecLinearFragments(
    const std::vector<T>& vec) {
  std::vector<std::vector<T>> ans;
  ans.emplace_back(std::vector<int>(0));
  for (auto i : iter::range<uint32_t>(1, vec.size())) {
    for (auto j : iter::range(vec.size())) {
      std::vector<T> currentFragment;
      if (j + i > vec.size()) {

      } else {
        currentFragment = getSubVector(vec, j, i);
        ans.push_back(currentFragment);
      }
    }
  }
  ans.push_back(vec);
  std::sort(ans.begin(), ans.end());
  return ans;
}
VecStr getAllLinearProteinFragments(const std::string& protein);
std::vector<std::vector<int>> getPossibleProteinsForSpectrum(
    const std::string& spectrum, bool useAllWeights = false);
bool trimCycle(std::vector<std::vector<int>>& nextCycle,
               std::vector<std::vector<int>>& matchesSpectrum,
               const std::map<int, int>& spectrumToWeights,
               const std::vector<int> specVec);
std::vector<std::vector<int>> growNextCycle(
    const std::vector<std::vector<int>>& previousCycle,
    const std::vector<int>& possibleWeights);
int scoreSpectrumAgreement(const std::map<int, int>& spectrum,
                           const std::map<int, int>& currentSpectrum);
std::multimap<int, std::vector<int>, std::greater<int>> growNextCycleScore(
    const std::multimap<int, std::vector<int>, std::greater<int>>&
        previousCycle,
    const std::vector<int>& possibleWeights,
    const std::map<int, int>& spectrumCounts, int parentMass, bool linear);
std::multimap<int, std::vector<int>, std::greater<int>> trimCycleScore(
    std::multimap<int, std::vector<int>, std::greater<int>>& nextCycle,
    std::multimap<int, std::vector<int>, std::greater<int>>& matchesSpectrum,
    int parentMass, int leaderBoardNumber, int& currentLeader);
std::multimap<int, std::vector<int>, std::greater<int>>
    getPossibleProteinsForSpectrum(const std::string& spectrum,
                                   int leaderBoardNumber, bool verbose,
                                   bool useAllWeights, bool convolution,
                                   int convolutionCutOff, bool linear);
std::map<int, int> getConvolutionWeights(std::vector<int> experimentalSpectrum,
                                         int multiplicityCutOff,
                                         int lowerBound = 57,
                                         int upperBound = 200);
std::vector<int> convolutionWeights(std::vector<int> experimentalSpectrum,
                                    int multiplicityCutOff, int lowerBound = 57,
                                    int upperBound = 200);
std::vector<int> topConvolutionWeights(std::vector<int> experimentalSpectrum,
                                       int mItems, int lowerBound = 57,
                                       int upperBound = 200);
int64_t getMinCoins(int64_t change, const std::vector<int64_t>& coins,
                    std::unordered_map<int64_t, int64_t>& cache);

template <typename T, typename CLUSTER>
std::vector<CLUSTER> roughCreateOTU(const std::vector<T>& reads,
                                    aligner& alignerObj, double percentCutOff,
                                    double gapCutOff, double queryCutoff, bool local,
                                    bool weighHomopolymers) {
  std::vector<CLUSTER> ans;
  uint32_t tenPercent = .1 * readVec::getTotalReadCount(reads);
  uint32_t count = 0;

  for (const auto& read : reads) {
    if (count % tenPercent == 0) {
      std::cout << "On " << count << " of "
                << readVec::getTotalReadCount(reads) << std::endl;
    }
    ++count;
    bool foundAnOtu = false;
    for (auto& compare : ans) {
      alignerObj.alignVec(compare, read, local);
      alignerObj.profilePrimerAlignment(compare, read, weighHomopolymers);
      /*
      std::cout << compare.seqBase_.name_ << std::endl;
      std::cout << alignerObj.alignObjectA_.seqBase_.seq_ << std::endl;
      std::cout << read.seqBase_.name_ << std::endl;
      std::cout << alignerObj.alignObjectB_.seqBase_.seq_ << std::endl;
      std::cout << "pi: " << alignerObj.percentIdentity_ <<"\tpg: "
      << alignerObj.percentageGaps_ << "\tqc: "<< alignerObj.queryCoverage_ <<
      std::endl;
      std::cout << "pi - pg " << alignerObj.percentIdentity_ -
      alignerObj.percentageGaps_;
      std::cout <<std::endl;*/
      // if (alignerObj.percentIdentity_ >=percentCutOff &&
      // alignerObj.percentageGaps_ <= gapCutOff) {
      if (alignerObj.comp_.distances_.percentIdentity_ >= percentCutOff &&
          gapCutOff >= alignerObj.comp_.distances_.percentageGaps_ &&
          alignerObj.comp_.distances_.queryCoverage_ >=queryCutoff) {
        // std::cout << alignerObj.percentIdentity_ << std::endl;
        // std::cout << alignerObj.alignObjectA_.seqBase_.seq_ << std::endl;
        // std::cout << alignerObj.alignObjectB_.seqBase_.seq_ << std::endl;
        // readObject tempObject = readObject(seqInfo(
        // read.getStubName(false), read.seqBase_.seq_, read.seqBase_.qual_));
        // tempObject.seqBase_.cnt_ = read.seqBase_.cnt_;
        // CLUSTER adding = CLUSTER(tempObject);
        CLUSTER adding = CLUSTER(readObject(read.seqBase_));
        compare.addRead(adding);
        foundAnOtu = true;
        break;
      }
    }
    if (!foundAnOtu) {
      readObject tempObject = readObject(seqInfo(
          read.getStubName(false), read.seqBase_.seq_, read.seqBase_.qual_));
      tempObject.seqBase_.cnt_ = read.seqBase_.cnt_;
      CLUSTER adding = CLUSTER(tempObject);
      ans.push_back(adding);
    }
  }
  return ans;
}
template <typename T>
std::vector<std::vector<T>> roughCreateOTULong(const std::vector<T>& reads,
                                               aligner& alignerObj,
                                               double percentCutOff, bool local,
                                               bool weighHomopolymers) {
  std::vector<std::vector<T>> otus;
  for (const auto& read : reads) {
    bool foundAnOtu = false;
    uint32_t otuPos = 0;
    for (auto& compare : otus) {
      for (auto& otu : compare) {
        alignerObj.alignVec(otu, read, local);
        alignerObj.profilePrimerAlignment(otu, read, weighHomopolymers);
        if (alignerObj.comp_.distances_.percentIdentity_ - alignerObj.comp_.distances_.percentageGaps_ >=
            percentCutOff) {
          foundAnOtu = true;
          break;
        }
      }
      if (foundAnOtu) {
        break;
      }
      ++otuPos;
    }
    if (!foundAnOtu) {
      otus.emplace_back(std::vector<T>{read});
    } else {
      otus[otuPos].emplace_back(read);
    }
  }
  return otus;
}
void processAlnInfoInput(aligner& alignerObj,
                         const std::string& alnInfoDirName);
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
std::unordered_map<std::string, bib::color> getColorsForNames(
    const std::vector<T>& reads, double sat, double lum) {
  std::unordered_map<std::string, bib::color> colorsForName;
  uint32_t count = 0;
  std::vector<bib::color> colors = bib::evenHuesAll(sat, lum, len(reads));
  for (const auto& read : reads) {
    colorsForName[read.getReadId()] = colors[count];
    ++count;
  }
  return colorsForName;
}




std::vector<uint32_t> getWindowQuals(const std::vector<uint32_t>& qual,
                                     uint32_t qWindowSize, uint32_t pos);
std::vector<double> likelihoodForBaseQ(
    const std::vector<uint32_t>& qual,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMeanQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMedianQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);
std::vector<double> likelihoodForMinQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods);

double getChangeInHydro(const char& firstAA, const char& secondAA);
template<typename T>
std::vector<double> getHydroChanges(const std::string& originalCodon,
                                    const VecStr& mutantCodons,
                                    const T& code) {
  std::vector<double> ans;
  if (originalCodon.size() != 3) {
    std::cout << "codon needs to be 3 bases" << std::endl;
    std::cout << originalCodon << std::endl;
    std::cout << originalCodon.size() << std::endl;
    return ans;
  }
  char originalAA = code.at(originalCodon);
  if (originalAA == '*') {
    return ans;
  }
  for (const auto& codon : mutantCodons) {
    char mutantAA = code.at(codon);
    if (mutantAA == '*') {
      continue;
    }
    ans.emplace_back(getChangeInHydro(originalAA, mutantAA));
  }
  return ans;
}
cluster createClusterFromSimiliarReads(std::vector<readObject> similiarReads,
                                       aligner& alignerObj);

template<typename READ>
std::vector<READ> vecStrToReadObjs(const VecStr & strs, const std::string & stubName){
	std::vector<READ> ans;
	for(const auto & strPos : iter::range(strs.size())){
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
	for(const auto & refPos : iter::range(len(refSeqs))){
		auto & ref = refSeqs[refPos];
		alignerObj.alignVec(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_, weightHomopolymers);
		if(alignerObj.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
			return foundAChimera;
		}
	}
	for(const auto & refPos : iter::range(len(refSeqs))){
		auto & ref = refSeqs[refPos];
		alignerObj.alignVec(ref.seqBase_, read.seqBase_, false);
		alignerObj.profilePrimerAlignment(ref.seqBase_, read.seqBase_, weightHomopolymers);
		if(alignerObj.mismatches_.empty()){
			foundAnExactMatch = true;
			foundAChimera = false;
		} else if(alignerObj.mismatches_.size() > 0 && !foundAnExactMatch){
			auto savedMismatches = alignerObj.mismatches_;
			//check front
			alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, 11, true, false, true, false,
					weightHomopolymers, 0,
					alignerObj.getAlignPosForSeqBPos(savedMismatches.begin()->second.seqBasePos));
			bool passFront = chiOverlap.passErrorProfile(alignerObj.comp_);

			//check back
			alignerObj.profileAlignment(ref.seqBase_, read.seqBase_, 11, true, false, true, false,
					weightHomopolymers,
								alignerObj.getAlignPosForSeqBPos(savedMismatches.rbegin()->second.seqBasePos + 1));
			bool passBack = chiOverlap.passErrorProfile(alignerObj.comp_);

			auto firstRefAlignA = alignerObj.alignObjectA_;
			auto firstRefAlignB = alignerObj.alignObjectB_;
			//std::cout << "pass front: " << passFront << std::endl;
			//std::cout << "pass back: " << passBack << std::endl;
			if(passFront){
				for(const auto & secondRefPos : iter::range(refPos + 1, len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}

					alignerObj.alignVec(secondRef.seqBase_, read.seqBase_, false);
					//check to see if from the mismatch on from the other one ref matches
					//to another ref
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers);
					if(alignerObj.mismatches_.size() < 1){
						continue;
					}
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers,
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
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, to_string(read.seqBase_.frac_),
							to_string(inflectionPoint),
							firstChiName, "1", to_string(inflectionPointPar1),
							secondChiName, "1", to_string(inflectionPointPar2)});
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
				for(const auto & secondRefPos : iter::range(refPos + 1,len(refSeqs))){
					auto & secondRef = refSeqs[secondRefPos];
					if(ref.seqBase_.name_ == secondRef.seqBase_.name_){
						continue;
					}
					alignerObj.alignVec(secondRef.seqBase_, read.seqBase_, false);
					//when comparing to ref's need to make sure it's not the sure sequence actually ahah
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers);
					if(alignerObj.mismatches_.size() < 1){
						continue;
					}
					alignerObj.profileAlignment(secondRef.seqBase_, read.seqBase_, 11, true, false, true, false, weightHomopolymers,
												0, alignerObj.getAlignPosForSeqBPos(savedMismatches.rbegin()->second.seqBasePos + 1));
					if(chiOverlap.passErrorProfile(alignerObj.comp_)){
						foundAChimera = true;
						firstChiName = ref.seqBase_.name_;
						secondChiName = secondRef.seqBase_.name_;
						inflectionPoint = savedMismatches.rbegin()->second.seqBasePos;
						inflectionPointPar1 = savedMismatches.rbegin()->second.refBasePos;
						std::string refTwoPortion = alignerObj.alignObjectA_.seqBase_.seq_.substr(0, alignerObj.getAlignPosForSeqAPos(savedMismatches.rbegin()->second.refBasePos) + 1);
						inflectionPointPar2 = removeCharReturn(refTwoPortion, '-').size() - 1;
						outInfo.content_.emplace_back(VecStr{read.seqBase_.name_, to_string(read.seqBase_.frac_),
							to_string(inflectionPoint),
							firstChiName, "1", to_string(inflectionPointPar1),
							secondChiName, "1", to_string(inflectionPointPar2)});
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

table getErrorFractions(const table & errorTab);
table getErrorFractionsCoded(const table & errorTab);
table getIndelSizeStats(const table & indelTab);
table getIndelDistribution(const table & indelTab);
table getErrorStats(const table & errorTab);
table getErrorDist(const table & errorTab);
table getSeqPosTab(const std::string & str);

uint32_t processCutOffStr(const std::string& runCutOffString,
  uint64_t readCount);






uint32_t getAlnPos(const std::string & seq, uint32_t realSeqPos);

uint32_t getRealPos(const std::string & seq, uint32_t seqAlnPos);


struct ExtractionInfo {
	ExtractionInfo(const table & allStatsTab, const table & allprofileTab,
			const table & profileBySampTab) :
			allStatsTab_(allStatsTab), allProfileTab_(allprofileTab), profileBySampTab_(
					profileBySampTab) {

	}
	ExtractionInfo(){}

	table allStatsTab_;
	table allProfileTab_;
	table profileBySampTab_;
};

ExtractionInfo collectExtractionInfo(const std::string & dirName, const std::string & indexToDir, const std::string & sampNames);
ExtractionInfo collectExtractionInfoDirectName(const std::string & dirName, const std::string & indexToDir, const std::string & sampNames);


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
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "seqToolsUtils.cpp"
#endif
