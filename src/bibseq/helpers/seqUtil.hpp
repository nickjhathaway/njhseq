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
//
//  seqUtil
//  ampliconCluster
//
//  Created by Nick Hathaway on 8/31/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//


#include "bibseq/utils.h"
#include "bibseq/objects/kmer/kmer.hpp"
#include "bibseq/seqToolsUtils/aminoAcidInfo.hpp"

#include "bibseq/objects/dataContainers/tables/table.hpp"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
#include "bibseq/objects/counters/charCounter.hpp"

namespace bibseq {

class seqUtil {

 public:
  // empty constructor, simply a class that holds some functions
  seqUtil() {}
  // reverse complement RNA or DNA strand
  static std::string reverseComplement(const std::string &seq,
                                       const std::string &seqType);
  static bool reversePalindrome(const std::string &seq,
                                const std::string &seqType);
  static std::map<size_t, size_t> findReversePalindromes(
      const std::string &seq, const std::string &seqType,
      size_t lowerSizeLimit = 4, size_t upperSizeLimit = 12);
  static void printOutReversePalindromes(const std::string &seq,
                                         const std::string &seqType,
                                         std::ostream &out,
                                         bool multipleAtSameSite,
                                         size_t lowerSizeLimit = 4,
                                         size_t upperSizeLimit = 12);
  // convert RNA or DNA strand to protein
  static std::string convertToProtein(const std::string &seq, size_t start = 0,
                                      bool forceStartM = false);
  static std::string convertOneCodon(const std::string &codon);
  // transcribe DNA
  static void transcribe(std::string &theWord);
  static void transcribeToRNA(std::string &theWord);
  static void convertToProteinFromcDNA(std::string &theWord, size_t start = 0,
                                       bool forceStartM = false);
  static std::string convertToProteinFromcDNAReturn(const std::string &theWord,
                                                    size_t start = 0,
                                                    bool forceStartM = false);
  static void changeToRNA(std::string &theWord);
  static std::string changeToRNAOut(const std::string &theWord);
  static int computeHammingDistance(const std::string &firstSeq,
                                    const std::string &secondSeq);
  // consensus
  // static std::string consensusFromVectorOfStrings(const VecStr &strings,
  // std::ostream &out);

  // translate phred score
  static char translatePhred(int qual);
  static int translatePhred(char qual);
  // read in primers
  static std::map<std::string, std::pair<std::string, std::string>> readPrimers(
      const std::string &idFileName, bool multiplex);
  static table readPrimers(const std::string &idFileName,
                           const std::string &fileDelim,
                           bool forceRead = false);
  static table readBarcodes(const std::string &idFileName,
                            const std::string &fileDelim, int &barcodeSize,
                            bool forceRead = false);
  // read in a scoring matrix
  static std::unordered_map<char, std::unordered_map<char, int>>
      readScoringMatrix(const std::string &fileName);
  static std::unordered_map<char, std::unordered_map<char, int>>
        readScoringMatrix(std::istream& in);

  static VecStr getAllProteins(const std::string &seq);
  static VecStr findAllOpenFrames(const VecStr &proteins);
  static VecStr findOpenFramesFromSeq(const std::string &seq);
  /**
   * @brief Convert a flow gram to sequence.
   *
   * @param flows The flows to be translated
   * @param flowSeq The flowOrder, defaults to TACG
   * @return The translated sequence
   */
  static std::string getSeqFromFlow(const std::vector<double> &flows,
                                    const std::string &flowSeq = "TACG");
  static VecStr getBarcodesInOrderTheyAppear(const std::string &fileName);

  // new
  static bool doesSequenceContainDegenerativeBase(const std::string &seq);
  static size_t checkQualityWindowPos(int windowSize, int minimumAverageQaul,
                                      int stepSize,
                                      const std::vector<uint32_t> &quality);
  static bool checkQualityWindow(int windowSize, int minimumAverageQaul,
                                 int stepSize,
                                 const std::vector<uint32_t> &quality);
  static void processQualityWindowString(const std::string &qualityWindowString,
                                         uint32_t &qualityWindowLength,
																				 uint32_t &qualityWindowStep,
																				 uint32_t &qualityWindowThres);
  static std::string findLongestSharedSubString(VecStr dnaStrings);
  static VecStr findLongestShardedMotif(VecStr dnaStrings);
  static double getAverageErrorRate(const std::vector<int> &qual);
  static std::map<std::string, int> makeDNAKmerCompMap(int kLength);
  static std::map<std::string, kmer> makeDNAKmerMap(int kLength);
  static double calculateWeightOfProteinDouble(const std::string &protein);
  static int calculateWeightOfProteinInt(const std::string &protein);
  static uint64_t getNumberOfPossibleDNAStrandsFromProtein(
      const std::string &protein);
  static std::string removeIntronsThenTranslate(std::string dna,
                                                const VecStr &introns);
  static double calculateTransitionTransversionRatio(const std::string &seq1,
                                                     const std::string &seq2);
  static bool isHomopolyer(const std::string &seq);
  static bool checkTwoEqualSeqs(const std::string &seq1,
                                const std::string &seq2,
                                int allowableMismatches);
  static std::map<std::string, kmer> adjustKmerCountsForMismatches(
      const std::map<kmer, int> &kmers, int allowableMismatches);
  static std::map<std::string, kmer> adjustKmerCountsForMismatches(
      const std::map<std::string, kmer> &kmers, int allowableMismatches);

  static std::unordered_map<std::string, uint32_t> getFuzzyKmerCount(const std::string &seq,
                                                      uint32_t kLength,
                                                      uint32_t allowableMutations,
                                                      bool checkComplement);

  static void removeLowerCase(std::string &sequence,
                              std::vector<uint32_t> &quality);

  static std::pair<std::string, std::vector<uint32_t>> removeLowerCaseReturn(
      std::string sequence, std::vector<uint32_t> quality);

  static int getCyclopeptideLengthFromSprectumLength(uint64_t length);
  static std::vector<std::vector<char>> getPossibleCyclopeptideFromSpretrum(
      const std::vector<int> &spectrum);
  static int getNumberOfPossibleLinearPeptides(uint64_t lengthOfProtein);
  static std::string removeGapsReturn(const std::string &seq);
  static void removeGaps(std::string &seq);
  static VecStr getFuzzySharedMotif(const VecStr &strs, uint32_t kLength,
  		uint32_t allowableMutations,
                                    bool checkComplement);
  static std::unordered_map<uint64_t, std::string> findMinimumHammingDistance(
      const std::string &seq, const std::string &subSeq, int kLength);
  static std::string createDegenerativeString(const VecStr &dnaString);
  static std::vector<uint32_t> rearrangeQuals(
      const std::vector<uint32_t> &qual,
      const std::vector<uint32_t> &positions);
  static std::vector<uint32_t> getQualPositions(const std::string &consensus,
                                                const std::string &compare);
  static uint32_t countMismatchesInAlignment(const std::string &ref,
                                             const std::string &compare,
                                             const char &ignore = '-');
  static void printQualCountsFiles(
      const std::string &workingDir, const std::string &seqName,
      std::map<std::string, std::map<double, uint32_t>> counts, bool overWrite,
      bool exitOnFailure);
  static void printMismatchQualCountsFiles(
      const std::string &workingDir, const std::string &seqName,
      std::map<std::string, std::map<double, uint32_t>> counts,
      std::map<std::string, std::map<double, uint32_t>> mismatchCounts,
      bool overWrite, bool exitOnFailure);
  static std::map<std::string,
                  std::unordered_map<std::string, std::vector<double>>>
      getCountsForModel(
          std::map<std::string, std::map<double, uint32_t>> counts,
          std::map<std::string, std::map<double, uint32_t>> mismatchCounts);

  static std::unordered_map<std::string, std::vector<double>>
      getCountsForSpecificModel(std::map<double, uint32_t> counts,
                                std::map<double, uint32_t> mismatchCounts);

  static std::map<std::string, std::unordered_map<double, double>>
      getTrueErrorRate(
          std::map<std::string, std::map<double, uint32_t>> counts,
          std::map<std::string, std::map<double, uint32_t>> mismatchCounts);
  static std::unordered_map<double, double> getTrueErrorRateSpecific(
      std::map<double, uint32_t> counts,
      std::map<double, uint32_t> mismatchCounts);

  static void rstripRead(std::string & str,
  		std::vector<uint32_t> & qual, char c);

  /*static void updateMismatchCounts(const std::string & consensus, const
     std::string & mutatant,
                                                                                                                                        const std::vector<uint32_t> & quals,  std::unordered_map<uint32_t, std::unordered_map<std::pair<char, char>, uint32_t>> & mutCounts);
  */
  /*
  static double probabilityOfKmer(const std::string & kmer,const
  std::unordered_map<uint, std::unordered_map<char, double>>& mapsOfProbs);
  static std::vector<kmer> mostProbableKmers(const std::string & seq, int
  kLength,const std::unordered_map<uint, std::unordered_map<char, double>>&
  mapsOfProbs);*/
}; // class seqUtils
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "seqUtil.cpp"
#endif
