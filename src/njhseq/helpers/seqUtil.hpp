#pragma once
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
//  seqUtil
//
//  Created by Nick Hathaway on 8/31/12.
//


#include "njhseq/utils.h"
#include "njhseq/objects/kmer/kmer.hpp"
#include "njhseq/seqToolsUtils/aminoAcidInfo.hpp"

#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "njhseq/objects/counters/charCounter.hpp"

namespace njhseq {

class seqUtil {

 public:
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
                                      const std::vector<uint8_t> &quality);
  static bool checkQualityWindow(int windowSize, int minimumAverageQaul,
                                 int stepSize,
                                 const std::vector<uint8_t> &quality);
  static void processQualityWindowString(const std::string &qualityWindowString,
                                         uint32_t &qualityWindowLength,
																				 uint32_t &qualityWindowStep,
																				 uint8_t &qualityWindowThres);
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


  static void removeLowerCase(std::string &sequence,
                              std::vector<uint8_t> &quality);

  static std::pair<std::string, std::vector<uint8_t>> removeLowerCaseReturn(
      std::string sequence, std::vector<uint8_t> quality);

  static std::string removeGapsReturn(const std::string &seq);
  static void removeGaps(std::string &seq);
  static std::unordered_map<uint64_t, std::string> findMinimumHammingDistance(
      const std::string &seq, const std::string &subSeq, int kLength);
  static std::string createDegenerativeString(const VecStr &dnaString);
  static std::string genMotifStrAccountDegenBase(const std::string & primer);



  static std::vector<uint32_t> rearrangeQuals(
      const std::vector<uint8_t> &qual,
      const std::vector<uint32_t> &positions);
  static std::vector<uint32_t> getQualPositions(const std::string &consensus,
                                                const std::string &compare);
  static uint32_t countMismatchesInAlignment(const std::string &ref,
                                             const std::string &compare,
                                             const char &ignore = '-');


  static void rstripRead(std::string & str,
  		std::vector<uint8_t> & qual, char c);


}; // class seqUtils
}  // namespace njhseq


