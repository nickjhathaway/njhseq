#pragma once
//
//  seqUtil
//  ampliconCluster
//
//  Created by Nick Hathaway on 8/31/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//


#include "bibseq/utils.h"
#include "bibseq/objects/counters/letterCounter.hpp"
#include "bibseq/objects/helperObjects/kmer.hpp"
#include "kmerCalculator.hpp"
#include "bibseq/objects/dataContainers/table.hpp"
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

  // translate phred score
  static char translatePhred(int qual);
  static int translatePhred(char qual);
  // read in primers
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




  static bool doesSequenceContainDegenerativeBase(const std::string &seq);
  static size_t checkQualityWindowPos(int windowSize, int minimumAverageQaul,
                                      int stepSize,
                                      const std::vector<uint32_t> &quality);
  static bool checkQualityWindow(int windowSize, int minimumAverageQaul,
                                 int stepSize,
                                 const std::vector<uint32_t> &quality);
  static void processQualityWindowString(const std::string &qualityWindowString,
                                         int &qualityWindowLength,
                                         int &qualityWindowStep,
                                         int &qualityWindowThres);



  static bool isHomopolyer(const std::string &seq);

  static void removeLowerCase(std::string &sequence,
                              std::vector<uint32_t> &quality);

  static std::pair<std::string, std::vector<uint32_t>> removeLowerCaseReturn(
      std::string sequence, std::vector<uint32_t> quality);

  static std::string removeGapsReturn(const std::string &seq);
  static void removeGaps(std::string &seq);
  static std::vector<uint32_t> rearrangeQuals(
      const std::vector<uint32_t> &qual,
      const std::vector<uint32_t> &positions);
  static std::vector<uint32_t> getQualPositions(const std::string &consensus,
                                                const std::string &compare);
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


}; // class seqUtils
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "seqUtil.cpp"
#endif
