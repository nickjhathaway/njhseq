#pragma once
//
//  readObject.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 7/24/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//

#include "bibseq/utils/utils.hpp"
#include "bibseq/objects/counters/letterCounter.hpp"
#include "bibseq/objects/counters/charCounter.hpp"
#include "bibseq/objects/seqObjects/baseReadObject.hpp"
#include "bibseq/helpers/seqUtil.hpp"
#include <cppitertools/range.hpp>
namespace bibseq {

class readObject : public baseReadObject {

 public:
  // constructors
  // empty constructor
  readObject() : baseReadObject() { initializeOthers(); }
  readObject(const seqInfo& seqBase, bool processed = false);


  // constructor helpers
  void initializeOthers();
  void processRead(bool processed);

  std::vector<double> qualMeans_;
  std::vector<double> qualProbMeans_;
  void setQualMeans(uint32_t windowSize);
  void setQualProbMeans(uint32_t windowSize,
                        const std::array<double, 100>& errorLookUp);

  std::vector<double> flowValues;
  std::vector<double> processedFlowValues;

  std::string seqClip;
  std::vector<uint32_t> qualityClip;

  std::string sampName;
  std::string expectsString;
  int numberOfFlows;

  //double cumulativeFraction;
  //double normalizedFraction;
  double averageErrorRate;

  uint32_t basesAboveQualCheck_;
  double fractionAboveQualCheck_;
  bool remove;
  //letterCounter counter_;
  charCounterArray counter_;

  std::string condensedSeq;
  std::vector<uint32_t> condensedSeqQual;
  std::vector<std::pair<uint32_t, uint32_t>> condensedSeqQualPos;
  std::vector<int> condensedSeqCount;

  void adjustHomopolyerRunQualities();

  letterCounter condensedCounter;

  void addQual(const std::string & stringQual);
  void addQual(const std::string & stringQual, uint32_t offSet);
  void addQual(const std::vector<uint32_t> & quals);

  void setFractionByCount(size_t totalNumberOfReads);
  //void setNormalizedFraction(size_t totalNumberOfSamples);

  void setBaseCountOnQualCheck(uint32_t qualCheck);
  void setLetterCount();
  void setLetterCount(const std::vector<char> & alph);
  double getGCContent();

  // creates a condensed sequence, meaning all the homopolymer runs are
  // collapsed down into singles
  void createCondensedSeq();
  void setCondensedCounter();
  // get the string for the various vectors of info
  // std::string getQualString() const ;
  // std::string getFastqQualString(int offset) const ;
  // get the quality clipings used with the sff file
  void setClip(size_t leftPos, size_t rightPos);
  void setClip(size_t rightPos);
  void setClip(const std::pair<int, int>& positions);
  void trimFront(size_t upToPosNotIncluding);
  void trimBack(size_t fromPositionIncluding);
  // get various information about the qualities with the sequence
  double getAverageQual() const;
  double getAverageErrorRate() const;
  uint32_t getSumQual() const;
  // change the name of the object
  virtual void updateName();
  virtual void setName(const std::string& newName);
  void appendName(const std::string& add);
  void setFractionName();
  //void setCumulativeFractionName();
  //void setNormalizedFractionName();
  //void setNormalizedFractionByCumulativeFraction(size_t totalNumberOfSamples);
  std::string getStubName(bool removeChiFlag) const;
  std::string getReadId() const;
  std::string getOtherReadSampName(const readObject& cr) const;
  std::string getOwnSampName() const;
  // find all occurences of the subsequence in the sequence
  std::vector<size_t> findSubsequenceOccurences(const std::string& findStr)
      const;
  void replace(const std::string& toBeReplaced, const std::string& replaceWith,
               bool allOccurences = true);

  uint32_t getNumberOfBasesFromFlow(const std::vector<double>& flows);
  int clipToNinetyPercentOfFlows(size_t cutOff);
  size_t findFirstNoisyFlow(const std::vector<double>& flows);
  bool flowNoiseProcess(size_t cutoff);

  // Protein conversion
  void convertToProteinFromcDNA(bool transcribeToRNAFirst, size_t start = 0,
                                bool forceStartM = false);
  std::string getProteinFromcDNA(bool transcribeToRNAFirst, size_t start = 0,
                                 bool forceStartM = false) const;
  //
  void updateQualCounts(std::map<uint32_t, uint32_t>& qualCounts) const;
  void updateQualCounts(
      std::map<std::string, std::map<double, uint32_t>>& counts,
      int qualWindowSize, const std::array<double, 100>& qualErrorLookUp) const;
  void updateBaseQualCounts(std::map<double, uint32_t>& baseCounts,
                            uint32_t pos) const;
  void updateBaseQualCounts(std::map<double, uint32_t>& baseCounts) const;

  void updateQualWindowCounts(
      uint32_t pos, std::map<std::string, std::map<double, uint32_t>>& counts,
      int qualWindowSize) const;

  void updateQualWindowCounts(
      std::map<std::string, std::map<double, uint32_t>>& counts,
      int qualWindowSize) const;

  void updateQaulCountsAtPos(
      uint32_t pos, std::map<std::string, std::map<double, uint32_t>>& counts,
      int qualWindowSize, const std::array<double, 100>& qualErrorLookUp) const;
  // various outputs
  // void outPutFastq(std::ostream &fastqFile)const;
  // void outPutSeq(std::ostream &fastaFile)const;
  // void outPutQual(std::ostream &qualFile)const;

  void outPutFlows(std::ostream& flowsFile) const;
  void outPutFlowsRaw(std::ostream& flowdataFile) const;
  void outPutPyroData(std::ostream& pyroNoiseFile) const;
  void outPutSeqClip(std::ostream& seqClipfile) const;
  void outPutQualClip(std::ostream& qualClipfile) const;
  void outPutCondensedSeq(std::ostream& condensedSeqfile) const;
  void outPutCondensedQual(std::ostream& condensedQualFile) const;
  // check to seq to qual, output base and quality right next to each other
  void checkSeqQual(std::ostream& outFile) const;
  // clear the readObject, reset everything to an empty readObject
  void clearObject();
  // comparing reads
  bool operator>(const readObject& otherRead) const;
  bool operator<(const readObject& otherRead) const;
  bool operator==(const readObject& otherRead) const;
  bool operator<=(const readObject& otherRead) const;
  bool operator>=(const readObject& otherRead) const;
  // description
  virtual void printDescription(std::ostream& out, bool deep = false) const;
  using size_type = baseReadObject::size_type;
};

template<>
inline readObject::size_type len(const readObject & read){
	return read.seqBase_.seq_.size();
}
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "readObject.cpp"
#endif
