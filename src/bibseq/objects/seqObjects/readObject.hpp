#pragma once
//
//  readObject.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 7/24/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
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
#include "bibseq/utils/utils.hpp"
#include "bibseq/objects/counters/charCounter.hpp"
#include "bibseq/objects/seqObjects/baseReadObject.hpp"
#include "bibseq/helpers/seqUtil.hpp"
#include <cppitertools/range.hpp>
namespace bibseq {

class readObject : public baseReadObject {

 public:
  // constructors
  // empty constructor
  readObject();
  readObject(const seqInfo& seqBase, bool processed = false);


  // constructor helpers
  void initializeOthers();
  void processRead(bool processed);





  std::string sampName;
  std::string expectsString;

  std::unordered_map<std::string, std::string> meta_;


  double averageErrorRate;

  //uint32_t basesAboveQualCheck_;
  double fractionAboveQualCheck_;
  bool remove;

  charCounter counter_;

  std::string condensedSeq;
  std::vector<uint32_t> condensedSeqQual;
  std::vector<std::pair<uint32_t, uint32_t>> condensedSeqQualPos;
  std::vector<int> condensedSeqCount;

	template<typename T>
	void addMeta(const std::string & key, const T & val, bool replace) {
		if (containsMeta(key) && !replace) {
			std::stringstream ss;
			ss << "Error in " << bib::bashCT::boldBlack(__PRETTY_FUNCTION__)
					<< " attempting to add meta that's already in meta_, use replace = true to replace"
					<< std::endl;
			throw std::runtime_error { ss.str() };
		} else {
			meta_[key] = estd::to_string(val);
		}
	}
  bool containsMeta(const std::string & key) const;
	std::string getMeta(const std::string & key) const;
	void processNameForMeta();
	//should only be called when meta data is already present in name
	void resetMetaInName();
	template<typename T>
	T getMeta(const std::string & key) const {
		return bib::lexical_cast<T>(getMeta(key));
	}

	virtual Json::Value toJson() const;


  void adjustHomopolyerRunQualities();



  void addQual(const std::string & stringQual);
  void addQual(const std::string & stringQual, uint32_t offSet);
  void addQual(const std::vector<uint32_t> & quals);

  void setFractionByCount(size_t totalNumberOfReads);

  virtual double getQualCheck(uint32_t qualCutOff)const ;
  virtual void setBaseCountOnQualCheck(uint32_t qualCheck);
  virtual void setLetterCount();
  virtual void setLetterCount(const std::vector<char> & alph);
  virtual double getGCContent();

  virtual void createCondensedSeq();

  // get the quality clipings used with the sff file
  void setClip(size_t leftPos, size_t rightPos);
  void setClip(size_t rightPos);
  void setClip(const std::pair<int, int>& positions);
  void trimFront(size_t upToPosNotIncluding);
  void trimBack(size_t fromPositionIncluding);
  // get various information about the qualities with the sequence
  double getAverageQual() const;
  virtual double getAverageErrorRate() const;
  uint32_t getSumQual() const;
  virtual void updateName();
  virtual void setName(const std::string& newName);
  void appendName(const std::string& add);;
  std::string getStubName(bool removeChiFlag) const;
  std::string getReadId() const;
  std::string getOtherReadSampName(const readObject& cr) const;
  std::string getOwnSampName() const;
  // find all occurences of the subsequence in the sequence
  std::vector<size_t> findSubsequenceOccurences(const std::string& findStr)
      const;
  void replace(const std::string& toBeReplaced, const std::string& replaceWith,
               bool allOccurences = true);



  // Protein conversion
  void translate(bool complement, bool reverse, size_t start = 0);
  readObject translateRet(bool complement, bool reverse, size_t start = 0) const;
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



  void outPutCondensedSeq(std::ostream& condensedSeqfile) const;
  void outPutCondensedQual(std::ostream& condensedQualFile) const;
  // check to seq to qual, output base and quality right next to each other
  void checkSeqQual(std::ostream& outFile) const;
  // clear the readObject, reset everything to an empty readObject
  void clearObject();
  // comparing reads
  bool operator>(const readObject& otherRead) const;
  bool operator<(const readObject& otherRead) const;
  //bool operator==(const readObject& otherRead) const;
  bool operator<=(const readObject& otherRead) const;
  bool operator>=(const readObject& otherRead) const;
  // description
  using size_type = baseReadObject::size_type;


  virtual ~readObject();
};

template<>
inline readObject::size_type len(const readObject & read){
	return read.seqBase_.seq_.size();
}
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "readObject.cpp"
#endif
