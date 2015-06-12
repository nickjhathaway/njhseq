#pragma once
//
//  sffObject.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 4/26/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/readObject.hpp"

namespace bibseq {

class sffObject : public readObject {

 public:
  //sffObject() : readObject() { sffAdditionalInitialation(); }
  sffObject(const seqInfo &seqBase, bool processed = false)
      : readObject(seqBase) {
    processRead(processed);
    initializeOthers();
  }
  sffObject() : headerLength_(0), nameLength_(0), numBases_(0),
			clipQualLeft_(0), clipQualRight_(0), clipAdapterLeft_(0), clipAdapterRight_(0) { }
	~sffObject() { }

  void sffAdditionalInitialation();
  void addSffInfo(std::unordered_map<std::string, std::string> &info);

  /*
  std::vector<uint32_t> flowIndexes;
  std::vector<double> selectFlowValues;
  int numberOfBases;
  int clipQualLeft;
  int clipQualRight;
  int clipAdapLeft;
  int clipAdapRight;*/

  std::string getFlowValuesString() const;
  std::string getFlowIndexesString() const;

  void outPutSelectFlowvalues(std::ostream &selectFlowdataFile) const;
  void outPutFlowIndexes(std::ostream &flowIndexesFile) const;
  void outPutFlowValuesPryo(std::ostream &flowValuesPryoFile) const;
  void setSffClip(size_t left, size_t right);
  void setSelectFlowValues();

  //

	//members
	uint16_t headerLength_;
	uint16_t nameLength_;
	uint32_t numBases_;
	uint16_t clipQualLeft_;
	uint16_t clipQualRight_;
	uint16_t clipAdapterLeft_;
	uint16_t clipAdapterRight_;
	//std::string name_; //length depends on nameLength
	std::string timestamp_;
	std::string region_;
	std::string xy_;
	std::vector<uint32_t> flowIndex_;
	std::vector<double> selectFlowValues;
	//functions
	virtual void printDescription(std::ostream & out, bool deep = false) const;
	void printHeaderSffTxtStyle(std::ofstream& out);
	void printSffTxtSeqData(std::ostream& out) ;
	void decodeName();
	bool sanityCheck();
};
uint32_t fromBase36(std::string base36);
//main header info for sffBinary file
struct sffBinaryHeader {
	uint32_t magicNumber;
	std::string version;
	uint64_t indexOffset;
	uint32_t indexLength;
	uint32_t numReads;
	uint16_t headerLength;
	uint16_t keyLength;
	uint16_t numFlowsPerRead;
	int32_t flogramFormatCode;
	std::string flowChars; //length depends on number flow reads
	std::string keySequence; //length depends on key length

	sffBinaryHeader(){ magicNumber=0; indexOffset=0; indexLength=0; numReads=0; headerLength=0; keyLength=0; numFlowsPerRead=0; flogramFormatCode='s'; }
	~sffBinaryHeader() { }

	int printSffTxtStyle(std::ostream & out);
	void printDescription(std::ostream & out, bool deep = false) const;
};

readObject convertSffObject(const sffObject& theOtherObject) ;

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "sffObject.cpp"
#endif
