#pragma once
//
//  readObjectIO.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 11/19/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/sffObject.hpp"
#include "bibseq/seqToolsUtils.h"
#include "bibseq/IO/textFileReader.hpp"
#include "bibseq/readVectorManipulation.h"
#include <api/BamReader.h>

namespace bibseq {

struct readObjectIOOptions {
  readObjectIOOptions() {}
  readObjectIOOptions(const std::string& fName, const std::string& sName,
                      const std::string& tName, const std::string& iFormat,
                      const std::string& oFormat, const std::string& oFilename,
                      bool nameProcessed, const std::string& lCaseBases,
                      bool removingGaps, bool overWrite,
                      bool exitOnWriteFailure,bool forceWrite, bool includeSpaceInName)
      : firstName_(fName),
        secondName_(sName),
        thirdName_(tName),
        inFormat_(iFormat),
        outFormat_(oFormat),
        outFilename_(oFilename),
        processed_(nameProcessed),
        lowerCaseBases_(lCaseBases),
        removeGaps_(removingGaps),
        overWriteFile_(overWrite),
        exitOnFailureToWrite_(exitOnWriteFailure),
        forceWrite_(forceWrite),
        includeWhiteSpaceInName_(includeSpaceInName) {};
  std::string firstName_;
  std::string secondName_;
  std::string thirdName_;
  std::string inFormat_;
  std::string outFormat_;
  std::string outFilename_;
  bool processed_;
  std::string lowerCaseBases_;
  bool removeGaps_;
  bool overWriteFile_;
  bool exitOnFailureToWrite_;
  bool forceWrite_;
  bool includeWhiteSpaceInName_;
};



class readObjectIO {

 public:
  readObjectIO() {
  	lineBuffer_ = VecStr(10);
  	//includeSpaceInNames = true;
  }
  std::vector<readObject> reads;
  std::vector<sffObject> sffReads;
  std::vector<std::pair<size_t, size_t>> index;
  VecStr lineBuffer_;
  uint32_t bufferPos_ = UINT32_MAX;
  uint32_t lineNum_ = 0;
  uint32_t bufferMax_ = 10;
  bool refillBuffer(std::istream & is);
  bool includeSpaceInNames = true;
  void readFastaStream(std::istream& is, bool processed, bool add);

  bool readNextFastaStream(std::istream& is, readObject& read, bool processed);

  void readFastqStream(std::istream& is, uint32_t offSet,
                                     bool processed, bool add);
  bool readNextFastqStream(std::istream& is, uint32_t offSet, readObject& read,
                           bool processed);
  void readBam(std::string filename, bool processed);
  bool readNextBam(BamTools::BamReader & bReader, readObject& read,
  		BamTools::BamAlignment & aln, bool processed);
  // single Readers
  void read(const std::string& fileType, std::string filename,
            bool processed = false);
  // combo readers
  void read(const std::string& fileType, std::string filename,
            std::string secondName, bool processed = false);
  // three file readers, only one is fasta, qual, flow reader as of 8/6/13
  void read(const std::string& fileType, std::string filename,
            std::string secondName, std::string thirdName,
            bool processed = false);

  void read(const readObjectIOOptions& options);

  /*! \brief Internal Write
   *
   *
   *  Write the internal vector reads and write it to file
   */
  void write(std::string filename, std::string format, bool overWrite,
             bool exitOnFailure) {
    // seqWriter writer;
    if (format == "sff") {
      format = "fastq";
    }
    writeVector(reads, filename, format, overWrite, exitOnFailure);
  }
  /*! \brief Internal Write with optioins
   *
   *
   *  Write the internal vector reads and write it to file
   */
  void write(const readObjectIOOptions& options) {
    // seqWriter writer;
    if (options.outFormat_ == "sff") {
      writeVector(reads, options.outFilename_, "fastq", options.overWriteFile_,
                  options.exitOnFailureToWrite_);
    } else {
      writeVector(reads, options.outFilename_, options.outFormat_,
                  options.overWriteFile_, options.exitOnFailureToWrite_);
    }
  }
  void write(const std::string & filename, const readObjectIOOptions& options) {
    // seqWriter writer;
    if (options.outFormat_ == "sff") {
      writeVector(reads, filename, "fastq", options.overWriteFile_,
                  options.exitOnFailureToWrite_);
    } else {
      writeVector(reads, filename, options.outFormat_,
                  options.overWriteFile_, options.exitOnFailureToWrite_);
    }
  }
  /*! \brief External Write
   *
   *
   *  This function takes a vector of any class of read (readObject, cluster, or
   *sffObject) writes to file in given format
   */
  template <typename T>
  static void write(const std::vector<T>& inputReads, std::string filename,
                    std::string format, bool overWrite, bool exitOnFailure,
                    int extra = 0) {
    // seqWriter writer;
    if (format == "sff") {
      format = "fastq";
    }
    writeVector(inputReads, filename, format, overWrite, exitOnFailure, extra);
  }
  /*! \brief External Write With Options
   *
   *
   *  This function takes a vector of any class of read (readObject, cluster, or
   *sffObject) writes to file in given format
   */
  template <typename T>
  static void write(const std::vector<T>& inputReads,
                    readObjectIOOptions options, int extra = 0) {
    // seqWriter writer;
    if (options.outFormat_ == "sff") {
      options.outFormat_ = "fastq";
    }
    writeVector(inputReads, options.outFilename_, options.outFormat_,
                options.overWriteFile_, options.exitOnFailureToWrite_, extra);
  }
  template <typename T>
	static void write(const std::vector<T>& inputReads, const std::string & filename,
										readObjectIOOptions options, int extra = 0) {
		// seqWriter writer;
		if (options.outFormat_ == "sff") {
			options.outFormat_ = "fastq";
		}
		writeVector(inputReads, filename, options.outFormat_,
								options.overWriteFile_, options.exitOnFailureToWrite_, extra);
	}

  static std::vector<readObject> getReferenceSeq(const std::string& refFilename,
                                                 const std::string& refFormat,
                                                 bool refProcessed,
                                                 uint64_t& maxLength);

 public:
  template <class T>
  static void writeVector(const std::vector<T>& reads, std::string stubName,
                          const std::string& format, bool overWrite,
                          bool exitOnFailure, int extra = 0) {
    if (format == "fastq") {
      writeFastqFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fasta") {
      writeFastaFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "flow") {
      writeFlowFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "flowRaw") {
      writeFlowRawFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fastaQual") {
      writeFastaQualFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fastaFlow") {
      writeFastaFlowFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fastqFlow") {
      writeFastqFlowFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fastaQualFlow") {
      writeFastaQualFlowFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "pyroData") {
      writePyroDataFile(reads, stubName, overWrite, exitOnFailure, extra);
    } else if (format == "clipFastaQual") {
      writeClipFastaQualFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "condensedFastaQual") {
      writeCondensedFastaQualFile(reads, stubName, overWrite, exitOnFailure);
    } else {
      std::cout << "Unrecognized format " << format << " while writing "
                << stubName << std::endl;
    }
  }

  template <class T>
  static void writeSimpleVector(const std::vector<T>& reads,
                                std::string stubName, const std::string& format,
                                bool overWrite, bool exitOnFailure) {
    if (format == "fastq") {
      writeFastqFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fasta") {
      writeFastaFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fastaQual") {
      writeFastaQualFile(reads, stubName, overWrite, exitOnFailure);
    } else {
      std::cout << "Unrecognized format " << format << " while writing "
                << stubName << std::endl;
    }
  }
  template <class T>
  static void writeSimpleVector(const std::vector<T>& reads,
                                const readObjectIOOptions& options) {
    writeSimpleVector(reads, options.outFilename_, options.outFormat_,
                      options.overWriteFile_, options.exitOnFailureToWrite_);
    return;
  }
  static const uint32_t IlluminaQualOffset;
  static const uint32_t SangerQualOffset;

 private:
  // void openTextFile(std::ofstream & file,std::string filename, std::string
  // fileExtention, bool overWrite, bool exitOnFailure);

  template <class T>
  static void writeFastqFile(const std::vector<T>& reads,
                             std::string fastqFileName, bool overWrite,
                             bool exitOnFailure) {
    std::ofstream fastqFile;
    openTextFile(fastqFile, fastqFileName, ".fastq", overWrite, exitOnFailure);

    for (const auto& read : reads) {
      read.seqBase_.outPutFastq(fastqFile);
    }
    fastqFile.close();
  }

  template <class T>
  static void writeFastaFile(const std::vector<T>& reads,
                             std::string fastaFilename, bool overWrite,
                             bool exitOnFailure) {
    std::ofstream fastaFile;
    openTextFile(fastaFile, fastaFilename, ".fasta", overWrite, exitOnFailure);
    for (const auto& read : reads) {
      read.seqBase_.outPutSeq(fastaFile);
    }
    fastaFile.close();
  }
  template <class T>
  static void writeFlowFile(const std::vector<T>& reads,
                            std::string flowValueFilename, bool overWrite,
                            bool exitOnFailure) {
    std::ofstream flowValuesFile;
    openTextFile(flowValuesFile, flowValueFilename, ".data", overWrite,
                 exitOnFailure);
    for (const auto & read : reads) {
      read.outPutFlows(flowValuesFile);
    }
    flowValuesFile.close();
  }
  template <class T>
  static void writeFlowRawFile(const std::vector<T>& reads,
                               std::string rawFlowValuesFilename,
                               bool overWrite, bool exitOnFailure) {

    std::ofstream rawFlowValuesFile;
    openTextFile(rawFlowValuesFile, rawFlowValuesFilename, ".raw", overWrite,
                 exitOnFailure);
    for (const auto & read : reads) {
      read.outPutFlowsRaw(rawFlowValuesFile);
    }
    rawFlowValuesFile.close();
  }
  template <class T>
  static void writeQualFile(const std::vector<T>& reads,
                            std::string qualFilename, bool overWrite,
                            bool exitOnFailure) {
    std::ofstream fastaQualFile;
    openTextFile(fastaQualFile, qualFilename, ".fasta.qual", overWrite,
                 exitOnFailure);
    for (const auto& read : reads) {
      read.seqBase_.outPutQual(fastaQualFile);
    }
    fastaQualFile.close();
  }
  template <class T>
  static void writeFastaQualFile(const std::vector<T>& reads,
                                 std::string stubname, bool overWrite,
                                 bool exitOnFailure) {
    writeFastaFile(reads, stubname, overWrite, exitOnFailure);
    writeQualFile(reads, stubname, overWrite, exitOnFailure);
  }
  template <class T>
  static void writeFastaFlowFile(const std::vector<T>& reads,
                                 std::string stubName, bool overWrite,
                                 bool exitOnFailure) {
    writeFastaFile(reads, stubName, overWrite, exitOnFailure);
    writeFlowFile(reads, stubName, overWrite, exitOnFailure);
  }
  template <class T>
  static void writeFastqFlowFile(const std::vector<T>& reads,
                                 std::string stubName, bool overWrite,
                                 bool exitOnFailure) {
    writeFastqFile(reads, stubName, overWrite, exitOnFailure);
    writeFlowFile(reads, stubName, overWrite, exitOnFailure);
  }
  template <class T>
  static void writeFastaQualFlowFile(const std::vector<T>& reads,
                                     std::string stubName, bool overWrite,
                                     bool exitOnFailure) {
    writeFastaQualFile(reads, stubName, overWrite, exitOnFailure);
    writeFlowFile(reads, stubName, overWrite, exitOnFailure);
  }
  template <class T>
  static void writePyroDataFile(const std::vector<T>& reads,
                                std::string pryoDataFilename, bool overWrite,
                                bool exitOnFailure, int extra) {
    std::ofstream pyroDataFile;
    openTextFile(pyroDataFile, pryoDataFilename, ".dat", overWrite,
                 exitOnFailure);
    pyroDataFile << reads.size() << " " << extra << std::endl;
    for (const auto & read : reads) {
      read.outPutPyroData(pyroDataFile);
    }
    pyroDataFile.close();
  }
  template <class T>
  static void writeClipFastaFile(const std::vector<T>& reads,
                                 std::string clipFastaFilename, bool overWrite,
                                 bool exitOnFailure) {
    std::ofstream clipFastaFile;
    openTextFile(clipFastaFile, clipFastaFilename, ".clip.fasta", overWrite,
                 exitOnFailure);
    for (const auto & read : reads) {
      read.outPutSeqClip(clipFastaFile);
    }
    clipFastaFile.close();
  }
  template <class T>
  static void writeClipQualFile(const std::vector<T>& reads,
                                std::string clipFastaQualFilename,
                                bool overWrite, bool exitOnFailure) {
    std::ofstream clipFastaQualFile;
    openTextFile(clipFastaQualFile, clipFastaQualFilename, ".clip.fasta.qual",
                 overWrite, exitOnFailure);
    for (const auto & read : reads) {
      read.outPutQualClip(clipFastaQualFile);
    }
    clipFastaQualFile.close();
  }
  template <class T>
  static void writeClipFastaQualFile(const std::vector<T>& reads,
                                     std::string stubName, bool overWrite,
                                     bool exitOnFailure) {
    writeClipFastaFile(reads, stubName, overWrite, exitOnFailure);
    writeQualFile(reads, stubName, overWrite, exitOnFailure);
  }
  template <class T>
  static void writeCondensedFastaFile(const std::vector<T>& reads,
                                      std::string fastaFilename, bool overWrite,
                                      bool exitOnFailure) {
    std::ofstream fastaFile;
    openTextFile(fastaFile, fastaFilename, ".condensed.fasta", overWrite,
                 exitOnFailure);
    for (const auto & read : reads) {
      read.outPutCondensedSeq(fastaFile);
    }
    fastaFile.close();
  }
  template <class T>
  static void writeCondensedQualFile(const std::vector<T>& reads,
                                     std::string qualFilename, bool overWrite,
                                     bool exitOnFailure) {
    std::ofstream qualFile;
    openTextFile(qualFile, qualFilename, ".condensed.fasta.qual", overWrite,
                 exitOnFailure);
    for (const auto & read : reads) {
      read.outPutCondensedQual(qualFile);
    }
    qualFile.close();
  }
  template <class T>
  static void writeCondensedFastaQualFile(const std::vector<T>& reads,
                                          std::string stubName, bool overWrite,
                                          bool exitOnFailure) {
    writeCondensedFastaFile(reads, stubName, overWrite, exitOnFailure);
    writeCondensedQualFile(reads, stubName, overWrite, exitOnFailure);
  }

 private:
  // individual readObject readers
  void readNextSeq(FILE* fastaFile, char& c, std::string& name,
                   std::string& seq);
  void readNextQual(FILE* file, char& c, std::string& name, std::string& qual);
  void readNextData(FILE* dataFile, char& c, std::string& name,
                    std::vector<double>& flows);
  void readNextFastq(FILE* file, char& c, std::string& name, std::string& seq,
                     std::string& qual);
  sffObject readNextSff(std::string& currentLine, std::fstream& sffFile);
  // read a name line (e.g. >NAME)
  void readNextName(FILE* file, std::string& name, char& c);

  // single files in
  void readQual(std::string filename, bool processed, bool add);
  void readFlow(std::string filename, bool processed, bool add);
  void readRaw(std::string filename, bool processed, bool add);

  // multiple files
  // two files
  void readFastaQual(std::string filename, std::string qualName,
                     bool processed);
  void readFastqFlow(std::string fastqFilename, std::string flowFilename,
                     bool processed);
  void readFastaFlow(std::string fastaFilename, std::string flowFilename,
                     bool processed);
  void readFlowQual(std::string flowFilename, std::string qualFilename,
                    bool processed);
  // triple files
  void readFastaQualFlow(std::string fastaFilename, std::string qualFilename,
                         std::string flowFilename, bool processed);
  // special cases
  void readPyroData(std::string filename, bool processed);
  void readSff(std::string filename);


  void readClustal(std::string filename, bool processed);
  void readClustalNg(std::string filename, bool processed);
  void readShorahOld(const std::string& shorahFilename);
  void readShorah(const std::string& shorahFilename);

 public:
  //reading sff binary info
  void readSffbin(const std::string & filename);
  bool readNextSffBin(std::ifstream& in, sffObject& read, int32_t numFlowReads);
  void readHeader(std::ifstream& in, sffBinaryHeader& header);
};


}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "readObjectIO.cpp"
#endif
