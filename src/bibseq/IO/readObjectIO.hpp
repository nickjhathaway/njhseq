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
//  readObjectIO.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 11/19/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/sffObject.hpp"
#include "bibseq/readVectorManipulation.h"
#include <api/BamReader.h>
#include <bibcpp/jsonUtils.h>
#include "bibseq/IO/readObjectIOOptions.hpp"
#include "bibseq/IO/cachedReader.hpp"
namespace bibseq {






class readObjectIO {
public:
  readObjectIO() {}
  std::vector<readObject> reads;
  std::vector<sffObject> sffReads;
  std::vector<std::pair<size_t, size_t>> index;
  bool includeSpaceInNames = true;

  void readFastaStream(std::istream& is, bool processed, bool add);
  bool readNextFastaStream(cachedReader& cReader, readObject& read, bool processed);
  //bool readNextFastaStream(std::istream& is, readObject& read, bool processed);

  void readFastaQualStream(std::istream& fastaReader,std::istream& qualReader, bool processed, bool add);
  bool readNextFastaQualStream(cachedReader& fastaReader, cachedReader& qualReader, readObject & read,bool processed);
  //bool readNextFastaQualStream(std::istream& fastaReader, std::istream& qualReader, readObject & read,bool processed);
  bool readNextQualStream(cachedReader& cReader, std::vector<uint32_t>& quals, std::string & name);
  //bool readNextQualStream(std::istream& is, std::vector<uint32_t>& quals, std::string & name);

  bool readNextFastqStream(std::istream& is, uint32_t offSet, readObject& read,
                           bool processed);
  void readFastqStream(std::istream& is, uint32_t offSet,
                                     bool processed, bool add);

  void readFastqStreamGz(std::string filename, uint32_t offSet,
                                     bool processed, bool add);

  void readBam(std::string filename, bool processed);
  bool readNextBam(BamTools::BamReader & bReader, readObject& read,
  		BamTools::BamAlignment & aln, bool processed);
  // single Readers
  void read(const std::string& fileType, std::string filename,
            bool processed = false);
  // combo readers
  void read(const std::string& fileType, std::string filename,
            std::string secondName, bool processed = false);


  void read(const readObjectIOOptions& options);


  /*! \brief Internal Write with optioins
   *
   *
   *  Write the internal vector reads and write it to file
   */
  void write(readObjectIOOptions options) {
  	write(reads, options);
  }
  void write(const std::string & filename, readObjectIOOptions options) {
  	write(reads,filename, options);
  }

  /*! \brief External Write With Options
   *
   *
   *  This function takes a vector of any class of read (readObject, cluster, or
   *sffObject) writes to file in given format
   */
  template <typename T>
  static void write(const std::vector<T>& inputReads, readObjectIOOptions options) {
    // seqWriter writer;
    if (options.outFormat_ == "sff") {
      options.outFormat_ = "fastq";
    }
    writeVector(inputReads, options);
  }
  template <typename T>
	static void write(const std::vector<T>& inputReads, const std::string & filename,readObjectIOOptions options) {
		// seqWriter writer;
		if (options.outFormat_ == "sff") {
			options.outFormat_ = "fastq";
		}
		options.outFilename_ = filename;
		writeVector(inputReads, options);
	}

  static std::vector<readObject> getReferenceSeq(const std::string& refFilename,
                                                 const std::string& refFormat,
                                                 bool refProcessed,
                                                 uint64_t& maxLength);

 public:
	template<class T>
	static void writeVector(const std::vector<T>& reads,
			readObjectIOOptions options) {
		writeVector(reads, options.outFilename_, options.outFormat_,
				options.overWriteFile_, options.append_, options.exitOnFailureToWrite_,
				options.extra_);
	}
  template <class T>
  static void writeVector(const std::vector<T>& reads, std::string stubName,
                          const std::string& format, bool overWrite, bool append,
                          bool exitOnFailure, int extra) {
    if (format == "fastq") {
      writeFastqFile(reads, stubName, overWrite, append,exitOnFailure);
    } else if (format == "fasta") {
      writeFastaFile(reads, stubName, overWrite, append, exitOnFailure);
    } else if (format == "flow") {
      writeFlowFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "flowRaw") {
      writeFlowRawFile(reads, stubName, overWrite, exitOnFailure);
    } else if (format == "fastaQual") {
      writeFastaQualFile(reads, stubName, overWrite, append, exitOnFailure);
    } else if (format == "pyroData") {
      writePyroDataFile(reads, stubName, overWrite, exitOnFailure, extra);
    } else if (format == "mothurData") {
      writeMothurDataFile(reads, stubName, overWrite, exitOnFailure, extra);
    } else if (format == "condensedFastaQual") {
      writeCondensedFastaQualFile(reads, stubName, overWrite, exitOnFailure);
    } else {
    	std::stringstream ss;
      ss << bib::bashCT::bold << bib::bashCT::red << "Unrecognized format " << format << " while writing "
                << stubName << bib::bashCT::reset;
      throw std::runtime_error{ss.str()};
    }
  }

  template <class T>
  static void writeSimpleVector(const std::vector<T>& reads,
                                std::string stubName, const std::string& format,
                                bool overWrite, bool append, bool exitOnFailure) {
    if (format == "fastq") {
      writeFastqFile(reads, stubName, overWrite, append, exitOnFailure);
    } else if (format == "fasta") {
      writeFastaFile(reads, stubName, overWrite, append, exitOnFailure);
    } else if (format == "fastaQual") {
      writeFastaQualFile(reads, stubName, overWrite, append, exitOnFailure);
    } else {
    	std::stringstream ss;
      ss << bib::bashCT::bold << bib::bashCT::red << "Unrecognized format " << format << " while writing "
                << stubName << bib::bashCT::reset;
      throw std::runtime_error{ss.str()};
    }
  }
  template <class T>
  static void writeSimpleVector(const std::vector<T>& reads,
                                const readObjectIOOptions& options) {
    writeSimpleVector(reads, options.outFilename_, options.outFormat_,
                      options.overWriteFile_, options.append_, options.exitOnFailureToWrite_);
    return;
  }
  static const uint32_t IlluminaQualOffset;
  static const uint32_t SangerQualOffset;

 private:
  // void openTextFile(std::ofstream & file,std::string filename, std::string
  // fileExtention, bool overWrite, bool exitOnFailure);

  template <class T>
  static void writeFastqFile(const std::vector<T>& reads,
                             std::string fastqFileName, bool overWrite,bool append,
                             bool exitOnFailure) {
    std::ofstream fastqFile;
    bib::files::openTextFile(fastqFile, fastqFileName, ".fastq", overWrite, append, exitOnFailure);

    for (const auto& read : reads) {
      read.seqBase_.outPutFastq(fastqFile);
    }
    fastqFile.close();
  }

  template <class T>
  static void writeFastaFile(const std::vector<T>& reads,
                             std::string fastaFilename, bool overWrite,bool append,
                             bool exitOnFailure) {
    std::ofstream fastaFile;
    bib::files::openTextFile(fastaFile, fastaFilename, ".fasta", overWrite, append, exitOnFailure);
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
                            std::string qualFilename, bool overWrite,bool append,
                            bool exitOnFailure) {
    std::ofstream fastaQualFile;
    bib::files::openTextFile(fastaQualFile, qualFilename, ".fasta.qual", overWrite, append, exitOnFailure);
    for (const auto& read : reads) {
      read.seqBase_.outPutQual(fastaQualFile);
    }
    fastaQualFile.close();
  }
  template <class T>
  static void writeFastaQualFile(const std::vector<T>& reads,
                                 std::string stubname, bool overWrite,bool append,
                                 bool exitOnFailure) {
    writeFastaFile(reads, stubname, overWrite,append, exitOnFailure);
    writeQualFile(reads, stubname, overWrite, append, exitOnFailure);
  }


  template <class T>
  static void writePyroDataFile(const std::vector<T>& reads,
                                std::string pryoDataFilename, bool overWrite,
                                bool exitOnFailure, int extra) {
    std::ofstream pyroDataFile;
    openTextFile(pyroDataFile, pryoDataFilename, ".dat", overWrite,
                 exitOnFailure);
    pyroDataFile << reads.size() << " " << extra << "\n";
    for (const auto & read : reads) {
      read.outPutPyroData(pyroDataFile);
    }
    pyroDataFile.close();
  }
  template <class T>
  static void writeMothurDataFile(const std::vector<T>& reads,
                                std::string pryoDataFilename, bool overWrite,
                                bool exitOnFailure, int extra) {
    std::ofstream pyroDataFile;
    openTextFile(pyroDataFile, pryoDataFilename, ".flow", overWrite,
                 exitOnFailure);
    pyroDataFile << extra << "\n";
    for (const auto & read : reads) {
      read.outPutPyroData(pyroDataFile);
    }
    pyroDataFile.close();
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
  sffObject readNextSff(std::string& currentLine, std::fstream& sffFile);
  // multiple files
  void readFastaQual(std::string filename, std::string qualName,
                     bool processed);

  // special cases
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
