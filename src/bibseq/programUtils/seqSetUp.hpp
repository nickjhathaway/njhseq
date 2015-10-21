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
//  seqSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
#include "bibseq/IO.h"
#include "bibseq/alignment.h"
#include <bibcpp/progutils.h>
namespace bibseq {

class seqSetUp : public bib::progutils::programSetUp {
 public:
  seqSetUp(int argc, char* argv[]);
  seqSetUp(const bib::progutils::commandLineArguments& inputCommands);
  seqSetUp(const MapStrStr& inputCommands);
  void initializeDefaults();

  // seq read in names
  readObjectIOOptions ioOptions_;
  std::string seq_;
  readObject seqObj_;
  // directory name
  std::string directoryName_;
  bool overWriteDir_;

  //reporting
  bool verbose_;
  bool debug_;
  bool quiet_;

  // reference filename
  std::string refFilename_;
  std::string refSecondName_;
  std::string refFormat_;
  bool refProcessed_;

  // alignment for ref Info;
  std::string gapRef_;
  std::string gapLeftRef_;
  std::string gapRightRef_;
  gapScoringParameters gapInfoRef_;

  // alignmentInfo;
  std::string gap_;
  std::string gapLeft_;
  std::string gapRight_;
  gapScoringParameters gapInfo_;

  bool local_;
  // scoring matrix
  substituteMatrix scoring_;
  int32_t generalMatch_;
  int32_t generalMismatch_;

  // alignment profiling
  bool countEndGaps_;
  bool weightHomopolymers_;
  std::string qualThres_;
  bool eventBased_;
  uint32_t primaryQual_;
  uint32_t secondaryQual_;
  uint32_t qualThresWindow_;

  // kmer options
  uint32_t kLength_;
  // kmer profiling
  std::string runCutOffString_;
  uint32_t runCutoff_;
  bool kmersByPosition_;
  bool expandKmerPos_;
  uint32_t expandKmerSize_;

  //alnment caching
  std::string alnInfoDirName_;
  std::string outAlnInfoDirName_;
  bool writingOutAlnInfo_;

  //general clustering options
  bool skipOnLetterCounterDifference_;
  double fractionDifferenceCutOff_;
  bool adjustHomopolyerRuns_;
  bool largestFirst_;
  bool firstMatch_;
  uint32_t bestMatchCheck_;

  // private:


  bool processDefaultReader(bool readInNamesRequired = true);
  bool processReadInNames(bool required = true);
  void processGap();
  void processGapRef();
  void processQualThres();
  void processIteratorMap(std::string& parametersFile,
                          std::map<int, std::vector<double>>& iteratorMap);
  void processIteratorMapOnPerId(std::string& parametersFile,
                          std::map<int, std::vector<double>>& iteratorMap);
  void processKmerLenOptions();
  void processKmerProfilingOptions();
  void processScoringPars();
  void processAlignerDefualts();
  void processDirectoryOutputName(const std::string& defaultName,
                                  bool mustMakeDirectory);
  void processDirectoryOutputName(bool mustMakeDirectory);
  void processWritingOptions();
  bool processRefFilename(bool required = false);
  bool processSeq(bool required = false);
  bool processSeq(std::string& inputSeq, const std::string& flag,
                  const std::string& parName, bool required = false);
  bool processVerbose();
  bool processDebug();
  bool processQuiet();
  void processAlnInfoInput();

  void processClusteringOptions();
  void processSkipOnNucComp();
  void processAdjustHRuns();
  // usage prints
  void printInputUsage(std::ostream& out);
  void printAdditionaInputUsage(std::ostream& out,
                                const std::string& lowerRemove);

  void printGapUsage(std::ostream & out)const;

  void printKmerProfilingUsage(std::ostream& out);
  void printQualThresUsage(std::ostream& out);
  void printAlignmentUsage(std::ostream& out);
  void printReferenceComparisonUsage(std::ostream& out);
  void printFileWritingUsage(std::ostream& out, bool all);
  void printAlnInfoDirUsage(std::ostream& out);
  void printAdditionalClusteringUsage(std::ostream& out);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "seqSetUp.cpp"
#endif
