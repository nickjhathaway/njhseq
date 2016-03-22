#pragma once
//
//  seqSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
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
#include "bibseq/programUtils/SeqSetUpPars.hpp"

#include <bibcpp/progutils.h>
namespace bibseq {

class seqSetUp : public bib::progutils::programSetUp {
 public:
	using bib::progutils::programSetUp::programSetUp;

  SeqSetUpPars pars_;
  VecStr readInFormatsAvailable_ {"-sff", "-sffBin", "-fasta", "-fastq",
  	"-bam", "-fastqgz", "-fastagz", "-fastq1", "-fastq2"};
  void processQualityFiltering();
  bool processDefaultReader(bool readInNamesRequired = true);
  bool processDefaultReader(const VecStr & formats, bool readInNamesRequired = true);
  bool processReadInNames(const VecStr & formats, bool required = true);
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

  void processComparison(comparison & comp);

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
