#pragma once
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
  seqSetUp(int argc, char* argv[]) : bib::progutils::programSetUp(argc, argv) {
    initializeDefaults();
  }
  seqSetUp(const bib::progutils::commandLineArguments& inputCommands)
      : bib::progutils::programSetUp(inputCommands) {
    initializeDefaults();
  }
  seqSetUp(const MapStrStr& inputCommands) : bib::progutils::programSetUp(inputCommands) {
    initializeDefaults();
  }

  void initializeDefaults() {
    ioOptions_.firstName_ = "";
    ioOptions_.secondName_ = "";
    ioOptions_.thirdName_ = "";
    ioOptions_.inFormat_ = "";
    ioOptions_.outFormat_ = "";
    ioOptions_.outFilename_ = "";
    ioOptions_.processed_ = false;
    ioOptions_.forceWrite_ = false;
    //
    seq_ = "";
    seqObj_ = readObject(seqInfo( "", seq_	));
    //
    ioOptions_.removeGaps_ = false;
    ioOptions_.lowerCaseBases_ = "nothing";
    //
    ioOptions_.overWriteFile_ = false;
    ioOptions_.exitOnFailureToWrite_ = false;
    //
    ioOptions_.includeWhiteSpaceInName_ = true;
    //
    directoryName_ = "";
    //
    refFilename_ = "";
    refSecondName_ = "";
    refFormat_ = "";
    refProcessed_ = false;
    //
    gap_ = "7.0,1";
    gapInfo_.processGapStr(gap_, gapInfo_.gapOpen_, gapInfo_.gapExtend_);

    gapLeft_ = "7.0,1";
    gapInfo_.processGapStr(gapLeft_, gapInfo_.gapLeftOpen_, gapInfo_.gapLeftExtend_);

    gapRight_ = "0.0,0.0";
    gapInfo_.processGapStr(gapRight_, gapInfo_.gapRightOpen_, gapInfo_.gapRightExtend_);
    gapInfo_.setIdentifer();
    verbose_ = false;
    debug_ = false;
    local_ = false;
    countEndGaps_ = false;
    weightHomopolymers_ = true;
    //

    qualThres_ = "20,15";
    primaryQual_ = 20;
    secondaryQual_ = 15;

    qualThesLowKmer_ = "30,25";
    primaryQualLowKmer_ = 30;
    secondaryQualLowKmer_ = 25;

    qualThresWindow_ = 5;

    eventBased_ = false;
    //
    alnInfoDirName_ = "";
    outAlnInfoDirName_ = "";
    writingOutAlnInfo_ = false;

    //
    runCutOffString_ = "1";
    runCutoff_ = 1;
    qualRunCutOffString_ = "1";
    qualRunCutoff_ = 1;
    kLength_ = 25;
    kmersByPosition_ = true;
    checkKmers_ = true;
    expandKmerPos_ = false;
    expandKmerSize_ = 5;

    //
    generalMatch_ = 2;
    generalMismatch_ = -2;
    scoring_ = substituteMatrix::createDegenScoreMatrix(1,-1);
    //general clustering
    regKmers_ = false;
    skipOnLetterCounterDifference_ = false;
    fractionDifferenceCutOff_ = 0.05;
    condensedCollapse_ = false;
    adjustHomopolyerRuns_ = false;

    quiet_ = false;
  }

  // seq read in names
  readObjectIOOptions ioOptions_;
  std::string seq_;
  readObject seqObj_;
  // main directoryName
  std::string directoryName_;
  // reference filename
  std::string refFilename_;
  std::string refSecondName_;
  std::string refFormat_;
  bool refProcessed_;
  // alignmentInfo;
  std::string gap_;
  std::string gapLeft_;
  std::string gapRight_;
  gapScoringParameters gapInfo_;

  bool verbose_;
  bool debug_;
  bool quiet_;
  bool local_;
  bool countEndGaps_;
  bool weightHomopolymers_;

  std::string qualThres_;
  bool eventBased_;
  int primaryQual_;
  int secondaryQual_;

  std::string qualThesLowKmer_;
  int primaryQualLowKmer_;
  int secondaryQualLowKmer_;

  int qualThresWindow_;

  //
  std::string alnInfoDirName_;
  std::string outAlnInfoDirName_;
  bool writingOutAlnInfo_;

  // kmer options
  std::string runCutOffString_;
  int runCutoff_;

  std::string qualRunCutOffString_;
  int qualRunCutoff_;

  int kLength_;
  bool kmersByPosition_;
  bool checkKmers_;
  bool expandKmerPos_;
  uint32_t expandKmerSize_;

  // scoring matrix
  //std::unordered_map<char, std::unordered_map<char, int>> scoringMatrixMap_;
  substituteMatrix scoring_;
  int generalMatch_;
  int generalMismatch_;
  //general clusreing
  bool regKmers_;
  bool skipOnLetterCounterDifference_;
  double fractionDifferenceCutOff_;
  bool condensedCollapse_;
  bool adjustHomopolyerRuns_;

  // private:
  void processRegKmers();
  void processSkipOnNucComp();
  void processCondensedCollapse();
  void processAdjustHRuns();

  bool processDefaultReader(bool readInNamesRequired = true);
  bool processReadInNames(bool required = true);
  void processGap();
  void processQualThres();
  void processIteratorMap(std::string& parametersFile,
                          std::map<int, std::vector<double>>& iteratorMap);
  void processKmerOptions();
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
  // usage prints
  void printInputUsage(std::ostream& out);
  void printAdditionaInputUsage(std::ostream& out,
                                const std::string& lowerRemove);
  void printKmerUsage(std::ostream& out);
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
