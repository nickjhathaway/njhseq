/*
 * SeqSetUpPars.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */
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
#include "SeqSetUpPars.hpp"

namespace bibseq {

SeqSetUpPars::SeqSetUpPars(){
	initializeDefaults();
}

void SeqSetUpPars::initializeDefaults() {
  ioOptions_.firstName_ = "";
  ioOptions_.secondName_ = "";
  ioOptions_.inFormat_ = SeqIOOptions::inFormats::NOFORMAT;
  ioOptions_.outFormat_ = SeqIOOptions::outFormats::NOFORMAT;
  ioOptions_.out_.outFilename_ = "out";
  ioOptions_.processed_ = false;
  ioOptions_.out_.append_ = false;
  //
  ioOptions_.removeGaps_ = false;
  ioOptions_.lowerCaseBases_ = "nothing";
  //
  ioOptions_.out_.overWriteFile_ = false;
  ioOptions_.out_.exitOnFailureToWrite_ = true;
  //
  ioOptions_.includeWhiteSpaceInName_ = true;
  //
  seq_ = "";
  seqObj_ = readObject(seqInfo("", seq_));

  //
  directoryName_ = "";
  overWriteDir_ = false;
  //


  //
  verbose_ = false;
  debug_ = false;
  quiet_ = false;

  //
  refIoOptions_.firstName_ = "";
  refIoOptions_.secondName_ = "";
  refIoOptions_.inFormat_ = SeqIOOptions::inFormats::NOFORMAT;
  refIoOptions_.processed_ = false;

  //
  gapRef_ = "5,1";
  gapInfoRef_.processGapStr(gapRef_, gapInfoRef_.gapOpen_, gapInfoRef_.gapExtend_);

  gapLeftRef_ = "0,0";
  gapInfoRef_.processGapStr(gapLeftRef_, gapInfoRef_.gapLeftOpen_, gapInfoRef_.gapLeftExtend_);

  gapRightRef_ = "0,0";
  gapInfoRef_.processGapStr(gapRightRef_, gapInfoRef_.gapRightOpen_, gapInfoRef_.gapRightExtend_);
  gapInfoRef_.setIdentifer();
  //
  gap_ = "7,1";
  gapInfo_.processGapStr(gap_, gapInfo_.gapOpen_, gapInfo_.gapExtend_);

  gapLeft_ = "7,1";
  gapInfo_.processGapStr(gapLeft_, gapInfo_.gapLeftOpen_, gapInfo_.gapLeftExtend_);

  gapRight_ = "0,0";
  gapInfo_.processGapStr(gapRight_, gapInfo_.gapRightOpen_, gapInfo_.gapRightExtend_);
  gapInfo_.setIdentifer();
  //
  local_ = false;
  generalMatch_ = 2;
  generalMismatch_ = -2;
  scoring_ = substituteMatrix::createDegenScoreMatrix(1,-1);
  //
  countEndGaps_ = false;
  weightHomopolymers_ = true;

  qualThres_ = "20,15";
  qScorePars_.primaryQual_ = 20;
  qScorePars_.secondaryQual_ = 15;
  qScorePars_.qualThresWindow_ = 2;

  eventBased_ = false;
  //
  alnInfoDirName_ = "";
  outAlnInfoDirName_ = "";
  writingOutAlnInfo_ = false;

  //
  runCutOffString_ = ".2%";
  runCutoff_ = 1;

  kLength_ = 9;
  kmersByPosition_ = true;
  expandKmerPos_ = false;
  expandKmerSize_ = 5;

  //general clustering
  skipOnLetterCounterDifference_ = false;
  fractionDifferenceCutOff_ = 0.05;
  adjustHomopolyerRuns_ = false;
  largestFirst_ = false;
  firstMatch_ = false;
  bestMatchCheck_ = 10;

  //quality filtering
  qFilPars_.checkingQWindow = false;
  qFilPars_.qualWindow_ = "50,5,25";
  qFilPars_.qualityWindowLength_ = 50;
  qFilPars_.qualityWindowStep_ = 5;
  qFilPars_.qualityWindowThres_ = 25;

  qFilPars_.checkingQFrac_ = false;
  qFilPars_.qualCheck_ = 30;
  qFilPars_.qualCheckCutOff_ = 0.75;

}

}  // namespace bibseq
