/*
 * SeqSetUpPars.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "SeqSetUpPars.hpp"

namespace njhseq {




Json::Value QualFilteringPars::toJson() const{
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["checkingQWindow"] = njh::json::toJson(checkingQWindow);

	ret["qualWindow_"] = njh::json::toJson(qualWindow_);
	ret["qualityWindowLength_"] = njh::json::toJson(qualityWindowLength_);
	ret["qualityWindowStep_"] = njh::json::toJson(qualityWindowStep_);
	ret["qualityWindowThres_"] = njh::json::toJson(qualityWindowThres_);
	ret["checkingQFrac_"] = njh::json::toJson(checkingQFrac_);
	ret["qualCheck_"] = njh::json::toJson(qualCheck_);
	ret["qualCheckCutOff_"] = njh::json::toJson(qualCheckCutOff_);
	ret["trimAtQualCutOff_"] = njh::json::toJson(trimAtQualCutOff_);
	ret["trimAtQual_"] = njh::json::toJson(trimAtQual_);

	return ret;
}

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
  gapInfoRef_.processGapStr(gapLeftRef_, gapInfoRef_.gapLeftQueryOpen_, gapInfoRef_.gapLeftQueryExtend_);
  gapInfoRef_.processGapStr(gapLeftRef_, gapInfoRef_.gapLeftRefOpen_, gapInfoRef_.gapLeftRefExtend_);

  gapRightRef_ = "0,0";
  gapInfoRef_.processGapStr(gapRightRef_, gapInfoRef_.gapRightQueryOpen_, gapInfoRef_.gapRightQueryExtend_);
  gapInfoRef_.processGapStr(gapRightRef_, gapInfoRef_.gapRightRefOpen_, gapInfoRef_.gapRightRefExtend_);
  gapInfoRef_.setIdentifer();
  //
  gap_ = "5,1";
  gapInfo_.processGapStr(gap_, gapInfo_.gapOpen_, gapInfo_.gapExtend_);

  gapLeft_ = "5,1";
  gapInfo_.processGapStr(gapLeft_, gapInfo_.gapLeftQueryOpen_, gapInfo_.gapLeftQueryExtend_);
  gapInfo_.processGapStr(gapLeft_, gapInfo_.gapLeftRefOpen_, gapInfo_.gapLeftRefExtend_);

  gapRight_ = "0,0";
  gapInfo_.processGapStr(gapRight_, gapInfo_.gapRightQueryOpen_, gapInfo_.gapRightQueryExtend_);
  gapInfo_.processGapStr(gapRight_, gapInfo_.gapRightRefOpen_, gapInfo_.gapRightRefExtend_);
  gapInfo_.setIdentifer();
  //
  local_ = false;
  generalMatch_ = 2;
  generalMismatch_ = -2;
  scoring_ = substituteMatrix::createDegenScoreMatrix(1,-1);
	degenScoring_ = false;
	caseInsensitiveScoring_ = false;
	lessNScoring_ = false;
  //
  colOpts_.alignOpts_.countEndGaps_ = false;
  colOpts_.iTOpts_.weighHomopolyer_ = true;

  qualThres_ = "20,15";
  qScorePars_.primaryQual_ = 20;
  qScorePars_.secondaryQual_ = 15;
  qScorePars_.qualThresWindow_ = 2;

  colOpts_.alignOpts_.eventBased_ = true;
  //
  alnInfoDirName_ = "";
  outAlnInfoDirName_ = "";
  writingOutAlnInfo_ = false;

  //
  colOpts_.kmerOpts_.runCutOffString_ = ".2%,10";
  colOpts_.kmerOpts_.runCutOff_ = 10;


  colOpts_.kmerOpts_.kLength_ = 9;
  colOpts_.kmerOpts_.kmersByPosition_ = true;
  expandKmerPos_ = false;
  expandKmerSize_ = 5;

  //general clustering
  colOpts_.skipOpts_.skipOnLetterCounterDifference_ = false;
  colOpts_.skipOpts_.fractionDifferenceCutOff_ = 0.05;
  colOpts_.iTOpts_.adjustHomopolyerRuns_ = false;
  colOpts_.bestMatchOpts_.findingBestMatch_= true;
  colOpts_.bestMatchOpts_.bestMatchCheck_ = 10;

  //quality filtering
  qFilPars_.checkingQWindow = false;
  qFilPars_.qualWindow_ = "50,5,20";
  qFilPars_.qualityWindowLength_ = 50;
  qFilPars_.qualityWindowStep_ = 5;
  qFilPars_.qualityWindowThres_ = 20;

  qFilPars_.checkingQFrac_ = false;
  qFilPars_.qualCheck_ = 30;
  qFilPars_.qualCheckCutOff_ = 0.75;

}

}  // namespace njhseq
