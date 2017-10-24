#pragma once
/*
 * SeqSetUpPars.hpp
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

#include "bibseq/IO/SeqIO/SeqIOOptions.hpp"
#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/alignment/alignerUtils.h"
#include "bibseq/objects/collapseObjects/opts/CollapserOpts.hpp"

namespace bibseq {
struct QualFilteringPars {

  bool checkingQWindow = false;
  std::string qualWindow_ = "50,5,25";
  uint32_t qualityWindowLength_ = 0;
  uint32_t qualityWindowStep_ = 0;
  uint32_t qualityWindowThres_ = 0;

  bool checkingQFrac_ = false;
  uint32_t qualCheck_ = 30;
  double qualCheckCutOff_ = 0.75;

  uint32_t trimAtQualCutOff_ = 2;
  bool trimAtQual_ = false;

  Json::Value toJson() const;

};
class SeqSetUpPars {
public:
	SeqSetUpPars();
  void initializeDefaults();

  // seq read in names
  SeqIOOptions ioOptions_;
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
  SeqIOOptions refIoOptions_;

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
  std::string qualThres_;
  QualScorePars qScorePars_;

  // kmer profiling
  bool expandKmerPos_;
  uint32_t expandKmerSize_;

  //alnment caching
  std::string alnInfoDirName_;
  std::string outAlnInfoDirName_;
  bool writingOutAlnInfo_;

  //general clustering options
  CollapserOpts colOpts_;

  ChimeraOpts chiOpts_;

  //quality filtering
  QualFilteringPars qFilPars_;

};


}  // namespace bibseq


