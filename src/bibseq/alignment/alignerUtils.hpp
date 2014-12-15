#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "bibseq/utils.h"
#include "bibseq/IO/fileUtils.hpp"
#include "bibseq/alignment/substituteMatrix.hpp"
namespace bibseq {

struct gapScoringParameters {
  // Constructors
  gapScoringParameters(int32_t gOpen, int32_t gExtend, int32_t gLeftOpen, int32_t gLeftExtend,
                       int32_t gRightOpen, int32_t gRightExtend)
      : gapOpen_(gOpen),
        gapExtend_(gExtend),
        gapRightOpen_(gRightOpen),
        gapRightExtend_(gRightExtend),
        gapLeftOpen_(gLeftOpen),
        gapLeftExtend_(gLeftExtend) {
    setIdentifer();
  };
  gapScoringParameters()
      : gapOpen_(7),
        gapExtend_(1),
        gapRightOpen_(7),
        gapRightExtend_(1),
        gapLeftOpen_(7),
        gapLeftExtend_(1) {
    setIdentifer();
  }
  gapScoringParameters(int32_t gapOpen, int32_t gapExtend)
      : gapOpen_(gapOpen),
        gapExtend_(gapExtend),
        gapRightOpen_(gapOpen),
        gapRightExtend_(gapExtend),
        gapLeftOpen_(gapOpen),
        gapLeftExtend_(gapExtend) {
    setIdentifer();
  }
  gapScoringParameters(const std::string& gapAll) {
    processGapStr(gapAll, gapOpen_, gapExtend_);
    processGapStr(gapAll, gapRightOpen_, gapRightExtend_);
    processGapStr(gapAll, gapLeftOpen_, gapLeftExtend_);
    setIdentifer();
  }
  gapScoringParameters(const std::string& gap, const std::string& gapLeft,
                       const std::string& gapRight) {
    processGapStr(gap, gapOpen_, gapExtend_);
    processGapStr(gapRight, gapRightOpen_, gapRightExtend_);
    processGapStr(gapLeft, gapLeftOpen_, gapLeftExtend_);
    setIdentifer();
  }
  virtual ~gapScoringParameters(){}

  // members
  int32_t gapOpen_;
  int32_t gapExtend_;
  int32_t gapRightOpen_;
  int32_t gapRightExtend_;
  int32_t gapLeftOpen_;
  int32_t gapLeftExtend_;
  std::string uniqueIdentifer_;
  // functions
  void setIdentifer() { uniqueIdentifer_ = getIdentifer(); }
  std::string getIdentifer() const {
    std::stringstream tempStream;
    tempStream << gapOpen_ << "," << gapExtend_ << "," << gapLeftOpen_ << ","
               << gapLeftExtend_ << "," << gapRightOpen_ << ","
               << gapRightExtend_;
    return tempStream.str();
  }
  void writePars(std::ostream& out) const {
    out << gapOpen_ << "," << gapExtend_ << "," << gapLeftOpen_ << ","
        << gapLeftExtend_ << "," << gapRightOpen_ << "," << gapRightExtend_
        << std::endl;
  }
  static void processGapStr(const std::string & gapStr, int32_t & open, int32_t & extend){
    auto gapToks = tokenizeString(gapStr, ",");
    if(gapToks.size() <2){
    	std::cout << "Gap String needs to be two numbers seperated by a comma, eg. 7,2" << std::endl
    			<< "first number is gapOpen pen and second is gapExtend pen" << std::endl;
    	exit(1);
    }
    open = std::stod(gapToks[0]);
    extend = std::stod(gapToks[1]);
  }
  bool operator==(const gapScoringParameters& otherPars) const {
    return (gapOpen_ == otherPars.gapOpen_ &&
            gapExtend_ == otherPars.gapExtend_ &&
            gapRightOpen_ == otherPars.gapRightOpen_ &&
            gapRightExtend_ == otherPars.gapRightExtend_ &&
            gapLeftOpen_ == otherPars.gapLeftOpen_ &&
            gapLeftExtend_ == otherPars.gapLeftExtend_);
  }

  bool operator!=(const gapScoringParameters& otherPars) const {
    return (gapOpen_ != otherPars.gapOpen_ ||
            gapExtend_ != otherPars.gapExtend_ ||
            gapRightOpen_ != otherPars.gapRightOpen_ ||
            gapRightExtend_ != otherPars.gapRightExtend_ ||
            gapLeftOpen_ != otherPars.gapLeftOpen_ ||
            gapLeftExtend_ != otherPars.gapLeftExtend_);
  }

  bool operator>(const gapScoringParameters& otherPars) const {
    return (gapOpen_ > otherPars.gapOpen_);
  }
  bool operator<(const gapScoringParameters& otherPars) const {
    return (gapOpen_ < otherPars.gapOpen_);
  }
  virtual void printDescription(std::ostream& out, bool deep = false) const {
    out << "gapScoringParameters{" << std::endl << "gapOpen_:" << gapOpen_
        << std::endl << "gapExtend_:" << gapExtend_ << std::endl
        << "gapLeftOpen_:" << gapLeftOpen_ << std::endl
        << "gapLeftExtend_:" << gapLeftExtend_ << std::endl
        << "gapRightOpen_:" << gapRightOpen_ << std::endl
        << "gapRightExtend_:" << gapRightExtend_ << std::endl
        << "uniqueIdentifer_:" << uniqueIdentifer_ << std::endl;
    out << "}" << std::endl;
  }

};

class errorProfile {
public:
  errorProfile()
      : oneBaseIndel_(0.0),
        twoBaseIndel_(0.0),
        largeBaseIndel_(0.0),
        hqMismatches_(0),
        lqMismatches_(0),
        lowKmerMismatches_(0) {}
  void resetCounts() {
    oneBaseIndel_ = 0;
    twoBaseIndel_ = 0;
    largeBaseIndel_ = 0;
    hqMismatches_ = 0;
    lqMismatches_ = 0;
    lowKmerMismatches_ = 0;
  }
  double oneBaseIndel_;
  double twoBaseIndel_;
  double largeBaseIndel_;
  int hqMismatches_;
  int lqMismatches_;
  int lowKmerMismatches_;
  bool passErrorProfile(const errorProfile& generatedError) const {
    return (oneBaseIndel_ >= generatedError.oneBaseIndel_ &&
            twoBaseIndel_ >= generatedError.twoBaseIndel_ &&
            largeBaseIndel_ >= generatedError.largeBaseIndel_ &&
            hqMismatches_ >= generatedError.hqMismatches_ &&
            lqMismatches_ >= generatedError.lqMismatches_);
  }
  bool passErrorProfileLowKmer(const errorProfile& generatedError) const {
    return (oneBaseIndel_ >= generatedError.oneBaseIndel_ &&
            twoBaseIndel_ >= generatedError.twoBaseIndel_ &&
            largeBaseIndel_ >= generatedError.largeBaseIndel_ &&
            hqMismatches_ >= generatedError.hqMismatches_ &&
            lqMismatches_ >= generatedError.lqMismatches_ &&
            lowKmerMismatches_ >= generatedError.lowKmerMismatches_);
  }

  void printErrors(std::ostream& out) const {
    out << oneBaseIndel_ << "\t" << twoBaseIndel_ << "\t" << largeBaseIndel_
        << "\t" << hqMismatches_ << "\t" << lqMismatches_ << "\t"
        << lowKmerMismatches_ << std::endl;
  }
  virtual void printDescription(std::ostream& out, bool deep = false) const {
    out << "errorProfile{" << std::endl << "oneBaseIndel_:" << oneBaseIndel_
        << std::endl << "twoBaseIndel_:" << twoBaseIndel_ << std::endl
        << "largeBaseIndel_:" << largeBaseIndel_ << std::endl
        << "hqMismatches_:" << hqMismatches_ << std::endl
        << "lqMismatches_:" << lqMismatches_ << std::endl
        << "lowKmerMismatches_:" << lowKmerMismatches_ << std::endl;
    out << "}" << std::endl;
  }
  virtual ~errorProfile(){}
};

class runningParameters {
public:
  runningParameters() : stopCheck_(100), smallCheckStop_(0), iterNumber_(0) {}
  runningParameters(const std::vector<double>& parameter, uint32_t iterNumber,
  		uint32_t clusterCount);
  // procedure parameters
  double stopCheck_;
  uint32_t smallCheckStop_;
  uint32_t iterNumber_;
  // error parameters
  errorProfile errors_;

  void printIterInfo(std::ostream & out, bool colorFormat);
};

struct runningParametersOld {
  runningParametersOld()
      : stopCheck_(100),
        smallCheckStop(0),
        oneBaseIndel_(0),
        twoBaseIndel_(0),
        largeBaseIndel_(0),
        hqMismatches_(0),
        lqMismatches_(0),
        homopolymerScore(0) {}

  runningParametersOld(std::vector<double>& parameter, int clusterCount) {
    stopCheck_ = parameter[0];
    if (stopCheck_ < 1.00 && stopCheck_ > 0.00) {
      stopCheck_ = (int)(stopCheck_ * clusterCount);
    } else if (stopCheck_ < 0.00) {
      stopCheck_ = clusterCount;
    }
    smallCheckStop = parameter[1];
    // error parameters
    oneBaseIndel_ = parameter[2];
    twoBaseIndel_ = parameter[3];
    largeBaseIndel_ = parameter[4];
    hqMismatches_ = parameter[5];
    lqMismatches_ = parameter[6];
    homopolymerScore = parameter[7];
  }
  // procedure paramters
  double stopCheck_;
  int smallCheckStop;
  // error parameters
  int oneBaseIndel_;
  int twoBaseIndel_;
  int largeBaseIndel_;
  int hqMismatches_;
  int lqMismatches_;
  double homopolymerScore;
};

struct errorProfileOld {

 public:
  errorProfileOld() {}
  int oneBaseIndel_;
  int twoBaseIndel_;
  int largeBaseIndel_;
  int hqMismatches_;
  int lqMismatches_;
  double homopolymerScore;
  bool passErrorProfile(const errorProfileOld& generatedError) {
    return (oneBaseIndel_ > generatedError.oneBaseIndel_ &&
            twoBaseIndel_ > generatedError.twoBaseIndel_ &&
            largeBaseIndel_ > generatedError.largeBaseIndel_ &&
            hqMismatches_ > generatedError.hqMismatches_ &&
            lqMismatches_ > generatedError.lqMismatches_ &&
            homopolymerScore > generatedError.homopolymerScore);
  }
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "alignerUtils.cpp"
#endif
