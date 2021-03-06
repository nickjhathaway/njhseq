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
#include "gapScoring.hpp"

namespace njhseq {
// Constructors
gapScoringParameters::gapScoringParameters(int32_t gOpen, int32_t gExtend, int32_t gLeftOpen, int32_t gLeftExtend,
                     int32_t gRightOpen, int32_t gRightExtend)
    : gapOpen_(gOpen),
      gapExtend_(gExtend),
      gapRightQueryOpen_(gRightOpen),
      gapRightQueryExtend_(gRightExtend),
      gapRightRefOpen_(gRightOpen),
      gapRightRefExtend_(gRightExtend),

      gapLeftQueryOpen_(gLeftOpen),
      gapLeftQueryExtend_(gLeftExtend),
      gapLeftRefOpen_(gLeftOpen),
      gapLeftRefExtend_(gLeftExtend)  {
  setIdentifer();
};

gapScoringParameters::gapScoringParameters(int32_t gOpen, int32_t gExtend,
		int32_t gLeftOpenQuery, int32_t gLeftExtendQuery,
		int32_t gLeftOpenRef, int32_t gLeftExtendRef,
	int32_t gRightOpenQuery, int32_t gRightExtendQuery,
	int32_t gRightOpenRef, int32_t gRightExtendRef)    : gapOpen_(gOpen),
      gapExtend_(gExtend),
      gapRightQueryOpen_(gRightOpenQuery),
      gapRightQueryExtend_(gRightExtendQuery),
      gapRightRefOpen_(gRightOpenRef),
      gapRightRefExtend_(gRightExtendRef),

      gapLeftQueryOpen_(gLeftOpenQuery),
      gapLeftQueryExtend_(gLeftExtendQuery),
      gapLeftRefOpen_(gLeftOpenRef),
      gapLeftRefExtend_(gLeftExtendRef)  {
  setIdentifer();
};

gapScoringParameters::gapScoringParameters()
    : gapOpen_(7),
      gapExtend_(1),
      gapRightQueryOpen_(7),
      gapRightQueryExtend_(1),
      gapRightRefOpen_(7),
      gapRightRefExtend_(1),

      gapLeftQueryOpen_(7),
      gapLeftQueryExtend_(1),
      gapLeftRefOpen_(7),
      gapLeftRefExtend_(1)  {
  setIdentifer();
}
gapScoringParameters::gapScoringParameters(int32_t gapOpen, int32_t gapExtend)
    : gapOpen_(gapOpen),
      gapExtend_(gapExtend),
      gapRightQueryOpen_(gapOpen),
      gapRightQueryExtend_(gapExtend),
      gapRightRefOpen_(gapOpen),
      gapRightRefExtend_(gapExtend),

      gapLeftQueryOpen_(gapOpen),
      gapLeftQueryExtend_(gapExtend),
      gapLeftRefOpen_(gapOpen),
      gapLeftRefExtend_(gapExtend)  {
  setIdentifer();
}
gapScoringParameters::gapScoringParameters(const std::string& gapAll) {
  processGapStr(gapAll, gapOpen_, gapExtend_);
  processGapStr(gapAll, gapRightQueryOpen_, gapRightQueryExtend_);
  processGapStr(gapAll, gapRightRefOpen_, gapRightRefExtend_);
  processGapStr(gapAll, gapLeftQueryOpen_, gapLeftQueryExtend_);
  processGapStr(gapAll, gapLeftRefOpen_, gapLeftRefExtend_);
  setIdentifer();
}
gapScoringParameters::gapScoringParameters(const std::string& gap, const std::string& gapLeft,
                     const std::string& gapRight) {
  processGapStr(gap, gapOpen_, gapExtend_);
  processGapStr(gapRight, gapRightQueryOpen_, gapRightQueryExtend_);
  processGapStr(gapRight, gapRightRefOpen_, gapRightRefExtend_);
  processGapStr(gapLeft, gapLeftQueryOpen_, gapLeftQueryExtend_);
  processGapStr(gapLeft, gapLeftRefOpen_, gapLeftRefExtend_);
  setIdentifer();
}


// members

// functions
void gapScoringParameters::setIdentifer() {
	uniqueIdentifer_ = getIdentifer();
}

std::string gapScoringParameters::getIdentifer() const {
	std::stringstream tempStream;
	tempStream << gapOpen_ << "," << gapExtend_ << "," << gapLeftQueryOpen_ << ","
			<< gapLeftQueryExtend_ << "," << gapLeftRefOpen_ << ","
			<< gapLeftRefExtend_ << "," << gapRightQueryOpen_ << ","
			<< gapRightQueryExtend_ << "," << gapRightRefOpen_ << ","
			<< gapRightRefExtend_;
	return tempStream.str();
}

void gapScoringParameters::writePars(std::ostream& out) const {
  out <<  getIdentifer()
      << std::endl;
}
void gapScoringParameters::processGapStr(const std::string & gapStr,
		int32_t & open, int32_t & extend) {
	auto gapToks = tokenizeString(gapStr, ",");
	if (gapToks.size() < 2) {
		std::stringstream ss;
		ss << "Gap String needs to be two numbers seperated by a comma, eg. 7,2"
				<< std::endl
				<< "first number is gapOpen pen and second is gapExtend pen"
				<< std::endl;
		ss << "Cannot process : " << gapStr << std::endl;
		throw std::runtime_error { ss.str() };
	}
	open = std::stod(gapToks[0]);
	extend = std::stod(gapToks[1]);
}

bool gapScoringParameters::operator==(const gapScoringParameters& otherPars) const {
  return (gapOpen_ == otherPars.gapOpen_ &&
          gapExtend_ == otherPars.gapExtend_ &&
					gapRightQueryOpen_ == otherPars.gapRightQueryOpen_ &&
					gapRightQueryExtend_ == otherPars.gapRightQueryExtend_ &&
					gapRightRefOpen_ == otherPars.gapRightRefOpen_ &&
					gapRightRefExtend_ == otherPars.gapRightRefExtend_ &&
					gapLeftQueryOpen_ == otherPars.gapLeftQueryOpen_ &&
					gapLeftQueryExtend_ == otherPars.gapLeftQueryExtend_&&
					gapLeftRefOpen_ == otherPars.gapLeftRefOpen_ &&
					gapLeftRefExtend_ == otherPars.gapLeftRefExtend_);
}

bool gapScoringParameters::operator!=(const gapScoringParameters& otherPars) const {
  return (gapOpen_ != otherPars.gapOpen_ ||
          gapExtend_ != otherPars.gapExtend_ ||
					gapRightQueryOpen_ != otherPars.gapRightQueryOpen_ ||
					gapRightQueryExtend_ != otherPars.gapRightQueryExtend_ ||
					gapRightRefOpen_ != otherPars.gapRightRefOpen_ ||
					gapRightRefExtend_ != otherPars.gapRightRefExtend_ ||
					gapLeftQueryOpen_ != otherPars.gapLeftQueryOpen_ ||
					gapLeftQueryExtend_ != otherPars.gapLeftQueryExtend_ ||
					gapLeftRefOpen_ != otherPars.gapLeftRefOpen_ ||
					gapLeftRefExtend_ != otherPars.gapLeftRefExtend_);
}

bool gapScoringParameters::operator>(const gapScoringParameters& otherPars) const {
  return (gapOpen_ > otherPars.gapOpen_);
}
bool gapScoringParameters::operator<(const gapScoringParameters& otherPars) const {
  return (gapOpen_ < otherPars.gapOpen_);
}


Json::Value gapScoringParameters::toJson() const {
	Json::Value ret;
	ret["class"] = "njhseq::gapScoringParameters";
	ret["gapOpen_"] = njh::json::toJson(gapOpen_);
	ret["gapExtend_"] = njh::json::toJson(gapExtend_);
	ret["gapRightQueryOpen_"] = njh::json::toJson(gapRightQueryOpen_);
	ret["gapRightQueryExtend_"] = njh::json::toJson(gapRightQueryExtend_);
	ret["gapRightRefOpen_"] = njh::json::toJson(gapRightRefOpen_);
	ret["gapRightRefExtend_"] = njh::json::toJson(gapRightRefExtend_);
	ret["gapLeftQueryOpen_"] = njh::json::toJson(gapLeftQueryOpen_);
	ret["gapLeftQueryExtend_"] = njh::json::toJson(gapLeftQueryExtend_);
	ret["gapLeftRefOpen_"] = njh::json::toJson(gapLeftRefOpen_);
	ret["gapLeftRefExtend_"] = njh::json::toJson(gapLeftRefExtend_);
	ret["uniqueIdentifer_"] = njh::json::toJson(uniqueIdentifer_);
	return ret;
}


}  // namespace njhseq
