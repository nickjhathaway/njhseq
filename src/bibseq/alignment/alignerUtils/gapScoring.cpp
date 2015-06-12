#include "gapScoring.hpp"

namespace bibseq {
// Constructors
gapScoringParameters::gapScoringParameters(int32_t gOpen, int32_t gExtend, int32_t gLeftOpen, int32_t gLeftExtend,
                     int32_t gRightOpen, int32_t gRightExtend)
    : gapOpen_(gOpen),
      gapExtend_(gExtend),
      gapRightOpen_(gRightOpen),
      gapRightExtend_(gRightExtend),
      gapLeftOpen_(gLeftOpen),
      gapLeftExtend_(gLeftExtend) {
  setIdentifer();
};
gapScoringParameters::gapScoringParameters()
    : gapOpen_(7),
      gapExtend_(1),
      gapRightOpen_(7),
      gapRightExtend_(1),
      gapLeftOpen_(7),
      gapLeftExtend_(1) {
  setIdentifer();
}
gapScoringParameters::gapScoringParameters(int32_t gapOpen, int32_t gapExtend)
    : gapOpen_(gapOpen),
      gapExtend_(gapExtend),
      gapRightOpen_(gapOpen),
      gapRightExtend_(gapExtend),
      gapLeftOpen_(gapOpen),
      gapLeftExtend_(gapExtend) {
  setIdentifer();
}
gapScoringParameters::gapScoringParameters(const std::string& gapAll) {
  processGapStr(gapAll, gapOpen_, gapExtend_);
  processGapStr(gapAll, gapRightOpen_, gapRightExtend_);
  processGapStr(gapAll, gapLeftOpen_, gapLeftExtend_);
  setIdentifer();
}
gapScoringParameters::gapScoringParameters(const std::string& gap, const std::string& gapLeft,
                     const std::string& gapRight) {
  processGapStr(gap, gapOpen_, gapExtend_);
  processGapStr(gapRight, gapRightOpen_, gapRightExtend_);
  processGapStr(gapLeft, gapLeftOpen_, gapLeftExtend_);
  setIdentifer();
}


// members

// functions
void gapScoringParameters::setIdentifer() { uniqueIdentifer_ = getIdentifer(); }
std::string gapScoringParameters::getIdentifer() const {
  std::stringstream tempStream;
  tempStream << gapOpen_ << "," << gapExtend_ << "," << gapLeftOpen_ << ","
             << gapLeftExtend_ << "," << gapRightOpen_ << ","
             << gapRightExtend_;
  return tempStream.str();
}
void gapScoringParameters::writePars(std::ostream& out) const {
  out << gapOpen_ << "," << gapExtend_ << "," << gapLeftOpen_ << ","
      << gapLeftExtend_ << "," << gapRightOpen_ << "," << gapRightExtend_
      << std::endl;
}
void gapScoringParameters::processGapStr(const std::string & gapStr, int32_t & open, int32_t & extend){
  auto gapToks = tokenizeString(gapStr, ",");
  if(gapToks.size() <2){
  	std::stringstream ss;
  	ss << "Gap String needs to be two numbers seperated by a comma, eg. 7,2" << std::endl
  			<< "first number is gapOpen pen and second is gapExtend pen" << std::endl;
  	ss << "Cannot process : " << gapStr << std::endl;
  	throw std::runtime_error{ss.str()};
  }
  open = std::stod(gapToks[0]);
  extend = std::stod(gapToks[1]);
}
bool gapScoringParameters::operator==(const gapScoringParameters& otherPars) const {
  return (gapOpen_ == otherPars.gapOpen_ &&
          gapExtend_ == otherPars.gapExtend_ &&
          gapRightOpen_ == otherPars.gapRightOpen_ &&
          gapRightExtend_ == otherPars.gapRightExtend_ &&
          gapLeftOpen_ == otherPars.gapLeftOpen_ &&
          gapLeftExtend_ == otherPars.gapLeftExtend_);
}

bool gapScoringParameters::operator!=(const gapScoringParameters& otherPars) const {
  return (gapOpen_ != otherPars.gapOpen_ ||
          gapExtend_ != otherPars.gapExtend_ ||
          gapRightOpen_ != otherPars.gapRightOpen_ ||
          gapRightExtend_ != otherPars.gapRightExtend_ ||
          gapLeftOpen_ != otherPars.gapLeftOpen_ ||
          gapLeftExtend_ != otherPars.gapLeftExtend_);
}

bool gapScoringParameters::operator>(const gapScoringParameters& otherPars) const {
  return (gapOpen_ > otherPars.gapOpen_);
}
bool gapScoringParameters::operator<(const gapScoringParameters& otherPars) const {
  return (gapOpen_ < otherPars.gapOpen_);
}
void gapScoringParameters::printDescription(std::ostream& out, bool deep) const {
  out << "gapScoringParameters{" << std::endl << "gapOpen_:" << gapOpen_
      << std::endl << "gapExtend_:" << gapExtend_ << std::endl
      << "gapLeftOpen_:" << gapLeftOpen_ << std::endl
      << "gapLeftExtend_:" << gapLeftExtend_ << std::endl
      << "gapRightOpen_:" << gapRightOpen_ << std::endl
      << "gapRightExtend_:" << gapRightExtend_ << std::endl
      << "uniqueIdentifer_:" << uniqueIdentifer_ << std::endl;
  out << "}" << std::endl;
}


}  // namespace bibseq
