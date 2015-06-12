#include <bibseq/alignment/alignerUtils/comparison.hpp>

namespace bibseq {
comparison::comparison()
    : oneBaseIndel_(0.0),
      twoBaseIndel_(0.0),
      largeBaseIndel_(0.0),
      hqMismatches_(0),
      lqMismatches_(0),
      lowKmerMismatches_(0),
			highQualityMatches_(0),
			lowQualityMatches_(0){}

void comparison::resetCounts() {
  oneBaseIndel_ = 0;
  twoBaseIndel_ = 0;
  largeBaseIndel_ = 0;
  hqMismatches_ = 0;
  lqMismatches_ = 0;
  lowKmerMismatches_ = 0;

  highQualityMatches_ = 0;
  lowQualityMatches_ = 0;

  distances_.reset();
}



bool comparison::passErrorProfile(const comparison& generatedError) const {
  return (oneBaseIndel_ >= generatedError.oneBaseIndel_ &&
          twoBaseIndel_ >= generatedError.twoBaseIndel_ &&
          largeBaseIndel_ >= generatedError.largeBaseIndel_ &&
          hqMismatches_ >= generatedError.hqMismatches_ &&
          lqMismatches_ >= generatedError.lqMismatches_ &&
          lowKmerMismatches_ >= generatedError.lowKmerMismatches_);
}

bool comparison::passIdThreshold(const comparison& generatedError) const {
  return generatedError.distances_.eventBasedIdentity_ >= distances_.percentIdentity_;
}

void comparison::printErrors(std::ostream& out) const {
  out << oneBaseIndel_ << "\t" << twoBaseIndel_ << "\t" << largeBaseIndel_
      << "\t" << hqMismatches_ << "\t" << lqMismatches_ << "\t"
      << lowKmerMismatches_ << std::endl;
}
void comparison::printDescription(std::ostream& out, bool deep ) const {
  out << "errorProfile{" << std::endl << "oneBaseIndel_:" << oneBaseIndel_
      << std::endl << "twoBaseIndel_:" << twoBaseIndel_ << std::endl
      << "largeBaseIndel_:" << largeBaseIndel_ << std::endl
      << "hqMismatches_:" << hqMismatches_ << std::endl
      << "lqMismatches_:" << lqMismatches_ << std::endl
      << "lowKmerMismatches_:" << lowKmerMismatches_ << std::endl;
  out << "}" << std::endl;
}

}  // namespace bibseq
