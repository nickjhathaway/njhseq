#include "gaps.hpp"
#include <iostream>
#include <sstream>
#include "bibseq/utils.h"

namespace bibseq {

gap::gap(uint32_t startP,
		const std::string& seq,
		uint32_t firstQual,
		bool ref)
    : startPos_(startP),
			size_(seq.size()),
			gapedSequence_(seq),
			qualities_{firstQual},
			ref_(ref)
      {}


std::string gap::outputGapInfoSingleLine() const {
  std::stringstream ret;
  if (ref_) {
    ret << "ref"
        << "\t";
  } else {
    ret << "read"
        << "\t";
  }
  ret << startPos_ << "\t" << gapedSequence_ << "\t" << vectorToString(qualities_, ",") << "\t"
      << size_ << "\t" << homoploymerScore_;
  return ret.str();
}
}  // namespace bibseq
