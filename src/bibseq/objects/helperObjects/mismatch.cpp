//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "mismatch.hpp"

namespace bibseq {

bool mismatch::isMismatchTransition(const char& baseA, const char& baseB) {
  bool transition = false;
  // current fix for degenerative base, need a better way
  if (baseA == 'R' || baseA == 'S' || baseA == 'Y' || baseA == 'W' ||
      baseA == 'K' || baseA == 'M' || baseA == 'N') {
    return transition;
  }
  if (baseB == 'R' || baseB == 'S' || baseB == 'Y' || baseB == 'W' ||
      baseB == 'K' || baseB == 'M' || baseB == 'N') {
    return transition;
  }
  char upperBaseA = toupper(baseA);
  char upperBaseB = toupper(baseB);
  if (upperBaseA == 'G' || upperBaseA == 'A') {
    if (upperBaseB == 'G' || upperBaseB == 'A') {
      transition = true;
    } else if (upperBaseB == 'C' || upperBaseB == 'T') {
      transition = false;
    } else {
      std::cout << "Unrecognized base " << upperBaseB << std::endl;
    }
  } else if (upperBaseA == 'C' || upperBaseA == 'T') {
    if (upperBaseB == 'G' || upperBaseB == 'A') {
      transition = false;
    } else if (upperBaseB == 'C' || upperBaseB == 'T') {
      transition = true;
    } else {
      std::cout << "Unrecognized base " << upperBaseB << std::endl;
    }
  } else {
    std::cout << "Unrecognized base " << upperBaseA << std::endl;
  }
  return transition;
}

std::string mismatch::outputInfoString() const {
  std::stringstream out;

  if (transition) {
    out << "transition\t";
  } else {
    out << "transversion\t";
  }
  out << refBasePos << "\t" << refBase << "\t" << refQual << "\t"
      << vectorToString(refLeadingQual, ",") << "\t"
      << vectorToString(refTrailingQual, ",") << "\t" << seqBasePos << "\t"
      << seqBase << "\t" << seqQual << "\t"
      << vectorToString(seqLeadingQual, ",") << "\t"
      << vectorToString(seqTrailingQual, ",") << "\t" << kMerFreqByPos << "\t"
      << kMerFreq;
  return out.str();
}
Json::Value mismatch::outputJson()const{
	Json::Value ret;
	ret["refBase"] = bib::json::toJson(refBase);
	ret["refQual"] = bib::json::toJson(refQual);
	ret["refLeadingQual"] = bib::json::toJson(refLeadingQual);
	ret["refTrailingQual"] = bib::json::toJson(refTrailingQual);
	ret["refBasePos"] = bib::json::toJson(refBasePos);
	ret["seqBase"] = bib::json::toJson(seqBase);
	ret["seqQual"] = bib::json::toJson(seqQual);
	ret["seqLeadingQual"] = bib::json::toJson(seqLeadingQual);
	ret["seqTrailingQual"] = bib::json::toJson(seqTrailingQual);
	ret["seqBasePos"] = bib::json::toJson(seqBasePos);
	ret["kMerFreqByPos"] = bib::json::toJson(kMerFreqByPos);
	ret["kMerFreq"] = bib::json::toJson(kMerFreq);
	ret["transition"] = bib::json::toJson(transition);
	ret["freq"] = bib::json::toJson(freq);
	ret["frac_"] = bib::json::toJson(frac_);
	return ret;
}


}  // namespace bib
