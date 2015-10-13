#include "mismatch.hpp"

namespace bibseq {

void mismatch::setTransitionTransverstion(){
	transition = isMismatchTransition(refBase, seqBase);
}

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
      std::cerr << "Unrecognized base " << upperBaseB << std::endl;
    }
  } else if (upperBaseA == 'C' || upperBaseA == 'T') {
    if (upperBaseB == 'G' || upperBaseB == 'A') {
      transition = false;
    } else if (upperBaseB == 'C' || upperBaseB == 'T') {
      transition = true;
    } else {
      std::cerr << "Unrecognized base " << upperBaseB << std::endl;
    }
  } else {
    std::cerr << "Unrecognized base " << upperBaseA << std::endl;
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
