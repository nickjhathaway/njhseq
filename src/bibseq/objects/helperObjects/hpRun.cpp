#include "hpRun.hpp"

namespace bibseq {

std::string hpRun::getStringInfo() {
  std::stringstream ans;
  ans << pos << "\t" << hpPosition << "\t" << base << "\t" << runSize << "\t"
      << vectorToString(quals, ",") << "\t" << vectorMean(quals);
  return ans.str();
}
}  // namespace bib
