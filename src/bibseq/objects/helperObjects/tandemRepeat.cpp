#include "tandemRepeat.hpp"
#include <iostream>

namespace bibseq {

void tandemRepeat::outPutInfo(std::ostream& out) const {
  out << "Number of repeats: " << numberOfRepeats << std::endl;
  out << "Start: " << startPos << " stop: " << stopPos << std::endl;
  out << "Alignment score: " << alignScore << std::endl;
  out << "Tandem: " << repeat << std::endl;
  out << "Size of consensus: " << repeat.size() << std::endl;
}
void tandemRepeat::outPutInfoFormated(std::ostream& out,
                                      const std::string& delim,
                                      bool first) const {
  // header is seq numberOfRepeats
  if (first) {
    out << "seq" << delim << "#ofRepeats" << delim << "seqSize" << delim << "alignScore" << delim << "start" << delim << "stop" << std::endl;
  }
  out << repeat << delim << numberOfRepeats << delim << repeat.size() << delim
      << alignScore;
  out << delim << startPos << delim << stopPos << std::endl;
}
void tandemRepeat::outPutInfoFormated(std::ostream& out,const std::string & name,
                                      const std::string& delim,
                                      bool first) const {
  // header is seq numberOfRepeats
  if (first) {
    out << "name" << delim << "seq" << delim << "#ofRepeats" << delim << "seqSize" << delim << "alignScore" << delim << "start" << delim << "stop" << std::endl;
  }
  out << name << delim << repeat << delim << numberOfRepeats << delim << repeat.size() << delim
      << alignScore;
  out << delim << startPos << delim << stopPos << std::endl;
}
}  // namespace bib
