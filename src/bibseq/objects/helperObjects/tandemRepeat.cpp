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
