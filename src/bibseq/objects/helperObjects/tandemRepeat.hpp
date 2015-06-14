#pragma once
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
//
//  tandemRepeat.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 11/27/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include <string>

namespace bibseq {

class tandemRepeat {

 public:
  // tandem repeat info
  std::string repeat;
  int numberOfRepeats;
  int alignScore;
  int startPos;
  int stopPos;

  // constructor
  tandemRepeat() {}
  tandemRepeat(const std::string& rep, int numberOfRep, int alignS, int startP,
               int stopP)
      : repeat(rep),
        numberOfRepeats(numberOfRep),
        alignScore(alignS),
        startPos(startP),
        stopPos(stopP) {}

  // output
  void outPutInfo(std::ostream& out) const;
  void outPutInfoFormated(std::ostream& out, const std::string& delim = "\t",
                          bool first = false) const;
  void outPutInfoFormated(std::ostream& out,const std::string & name, const std::string& delim = "\t",
                            bool first = false) const;
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "tandemRepeat.cpp"
#endif
