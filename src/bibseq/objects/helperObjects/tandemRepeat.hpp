#pragma once
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
