//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#pragma once

#include <iostream>
#include "bibseq/utils.h"
namespace bibseq {
namespace simulation {

class timeTracker {
 public:
  // Constructors
	timeTracker() : prefix_(""), printAtDeath_(true), fileName_("") {
		start_ = std::chrono::high_resolution_clock::now();
	}

	timeTracker(const std::string s)
      : prefix_(s), printAtDeath_(true), fileName_("") {
		start_ = std::chrono::high_resolution_clock::now();
	}

	timeTracker(const std::string s, bool printAtDeath)
      : prefix_(s), printAtDeath_(printAtDeath), fileName_("") {
		start_ = std::chrono::high_resolution_clock::now();
	}

	timeTracker(const std::string s, bool printAtDeath, const std::string& fileName)
      : prefix_(s), printAtDeath_(printAtDeath), fileName_(fileName) {
		start_ = std::chrono::high_resolution_clock::now();
	}

	timeTracker(const std::string s, const std::string& fileName)
      : prefix_(s), printAtDeath_(true), fileName_(fileName) {
		start_ = std::chrono::high_resolution_clock::now();
	}
  // deconstructor
  ~timeTracker() {
    if (printAtDeath_) {
      if (fileName_ != "") {
        std::ofstream outFile;
        outFile.open(fileName_, std::ios::app);
        print(prefix_, outFile, defaultIndent_, 6);
      } else {
        print(prefix_, std::cout, defaultIndent_, 6);
      }
    }
  }
  // members
  std::chrono::time_point<std::chrono::high_resolution_clock> start_;
  std::string prefix_;
  bool tab_ = false;
  bool printAtDeath_;
  std::string fileName_;
  uint32_t defaultIndent_ = 0;
  // functions

  double getRunTime();
  std::string getStringRunTime(bool wordy, int decPlaces);
  void print(const std::string& pre, std::ostream& out, uint32_t indentAmount,
             int decPlaces);
};
}  // simulation
}  // bib

#ifndef NOT_HEADER_ONLY
#include "timeTracker.cpp"
#endif
