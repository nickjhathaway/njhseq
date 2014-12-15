//
//  tickTock.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 2/12/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "timeTracker.hpp"

namespace bibseq {
namespace simulation {

double timeTracker::getRunTime() {

	return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start_).count() /static_cast<double>( std::chrono::high_resolution_clock::period::den);
}
std::string timeTracker::getStringRunTime(bool wordy, int decPlaces) {
  return getTimeFormat(getRunTime(), true, decPlaces);
}

void timeTracker::print(const std::string& pre, std::ostream& out,
                   uint32_t indentAmount, int decPlaces) {
  out << std::string(indentAmount, '\t');
  if (pre != "") {
    out << pre << " ";
  }
  out << "(" << getStringRunTime(true, decPlaces) << ")" << std::endl;
}
}  // simulation
}  // bib
