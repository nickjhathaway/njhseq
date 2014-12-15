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

#include "alignerUtils.hpp"
#include <bibcpp/bashUtils.h>
namespace bibseq {
namespace cto = bib::bashCT;
runningParameters::runningParameters(const std::vector<double>& parameter, uint32_t iterNumber,
		uint32_t clusterCount) {

  stopCheck_ = parameter[0];
  if (stopCheck_ < 1.00 && stopCheck_ > 0.00) {
    stopCheck_ = (int)(stopCheck_ * clusterCount);
  } else if (stopCheck_ < 0.00) {
    stopCheck_ = clusterCount;
  }
  smallCheckStop_ = parameter[1];
  iterNumber_ = iterNumber;
  // error parameters
  errors_.oneBaseIndel_ = parameter[2];
  errors_.twoBaseIndel_ = parameter[3];
  errors_.largeBaseIndel_ = parameter[4];
  errors_.hqMismatches_ = parameter[5];
  errors_.lqMismatches_ = parameter[6];
  if (parameter.size() > 7) {
    errors_.lowKmerMismatches_ = parameter[7];
  }
}
std::string centerText(const std::string& text, uint32_t maxWidth) {
	uint32_t spaceWidth = maxWidth - text.size();
	uint32_t front = round(spaceWidth/2.0);
	uint32_t back = spaceWidth - front;
  return std::string(front, ' ') + text + std::string(back, ' ');
}

void runningParameters::printIterInfo(std::ostream & out, bool colorFormat){
	std::cout << cto::bold
			<< cto::addBGColor(145) << " Iteration "
			<< cto::addBGColor(188) << " StopAfter "
			<< cto::addBGColor(145) << " SizeCutOff "
			<< cto::addBGColor(188) << " 1BaseIndels "
			<< cto::addBGColor(145) << " 2BaseIndels "
			<< cto::addBGColor(188) << " >2BaseIndels "
			<< cto::addBGColor(145) << " HQMismatches "
			<< cto::addBGColor(188) << " LQMismatches "
			<< cto::addBGColor(145) << " LKMismatches"  << cto::reset << std::endl;
	std::cout
	<< cto::addBGColor(145) << " " << centerText(to_string(iterNumber_), 9) << " "
	<< cto::addBGColor(188) << " " << centerText(to_string(stopCheck_), 9) << " "
	<< cto::addBGColor(145) << " " << centerText(to_string(smallCheckStop_), 10) << " "
	<< cto::addBGColor(188) << " " << centerText(to_string(errors_.oneBaseIndel_), 11) << " "
	<< cto::addBGColor(145) << " " << centerText(to_string(errors_.twoBaseIndel_), 11) << " "
	<< cto::addBGColor(188) << " " << centerText(to_string(errors_.largeBaseIndel_), 12) << " "
	<< cto::addBGColor(145) << " " << centerText(to_string(errors_.hqMismatches_), 12) << " "
	<< cto::addBGColor(188) << " " << centerText(to_string(errors_.lqMismatches_), 12) << " "
	<< cto::addBGColor(145) << " " << centerText(to_string(errors_.lowKmerMismatches_), 12)
	<< cto::reset << std::endl;
}

}  // namespace bib
