#include "tandemRepeat.hpp"
#include <iostream>

#include "bibseq/utils.h"

//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
namespace bibseq {

TandemRepeat::TandemRepeat(const std::string& rep, uint32_t numberOfRep,
		int alignS, uint32_t startP, uint32_t stopP) :
		repeat_(rep), numberOfRepeats_(numberOfRep), alignScore_(alignS), startPos_(
				startP), stopPos_(stopP) {

}

void TandemRepeat::outPutInfoFormated(std::ostream& out,
		const std::string & name, const std::string& delim) const {
	out
			<< bib::conToStr(
					toVecStr(name, repeat_, numberOfRepeats_, repeat_.size(), alignScore_,
							startPos_, stopPos_), delim) << std::endl;
}

void TandemRepeat::outPutInfoFormatedHeader(std::ostream& out,
		const std::string& delim) {
	out << bib::conToStr(VecStr { "name", "seq", "NumOfRepeats", "seqSize",
			"alignScore", "start", "stop" }, delim) << std::endl;
}

}  // namespace bib
