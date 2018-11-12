#include "tandemRepeat.hpp"
#include <iostream>

#include "njhseq/utils.h"

//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
namespace njhseq {

TandemRepeat::TandemRepeat(const std::string& rep, uint32_t numberOfRep,
		int alignS, uint32_t startP, uint32_t stopP) :
		repeat_(rep), numberOfRepeats_(numberOfRep), alignScore_(alignS), startPos_(
				startP), stopPos_(stopP) {

}

uint32_t TandemRepeat::getSize() const{
	return stopPos_ - startPos_;
}

void TandemRepeat::outPutInfoFormated(std::ostream& out,
		const std::string & name, const std::string& delim) const {
	out
			<< njh::conToStr(
					toVecStr(name, repeat_, numberOfRepeats_, repeat_.size(), alignScore_,
							startPos_, stopPos_, stopPos_ - startPos_), delim) << std::endl;
}

void TandemRepeat::outPutInfoFormatedHeader(std::ostream& out,
		const std::string& delim) {
	out << njh::conToStr(VecStr { "name", "seq", "NumOfRepeats", "seqSize",
			"alignScore", "start", "stop", "fullLength"}, delim) << std::endl;
}

}  // namespace njh
