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
/*

 * QualScorePars.cpp
 *
 *  Created on: Feb 7, 2016
 *      Author: nick
 */


#include "QualScorePars.hpp"

namespace njhseq {
QualScorePars::QualScorePars():primaryQual_(20), secondaryQual_(15), qualThresWindow_(
		2){}
QualScorePars::QualScorePars(uint8_t primaryQual, uint8_t secondaryQual,
		uint32_t qualThresWindow) :
		primaryQual_(primaryQual), secondaryQual_(secondaryQual), qualThresWindow_(
				qualThresWindow) {
}



}  // namespace njhseq

