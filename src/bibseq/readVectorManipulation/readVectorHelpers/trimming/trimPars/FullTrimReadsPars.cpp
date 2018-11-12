/*
 * FullTrimReadsPars.cpp
 *
 *  Created on: Dec 13, 2017
 *      Author: nick
 */
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

#include "FullTrimReadsPars.hpp"

namespace njhseq {

FullTrimReadsPars::FullTrimReadsPars(){
	tSeqPars_.includeSequence_ = false;
	tSeqPars_.sequenceToLowerCase_ = false;
	tSeqPars_.removePreviousSameBases_ = false;
	allowableErrors.distances_.query_.coverage_ = .75;
}

void FullTrimReadsPars::initForKSharedTrim(){
  allowableErrors.distances_.query_.coverage_ = .50;
  allowableErrors.hqMismatches_ = 2;
  allowableErrors.lqMismatches_ = 2;
  allowableErrors.oneBaseIndel_ = 2;
  allowableErrors.twoBaseIndel_ = 2;
  allowableErrors.largeBaseIndel_ = 2;
}


}  // namespace njhseq



