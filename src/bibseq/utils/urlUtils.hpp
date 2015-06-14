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
//  urlUtils.hpp
//
//  Created by Nick Hathaway on 5/27/15.
//  Copyright (c) 2015 University of Massachusetts Medical School. All rights
// reserved.
//

#include "bibseq/common.h"
#include "bibseq/utils/bitSwaps.hpp"

namespace bibseq {

size_t WriteCallback(char* contents, size_t size, size_t nmemb,
                     std::ostream* stream);
std::string GetURL(const std::string url);
void GetURLStream(const std::string url, std::ostream & out);


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "urlUtils.cpp"
#endif
