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
//  cachedReader.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 05/25/15.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"


namespace bibseq {


class cachedReader{

private:
  VecStr lineBuffer_;
  uint32_t bufferPos_ = std::numeric_limits<uint32_t>::max();
  uint32_t lineNum_ = 0;
  uint32_t bufferMax_ = 10;
  std::istream & is_;
  bool doneReadering_ = false;
public:
  cachedReader(std::istream & is);

  bool refillBuffer();

  const std::string & currentLine();

  bool setNextLine();

  bool done();

  void seek(uint64_t filePosition);

};


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "cachedReader.cpp"
#endif
