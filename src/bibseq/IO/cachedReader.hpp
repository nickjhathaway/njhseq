#pragma once
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
