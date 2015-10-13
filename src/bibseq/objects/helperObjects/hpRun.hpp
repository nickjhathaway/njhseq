#pragma once
//
//  hpRun.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 1/6/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include <vector>
#include "bibseq/utils.h"
namespace bibseq {

class hpRun {

 public:
  hpRun() {}
  hpRun(char hpRunBase, int position, int hpSize)
      : base(hpRunBase), pos(position), runSize(hpSize), count(1) {}
  hpRun(char hpRunBase, int position, int hpSize, int readCount)
      : base(hpRunBase), pos(position), runSize(hpSize), count(readCount) {}
  hpRun(char hpRunBase, int position, int hpSize, int readCount,
        const std::vector<int>& firstQual)
      : base(hpRunBase),
        pos(position),
        runSize(hpSize),
        count(readCount),
        quals(firstQual) {}
  hpRun(char hpRunBase, int position, int homopolyerPosition, int hpSize,
        int readCount, const std::vector<int>& firstQual)
      : base(hpRunBase),
        pos(position),
        hpPosition(homopolyerPosition),
        runSize(hpSize),
        count(readCount),
        quals(firstQual) {}
  char base;
  int pos;
  int hpPosition;
  int runSize;
  int count;
  std::vector<int> quals;
  void increaseCount(int runAdd, int qualAdd) {
    count += runAdd;
    quals.push_back(qualAdd);
  }
  hpRun& operator++() {
    ++count;
    return *this;
  }
  hpRun& operator+=(int add) {
    count += add;
    return *this;
  }
  std::string getStringInfo();
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "hpRun.cpp"
#endif
