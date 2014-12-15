//
//  randomGenerator.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "randomGenerator.hpp"

namespace bibseq {
std::map<int32_t, int32_t, std::greater<int>> generate(randomGenerator& gen,
                                               uint32_t start, uint32_t stop,
                                               uint32_t num, bool verbose) {
  std::map<int32_t, int32_t, std::greater<int32_t>> counts;
  for (auto i : iter::range(start, stop)) {
    counts[i] = 0;
  }
  auto randNums = gen.unifRandVector(start, stop, num);
  for (auto i : randNums) {
    ++counts[i];
  }
  if (verbose) {
    for (const auto& c : counts) {
      std::cout << c.first << " : " << c.second << std::endl;
    }
    std::cout << std::endl;
  }
  return counts;
}
}
