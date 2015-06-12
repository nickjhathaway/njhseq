#include "simulationCommon.hpp"
#include "bibseq/common.h"

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

namespace simulation {

std::array<double, 100> makeQualErrorArr() {
  std::array<double, 100> arr;
  for (auto i : iter::range(100)) {
    arr[i] = std::pow(10.0, (-i / 10.0));
  }
  return arr;
}
//randomGenerator.hpp
const std::array<double, 100> constants::QualErrorArr = makeQualErrorArr();

}  // simulation
}  // bib
