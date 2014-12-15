//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "simulationCommon.hpp"
#include "bibseq/common.h"

namespace bibseq {
namespace simulation {

std::array<double, 100> makeQualErrorArr() {
  std::array<double, 100> arr;
  for (auto i : iter::range(100)) {
    arr[i] = std::pow(10.0, (-i / 10.0));
  }
  return arr;
}

const std::array<double, 100> constants::QualErrorArr = makeQualErrorArr();

}  // simulation
}  // bib
