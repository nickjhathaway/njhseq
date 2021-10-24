/*
 * numUtils.cpp
 *
 *  Created on: Jun 22, 2014
 *      Author: nickhathaway
 */
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
#include "numUtils.hpp"

namespace njhseq {






std::vector<double> getRange(double start, double stop, uint32_t num) {
  double difference = stop - start;
  double step = difference / (num - 1);
  std::vector<double> ans;
  for (const auto i : iter::range<uint32_t>(0, num)) {
    ans.emplace_back((i * step) + start);
  }
  return ans;
}


} /* namespace njh */
