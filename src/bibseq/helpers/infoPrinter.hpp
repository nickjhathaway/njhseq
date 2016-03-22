#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//  infoPrinter.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/8/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//
#include "bibseq/objects/collapseObjects.h"
namespace bibseq {
class infoPrinter {

 public:
  static void printSampleCollapseInfo(
      const std::map<std::string, collapse::sampleCollapse>& sampCollapses,
      bool checkingExpected, const std::string fileName,
      const collapse::populationCollapse& popCollapse, bool population);

  static void printPopulationCollapseInfo(
      const collapse::populationCollapse& popCollapse,
      const std::string& fileName, bool checkingExpected);

  static void printInfoForSamps(
      const std::map<std::string, collapse::sampleCollapse>& sampCollapses,
      bool checkingExpected, std::ostream & sampOut, std::ostream & popOut,
      const collapse::populationCollapse& popCollapse, bool population, const VecStr & forSamps, const std::string & groupName,
			const std::set<uint32_t> & otherPopPositions);

};

}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "infoPrinter.cpp"
#endif
