#pragma once
//
//  ManipulateTableRunner.hpp
//
//  Created by Nicholas Hathaway on 10/16/13.
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

#include "ManipulateTableSetUp.hpp"

namespace njhseq {

class ManipulateTableRunner : public njh::progutils::ProgramRunner {
 public:
  ManipulateTableRunner();

  static int catOrganized(const njh::progutils::CmdArgs & inputCommands);

  static int addColumn(const njh::progutils::CmdArgs & inputCommands);
  static int changeDelim(const njh::progutils::CmdArgs & inputCommands);
  static int sortTable(const njh::progutils::CmdArgs & inputCommands);

  static int tableExtractColumns(const njh::progutils::CmdArgs & inputCommands);
  static int tableExtractElementsWithPattern(const njh::progutils::CmdArgs & inputCommands);
  static int tableExtractElementsStartingWith(const njh::progutils::CmdArgs & inputCommands);
  static int tableExtractColumnsStartsWith(const njh::progutils::CmdArgs & inputCommands);
  static int tableExtractColumnsWithPattern(const njh::progutils::CmdArgs & inputCommands);
  static int tableExtractCriteria(const njh::progutils::CmdArgs & inputCommands);

  static int trimContent(const njh::progutils::CmdArgs & inputCommands);
  static int getStats(const njh::progutils::CmdArgs & inputCommands);
  static int splitTable(const njh::progutils::CmdArgs & inputCommands);
  static int aggregateTable(const njh::progutils::CmdArgs & inputCommands);
  static int pivotTable(const njh::progutils::CmdArgs & inputCommands);
  static int rBind(const njh::progutils::CmdArgs & inputCommands);
  static int cBind(const njh::progutils::CmdArgs & inputCommands);

  static int countColumn(const njh::progutils::CmdArgs & inputCommands);
  static int countRowLengths(const njh::progutils::CmdArgs & inputCommands);
  static int extractColumnElementLength(const njh::progutils::CmdArgs & inputCommands);
  static int printCol(const njh::progutils::CmdArgs & inputCommands);

  static int splitColumnContainingMeta(const njh::progutils::CmdArgs & inputCommands);

  static int roughHistogramOfColumn(const njh::progutils::CmdArgs & inputCommands);

};
}  // namespace njhseq


