#pragma once
//
//  ManipulateTableRunner.hpp
//
//  Created by Nicholas Hathaway on 10/16/13.
//

// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "ManipulateTableSetUp.hpp"

namespace bibseq {

class ManipulateTableRunner : public bib::progutils::ProgramRunner {
 public:
  ManipulateTableRunner();

  static int catOrganized(const bib::progutils::CmdArgs & inputCommands);

  static int addColumn(const bib::progutils::CmdArgs & inputCommands);
  static int changeDelim(const bib::progutils::CmdArgs & inputCommands);
  static int sortTable(const bib::progutils::CmdArgs & inputCommands);

  static int tableExtractColumns(const bib::progutils::CmdArgs & inputCommands);
  static int tableExtractElementsWithPattern(const bib::progutils::CmdArgs & inputCommands);
  static int tableExtractElementsStartingWith(const bib::progutils::CmdArgs & inputCommands);
  static int tableExtractColumnsStartsWith(const bib::progutils::CmdArgs & inputCommands);
  static int tableExtractColumnsWithPattern(const bib::progutils::CmdArgs & inputCommands);
  static int tableExtractCriteria(const bib::progutils::CmdArgs & inputCommands);

  static int trimContent(const bib::progutils::CmdArgs & inputCommands);
  static int getStats(const bib::progutils::CmdArgs & inputCommands);
  static int splitTable(const bib::progutils::CmdArgs & inputCommands);
  static int aggregateTable(const bib::progutils::CmdArgs & inputCommands);
  static int pivotTable(const bib::progutils::CmdArgs & inputCommands);
  static int rBind(const bib::progutils::CmdArgs & inputCommands);
  static int cBind(const bib::progutils::CmdArgs & inputCommands);

  static int countColumn(const bib::progutils::CmdArgs & inputCommands);
  static int countRowLengths(const bib::progutils::CmdArgs & inputCommands);
  static int extractColumnElementLength(const bib::progutils::CmdArgs & inputCommands);
  static int printCol(const bib::progutils::CmdArgs & inputCommands);

  static int splitColumnContainingMeta(const bib::progutils::CmdArgs & inputCommands);

  static int roughHistogramOfColumn(const bib::progutils::CmdArgs & inputCommands);

};
}  // namespace bibseq


