#pragma once
//
//  ManipulateTableRunner.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/16/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
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


