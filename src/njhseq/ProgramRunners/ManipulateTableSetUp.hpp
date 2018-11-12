#pragma once
//
//  ManipulateTableSetUp.hpp
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

#include "njhseq/objects/dataContainers/tables/TableIOOpts.hpp"
#include <njhcpp/progutils.h>

namespace njhseq {

class ManipulateTableSetUp : public njh::progutils::ProgramSetUp {

 public:
  // constructors
  ManipulateTableSetUp(int argc, char *argv[]);
  ManipulateTableSetUp(const njh::progutils::CmdArgs &inputCommands);
  // common defaults defaults

  TableIOOpts ioOptions_;
  std::string directoryName_;
  std::string sortByColumn_;
  bool decending_;
  bool addHeader_ = false;
  bool verbose_ = false;
  bool debug_ = false;

	void initializeDefaults();



  // setUps for programs
  void setUpCatOrganized();
  void setUpChangeDelim();
  void setUpSortTable();
  void setUpExtract(bool &extractingByMultipleNames, VecStr &columns,
                    bool &extractingByColumnElements, std::string &column,
                    std::string &element, std::string &elementContains,
                    std::string &elementStartsWith,
                    bool &extractingColumnsLoose, std::string &colStartsWith,
                    std::string &colContains, bool &getUniqueRows);

  void setUpTrimContent(std::string &trimAt);
  void setUpGetStats(std::string &trimAt);
  void setUpSplitTable(std::string &column, bool &getStatsInstead,
                       bool &splitingLoose, std::string &splitLooseOccurrence);
  void setUpAggregateTable(std::string &functionName, bool &advance,
                           std::string &columnName);

  void setUpPivotTable(std::string &columName, std::string &matchColumn);

  void setUpRBind(std::string &contains, bool &recursive);
  void setUpCBind(std::string &contains, bool &recursive);
  // defaults setter for the various programs
  bool processVerbose();
  bool processDebug();
  bool processFileName(bool required = true);
  void processNonRquiredDefaults();
  bool processSorting(bool required = false);
  void processWriteOutOptions();
  bool processDefaultProgram(bool fileRequired = true);
  void processDirectoryOutputName(const std::string &defaultName,
                                  bool mustMakeDirectory);
  void printDefaultOptinalOptions(std::ostream &out);
};
}  // namespace njhseq


