#pragma once
//
//  ManipulateTableSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/16/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/dataContainers/tables/TableIOOpts.hpp"

namespace bibseq {

class ManipulateTableSetUp : public bib::progutils::ProgramSetUp {

 public:
  // constructors
  ManipulateTableSetUp(int argc, char *argv[]);
  ManipulateTableSetUp(const bib::progutils::CmdArgs &inputCommands);
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
}  // namespace bibseq


