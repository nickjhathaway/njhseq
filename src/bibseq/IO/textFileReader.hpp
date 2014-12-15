#pragma once
//
//  textFileReader.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 11/11/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/IO/fileUtils.hpp"
#include <iomanip>
#include "bibseq/objects/dataContainers/table.hpp"

namespace bibseq {

class textFileReader {

 public:
  // constructor
  textFileReader();
  textFileReader(const std::string &delimiter) : delim(delimiter) {
    if (delim == "tab") {
      delim = "\t";
    } else if (delim == "whitespace") {
      delim = " ";
    }
  }

  // read in a file
  void readFile(const std::string &filename, bool header = false);
  // file content
  table fileContent;
  void outPut(const outOptions &options);
  std::string delim;
  static std::string getFirstLine(const std::string &filename);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "textFileReader.cpp"
#endif
