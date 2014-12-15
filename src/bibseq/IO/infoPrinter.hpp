#pragma once
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
      const std::string& fileName, bool checkingExpected, bool protein);
};

}  // bib

#ifndef NOT_HEADER_ONLY
#include "infoPrinter.cpp"
#endif