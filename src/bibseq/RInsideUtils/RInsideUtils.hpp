#pragma once
//
//  RInsideUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 2/18/14.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include <RInside.h>
/////dealing with primers
// find super cluster from searching sub cluster
namespace bibseq {

/*
  void installLib (const std::string & lib,
  RInside & rSession);
  std::unordered_map<std::string, std::vector<double>> getModelQualError(RInside
  & R, const std::vector<double> & qualVec,
  const std::vector<double> & errorVec,
  const std::vector<double> & weightVec,
  uint32_t qualStart,
  uint32_t qualStop,
  double qualStep,
  const std::string & model);
  std::unordered_map<std::string, std::vector<double>> getModelQualError(RInside
  & R, const std::unordered_map<std::string, std::vector<double>>& data,
  uint32_t qualStart,
  uint32_t qualStop,
  double qualStep,
  const std::string & model);
  std::unordered_map<double, double> createLikelihoodMap(const
  std::unordered_map<std::string, std::vector<double>> & predOut);
*/

std::unordered_map<std::string, std::vector<double>> getModelQualError(
    const std::vector<double>& qualVec, const std::vector<double>& errorVec,
    const std::vector<double>& weightVec, uint32_t qualStart, uint32_t qualStop,
    double qualStep, const std::string& model);
std::unordered_map<std::string, std::vector<double>> getModelQualError(
    const std::unordered_map<std::string, std::vector<double>>& data,
    uint32_t qualStart, uint32_t qualStop, double qualStep,
    const std::string& model);
std::unordered_map<double, double> createLikelihoodMap(
    const std::unordered_map<std::string, std::vector<double>>& predOut);

class ownRInside {
 private:
  static std::unique_ptr<RInside> rSession_;

 public:
  ownRInside();
  ~ownRInside();

  RInside& get() { return *rSession_; }
  // installing libraries
  void installLib(const std::string& lib);
  void openPdfDevice(const std::string & pdfFilename, double width, double height);
  void closeCurrentDevice();
  // multiple commands
  void multipleParseEvalQ(const VecStr& cmds);
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "RInsideUtils.cpp"
#endif
