//
//  infoPrinter.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/8/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "infoPrinter.hpp"
namespace bibseq {

void infoPrinter::printSampleCollapseInfo(
    const std::map<std::string, collapse::sampleCollapse>& sampCollapses,
    bool checkingExpected, const std::string fileName,
    const collapse::populationCollapse& popCollapse, bool population) {
  std::ofstream infoFile;
  openTextFile(infoFile, fileName, ".txt", true, false);
  if (population) {
    infoFile
        << "Sample\tpopUID\tpopReadCntTot\tpopInputClusterCnt\thapPopFrac\thapS"
           "umSampPopFrac"
        << "\thapMeanSampFrac\thapReadFrac\thapSampCNT\thapSampFrac\thapReadCNT"
        << "\thapClusterCNT\tclusterNames\thapConsesus\t";
  }
  uint64_t maxRunCount = 0;
  for (auto& sampCollapse : sampCollapses) {
    if (maxRunCount < sampCollapse.second.collapsed_.infos_.size()) {
      maxRunCount = sampCollapse.second.collapsed_.infos_.size();
    }
  }
  infoFile << "sName\tsReadCntTot\tsInputCluster\tclusterID";
  infoFile << "\tcAveragedFrac\tcReadFrac\tcRepCnt\tcConsensus";
  infoFile << "\tReadCnt\tcChiReadCnt\tcChiClusCnt\tcChiRepCnt\tcInputNames";
  std::string templateRunSring =
      "\tRunNUM\tRNUMtotalCntExcluded\tRNUMtotalFracExcluded\tRNUMcntChiExcluded\tRNUMfracChiExcluded\tRNUM."
      "MapFrac\tRNUM.ReadCnt\tRNUM.ClusCnt\tRNUM.totalReadCnt";
  for (uint64_t i = 1; i <= maxRunCount; ++i) {
    infoFile << replaceString(templateRunSring, "NUM", std::to_string(i));
  }
  if (checkingExpected) {
    infoFile << "\tbestExpected";
  }
  infoFile << std::endl;
  for (auto& sampCollapse : sampCollapses) {
    std::cout << sampCollapse.second.sampName_ << std::endl;
    auto currentInfos = sampCollapse.second.getAllInfoMap(checkingExpected);
    for (const auto& info : currentInfos) {
      if (population) {
        std::cout << "\t" << info.first << std::endl;
        infoFile << sampCollapse.second.sampName_;
        infoFile
            << "\t"
            << popCollapse.collapsed_.clusters_
                   [popCollapse.collapsed_.subClustersPositions_.at(info.first)]
                       .getPopStandardInfo(
                            popCollapse.collapsed_.totalReadCount_,
                            popCollapse.collapsed_.numberOfClusters_,
                            len(popCollapse.collapsed_.infos_), false);
        infoFile << "\t";
      }
      infoFile << info.second << std::endl;
    }
  }
}
void infoPrinter::printPopulationCollapseInfo(
    const collapse::populationCollapse& popCollapse,
    const std::string& fileName, bool checkingExpected, bool protein) {

  std::ofstream popInfoOutFile;
  openTextFile(popInfoOutFile, fileName, ".txt", true, false);
  popInfoOutFile
      << "popUID\tpopReadCntTot\tpopInputClusterCnt\thapPopFrac\thapSumSampPopF"
         "rac"
      << "\thapMeanSampFrac\thapReadFrac\thapSampCNT\thapSampFrac\thapReadCNT"
      << "\thapClusterCNT\tclusterNames\thapConsesus";
  if (protein) {
    popInfoOutFile << "\tprotein";
  }
  if (checkingExpected) {
    popInfoOutFile << "\tbestExpected";
  }
  popInfoOutFile << std::endl;
  for (const auto& clus : iter::reverse(popCollapse.collapsed_.clusters_)) {
    popInfoOutFile << clus.getPopStandardInfo(
                          popCollapse.collapsed_.totalReadCount_,
                          popCollapse.collapsed_.numberOfClusters_,
                          len(popCollapse.collapsed_.infos_), protein);
    if (checkingExpected) {
      popInfoOutFile << "\t" << clus.expectsString;
    }
    popInfoOutFile << std::endl;
  }
}
}  // bib
