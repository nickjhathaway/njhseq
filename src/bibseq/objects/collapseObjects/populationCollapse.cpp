//
//  populationCollapse.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/2/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "populationCollapse.hpp"
namespace bibseq {
namespace collapse {
void populationCollapse::popCluster(
    collapser &collapserObj, std::map<int, std::vector<double>> iteratorMap,
    const std::string &sortBy, aligner &alignerObj) {

  collapsed_.clusters_ = collapserObj.collapseCluster(
      input_.clusters_, iteratorMap, sortBy, alignerObj);
  renameClusters("cumulativeFraction");
  updateCollapsedInfos(true);
}
void populationCollapse::writeFinal(const std::string &outDirectory,
                                    const std::string &outFormat,
                                    bool overWrite,
                                    bool exitOnFailureToWrite) const {
  collapsed_.writeClusters(outDirectory + populationName_, outFormat, overWrite,
                           exitOnFailureToWrite);
}
void populationCollapse::renameClusters(const std::string &sortBy) {
  readVecSorter::sortReadVector(collapsed_.clusters_, sortBy);
  renameReadNames(collapsed_.clusters_, populationName_, true, true);
  readVec::allUpdateName(collapsed_.clusters_);
}
void populationCollapse::writeFinalInitial(const std::string &outDirectory)
    const {
  for (const auto &clus : collapsed_.clusters_) {
    clus.writeOutClusters(outDirectory);
  }
}
void populationCollapse::renameToOtherPopNames(
    const std::vector<readObject> &previousPop, aligner &alignerObjIn) {
  // set up scoreing so gaps don't cost as much and end gaps don't cost anything
  gapScoringParameters refGapScore(3.0, 1, 0.0, 0.0, 0.0, 0.0);
  //auto currentGapScores = alignerObj.getGapScoring();
  //alignerObj.setGapScoring(refGapScore);
  aligner alignerObj(1000, refGapScore, substituteMatrix(2,-2));
  for (auto &clus : collapsed_.clusters_) {
    double bestScore = 0;
    uint32_t bestRefPos = 0;

    for (const auto &refPos : iter::range(len(previousPop))) {
      alignerObj.alignVec(previousPop[refPos], clus, false);
      if (alignerObj.parts_.score_ > bestScore) {
        bestScore = alignerObj.parts_.score_;
        bestRefPos = refPos;
      }
    }
    clus.seqBase_.name_ = previousPop[bestRefPos].seqBase_.name_;
    clus.updateName();
  }
  // set the scoring back to what it was
  //alignerObj.setGapScoring(currentGapScores);
}
}  // namespace collapse
}  // namespace bib
