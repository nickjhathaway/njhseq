#pragma once
//
//  bestDistGraph.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/23/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/alignment.h"
#include <bibcpp/graphics/color.hpp>
namespace bibseq {

class bestDistGraph {

 public:
  template <typename T>
  bestDistGraph(const std::vector<T>& reads, aligner& alignerObj, bool local,
                const std::string& otuName)
      : otuName_(otuName) {

    //std::cout << otuName << std::endl;
    double largestNode = 0;
    for (const auto& read : reads) {
      nodes_.push_back(node(readObject(read.seqBase_)));
      if (read.seqBase_.cnt_ > largestNode) {
        parentNode_ = nodes_.size() - 1;
        largestNode = read.seqBase_.cnt_;
      }
    }
    for (const auto& first : iter::range(nodes_.size())) {
      for (const auto& second : iter::range(nodes_.size())) {
        if (first == second) {
          continue;
        }
        alignerObj.alignVec(nodes_[first].read_, nodes_[second].read_, local);
        alignerObj.scoreAlignment(false);
        alignerObj.profilePrimerAlignment(nodes_[first].read_,
                                          nodes_[second].read_, true);
        ++bestDists_[alignerObj.comp_.hqMismatches_];
        uint32_t numOfGappedBases = 0;
        for(const auto & g : alignerObj.alignmentGaps_){
        	numOfGappedBases+= g.second.size_;
        }
        nodes_[first].children_[alignerObj.comp_.hqMismatches_].emplace_back(
            edge(second, alignerObj.comp_.distances_.ownDistance_, alignerObj.comp_.distances_.ownGapDistance_,
                 alignerObj.comp_.hqMismatches_, alignerObj.alignmentGaps_.size(),
                 numOfGappedBases, alignerObj.comp_.distances_.percentIdentity_,
                 alignerObj.comp_.distances_.percentageGaps_ > 0,
                 alignerObj.comp_.distances_.percentIdentity_ < 1.0));
      }
    }
  }
  struct edge {
    edge() {}
    edge(uint32_t childNodePos, double distance, double gapDistance,
         uint32_t misMatches, uint32_t numGaps, uint32_t numGappedBases, double identity, bool differsByIndels,
         bool differsByMismatches)
        : childNodePos_(childNodePos),
          distance_(distance),
          gapDistance_(gapDistance),
          misMatches_(misMatches),
          numGaps_(numGaps),
          numGappedBases_(numGappedBases),
          identity_(identity),
          differsByIndels_(differsByIndels),
          differsByMismatches_(differsByMismatches) {}
    uint32_t childNodePos_;
    double distance_;
    double gapDistance_;
    uint32_t misMatches_;
    uint32_t numGaps_;
    uint32_t numGappedBases_;
    double identity_;
    bool differsByIndels_;
    bool differsByMismatches_;
  };

  struct node {
    node() : added_(false) {}
    node(const readObject& read) : read_(read), added_(false) {}
    readObject read_;
    std::map<uint32_t, std::vector<edge>> children_;
    bool added_;
  };
  // memebers
  std::vector<node> nodes_;
  std::map<uint32_t, uint32_t> bestDists_;
  std::string otuName_;
  uint32_t parentNode_;
  std::unordered_map<std::string, std::vector<std::string>> mainClusterNames_;
  // funciotns
  bool allHaveBeenAdded();
  void printOutBestMatchInfos(std::ostream& out);
  void printOutBestMatch(
      std::ostream& out, std::ostream& gvOut, uint32_t nodePos,
      std::unordered_map<std::string, bool>& alreadyAdded, bool best,
      uint32_t misDist, bool recursive,
      const std::unordered_map<std::string, bib::color>& colsForName);

  // void printOutBestMatchForDist(std::ostream & out, std::ostream & gvOut,
  // uint32_t nodePos,
  //                            std::unordered_map<std::string,
  //                          bool> & alreadyAdded, uint32_t misDist);
  void createDotBestConnectedFile(
      const std::string& workingDir,
      const std::unordered_map<std::string, bib::color>& colorsForName,
      bool printTextFile);
  void printOutAllInfos(std::ostream& out);
  void printOutGV(std::ostream& out);
};

}  // bibseq

#ifndef NOT_HEADER_ONLY
#include "bestDistGraph.cpp"
#endif
