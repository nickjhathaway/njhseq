#pragma once
//
//  otuGraph.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/15/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/alignment.h"

namespace bibseq {
class otuGraph {

 public:
  otuGraph(const std::vector<readObject>& reads, aligner& alignerObj,
           bool local, const std::string& otuName)
      : otuName_(otuName) {
    double biggestRead = 0;
    for (const auto& read : reads) {
      nodes_.push_back(node(read));
      if (read.seqBase_.cnt_ > biggestRead) {
        biggestRead = read.seqBase_.cnt_;
        parentPos_ = nodes_.size() - 1;
      } else if (read.seqBase_.cnt_ == biggestRead) {
        if (read.getAverageErrorRate() <
            nodes_[parentPos_].read_.getAverageErrorRate()) {
          parentPos_ = nodes_.size() - 1;
        }
      }
    }
    uint32_t nodePos = 0;
    for (const auto& n : nodes_) {
      double bestIdentity = 0.0;
      double distance = 0.0;
      double gapDistance = 0.0;
      uint32_t misMatches = 0;
      bool differsByIndel = false;
      bool differsByMismatch = false;
      uint32_t bestPos = 0;
      uint32_t pos = 0;
      bool foundParent = false;
      for (const auto& secondN : nodes_) {
        if (n.read_.seqBase_.name_ != secondN.read_.seqBase_.name_ &&
            n.read_.seqBase_.cnt_ < secondN.read_.seqBase_.cnt_) {
          alignerObj.alignVec(secondN.read_, n.read_, local);
          alignerObj.scoreAlignment(false);
          alignerObj.profilePrimerAlignment(secondN.read_, n.read_, true);
          if (alignerObj.comp_.distances_.percentIdentity_ > bestIdentity) {
            bestIdentity = alignerObj.comp_.distances_.percentIdentity_;
            distance = alignerObj.comp_.distances_.ownDistance_;
            gapDistance = alignerObj.comp_.distances_.ownGapDistance_;
            misMatches = alignerObj.comp_.hqMismatches_;
            bestPos = pos;
            foundParent = true;
            if (alignerObj.comp_.distances_.percentageGaps_ > 0) {
              differsByIndel = true;
            } else {
              differsByIndel = false;
            }
            if (alignerObj.comp_.distances_.percentIdentity_ < 1.0) {
              differsByMismatch = true;
            } else {
              differsByMismatch = false;
            }
          }
          if (alignerObj.comp_.distances_.percentIdentity_ == bestIdentity) {
            if (nodes_[bestPos].read_.seqBase_.cnt_ <
                secondN.read_.seqBase_.cnt_) {
              bestPos = pos;
              distance = alignerObj.comp_.distances_.ownDistance_;
              gapDistance = alignerObj.comp_.distances_.ownGapDistance_;
              misMatches = alignerObj.comp_.hqMismatches_;
              foundParent = true;
              if (alignerObj.comp_.distances_.percentageGaps_ > 0) {
                differsByIndel = true;
              } else {
                differsByIndel = false;
              }
              if (alignerObj.comp_.distances_.percentIdentity_ < 1.0) {
                differsByMismatch = true;
              } else {
                differsByMismatch = false;
              }
            }
            if (nodes_[bestPos].read_.seqBase_.cnt_ ==
                secondN.read_.seqBase_.cnt_) {
              if (secondN.read_.getAverageErrorRate() <
                  nodes_[bestPos].read_.getAverageErrorRate()) {
                bestPos = pos;
                distance = alignerObj.comp_.distances_.ownDistance_;
                gapDistance = alignerObj.comp_.distances_.ownGapDistance_;
                misMatches = alignerObj.comp_.hqMismatches_;
                foundParent = true;
                if (alignerObj.comp_.distances_.percentageGaps_ > 0) {
                  differsByIndel = true;
                } else {
                  differsByIndel = false;
                }
                if (alignerObj.comp_.distances_.percentIdentity_ < 1.0) {
                  differsByMismatch = true;
                } else {
                  differsByMismatch = false;
                }
              }
            }
          }
        }
        ++pos;
      }
      if (foundParent) {
        nodes_[bestPos].children_.emplace_back(
            edge(nodePos, distance, gapDistance, misMatches, bestIdentity,
                 differsByIndel, differsByMismatch));
      }
      ++nodePos;
    }
  }

  struct edge {
    edge() {}
    edge(uint32_t childNodePos, double distance, double gapDistance,
         uint32_t misMatches, double identity, bool differsByIndels,
         bool differsByMismatches)
        : childNodePos_(childNodePos),
          distance_(distance),
          gapDistance_(gapDistance),
          misMatches_(misMatches),
          identity_(identity),
          differsByIndels_(differsByIndels),
          differsByMismatches_(differsByMismatches) {}
    uint32_t childNodePos_;
    double distance_;
    double gapDistance_;
    uint32_t misMatches_;
    double identity_;
    bool differsByIndels_;
    bool differsByMismatches_;
  };

  struct node {
    node() {}
    node(const readObject& read) : read_(read) {}
    readObject read_;
    std::vector<edge> children_;
  };
  // memebers
  std::vector<node> nodes_;
  uint32_t parentPos_;
  std::string otuName_;
  // funciotns

  // out
  void printParent(std::ostream& out);
  void printChildren(std::ostream& out, const node& currentNode);

  void printInfo(std::ostream& out);
  void printInfoGraphViz(std::ostream& out);
  void printChildrenGV(std::ostream& out, const node& currentNode);
  void printParentLevel(std::ostream& out);
  void printChildrenLevel(std::ostream& out, double x1, double y1,
                          double currentRadians, uint32_t currentNodePos);
  void printSimpleChildrenInfo(std::ostream& out, uint32_t nodePos);
  void getNames(VecStr& names, uint32_t nodePos, std::string currentName);

  std::string appendChildName(std::string& name, uint32_t currentNodePos);
};

}  // bib

#ifndef NOT_HEADER_ONLY
#include "otuGraph.cpp"
#endif
