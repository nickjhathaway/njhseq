
//
//  otuGraph.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/15/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "otuGraph.hpp"

namespace bibseq {

void otuGraph::printParent(std::ostream& out) {
  out << otuName_ << "\t" << nodes_[parentPos_].read_.seqBase_.name_ << "\t"
      << nodes_[parentPos_].read_.seqBase_.cnt_ << "\tparent\t0\t0"
      << std::endl;
}

void otuGraph::printChildren(std::ostream& out, const node& currentNode) {
  for (const auto& child : currentNode.children_) {
    out << otuName_ << "\t" << nodes_[child.childNodePos_].read_.seqBase_.name_
        << "\t" << nodes_[child.childNodePos_].read_.seqBase_.cnt_ << "\t"
        << currentNode.read_.seqBase_.name_ << "\t"
        << currentNode.read_.seqBase_.cnt_ << "\t" << child.distance_
        << std::endl;
    printChildren(out, nodes_[child.childNodePos_]);
  }
}
void otuGraph::printChildrenGV(std::ostream& out, const node& currentNode) {
  for (const auto& child : currentNode.children_) {
    out << currentNode.read_.getStubName(true) << " -- "
        << nodes_[child.childNodePos_].read_.getStubName(true)
        << " [len=" << child.distance_ << "];" << std::endl;
    printChildrenGV(out, nodes_[child.childNodePos_]);
  }
}
void otuGraph::printInfo(std::ostream& out) {
  printParent(out);
  printChildren(out, nodes_[parentPos_]);
}
void otuGraph::printInfoGraphViz(std::ostream& out) {
  out << "graph G  {" << std::endl;
  out << "overlap = false" << std::endl;
  out << "node [shape=circle,style=filled,width=.3, height=.3];" << std::endl;
  printChildrenGV(out, nodes_[parentPos_]);
  out << "}" << std::endl;
}
void otuGraph::printParentLevel(std::ostream& out) {
  double two_pi = 2 * (atan(1) * 4);
  out << "parent\t50\t50\t0\t" << 50 * nodes_[parentPos_].read_.seqBase_.frac_;
  out << "\t0\t0" << std::endl;
  double radsStep = two_pi / nodes_[parentPos_].children_.size();
  uint32_t childNum = 0;
  for (const auto& child : nodes_[parentPos_].children_) {
    out << "child\t50\t50\t";
    out << childNum* radsStep << "\t"
        << 50 * nodes_[child.childNodePos_].read_.seqBase_.frac_ << "\t"
        << child.distance_ << "\t"
        << 50 * nodes_[parentPos_].read_.seqBase_.frac_ << std::endl;
    printChildrenLevel(
        out, 50 + cos(childNum * radsStep) *
                      (50 * nodes_[parentPos_].read_.seqBase_.frac_ +
                       child.distance_ +
                       50 * nodes_[child.childNodePos_].read_.seqBase_.frac_),
        50 + sin(childNum * radsStep) *
                 (50 * nodes_[parentPos_].read_.seqBase_.frac_ +
                  child.distance_ +
                  50 * nodes_[child.childNodePos_].read_.seqBase_.frac_),
        childNum * radsStep, child.childNodePos_);
    ++childNum;
  }
}
void otuGraph::printChildrenLevel(std::ostream& out, double x1, double y1,
                                  double currentRadians,
                                  uint32_t currentNodePos) {
  if (nodes_[currentNodePos].children_.size() == 0) {
    return;
  }
  double pi = (atan(1) * 4);
  double radsStep = pi / nodes_[currentNodePos].children_.size();
  uint32_t childNum = 0;
  for (const auto& child : nodes_[currentNodePos].children_) {
    out << "child\t" << x1 << "\t" << y1 << "\t";
    double rads = 0;
    if (childNum % 2 == 0) {
      rads = currentRadians + childNum * radsStep;
    } else {
      rads = currentRadians - childNum * radsStep;
    }
    out << rads << "\t";
    out << 50 * nodes_[child.childNodePos_].read_.seqBase_.frac_ << "\t"
        << child.distance_ << "\t"
        << 50 * nodes_[currentNodePos].read_.seqBase_.frac_ << std::endl;
    printChildrenLevel(
        out,
        x1 + cos(rads) * (50 * nodes_[currentNodePos].read_.seqBase_.frac_ +
                          child.distance_ + 50 * nodes_[child.childNodePos_]
                                                     .read_.seqBase_.frac_),
        y1 + sin(rads) * (50 * nodes_[currentNodePos].read_.seqBase_.frac_ +
                          child.distance_ + 50 * nodes_[child.childNodePos_]
                                                     .read_.seqBase_.frac_),
        rads, child.childNodePos_);

    ++childNum;
  }
}

std::string otuGraph::appendChildName(std::string& name,
                                      uint32_t currentNodePos) {
  for (const auto& child : nodes_[currentNodePos].children_) {
    name.append("_" + nodes_[child.childNodePos_].read_.seqBase_.name_);
    appendChildName(name, child.childNodePos_);
  }
  return name;
}
void otuGraph::getNames(VecStr& names, uint32_t nodePos,
                        std::string currentName) {
  if (nodes_[nodePos].children_.empty()) {
    currentName.append(nodes_[nodePos].read_.seqBase_.name_);
    names.emplace_back(currentName);
  } else {
    currentName.append(nodes_[nodePos].read_.seqBase_.name_ + "->");
    for (const auto& child : nodes_[nodePos].children_) {
      getNames(names, child.childNodePos_, currentName);
    }
  }
}

void otuGraph::printSimpleChildrenInfo(std::ostream& out, uint32_t nodePos) {
  for (const auto& child : nodes_[nodePos].children_) {
    out << nodes_[child.childNodePos_].read_.seqBase_.name_ << "\t"
        << nodes_[nodePos].read_.seqBase_.name_ << "\t"
        << nodes_[child.childNodePos_].read_.seqBase_.frac_ << "\t"
        << nodes_[nodePos].read_.seqBase_.frac_ << "\t" << child.distance_
        << "\t" << child.gapDistance_ << "\t" << child.misMatches_ << "\t"
        << child.identity_ << "\t"
        << convertBoolToString(child.differsByIndels_) << "\t"
        << convertBoolToString(child.differsByMismatches_);
    out << std::endl;
    printSimpleChildrenInfo(out, child.childNodePos_);
  }
}
}