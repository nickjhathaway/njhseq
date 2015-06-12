
//
//  bestDistGraph.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/23/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bestDistGraph.hpp"
#include "bibseq/seqToolsUtils/seqToolsUtils.hpp"
#include <bibcpp/graphics/colorUtils.hpp>


namespace bibseq {

void bestDistGraph::printOutBestMatchInfos(std::ostream& out) {
  std::unordered_map<std::string, bool> alreadyAdded;
  for (const auto& n : nodes_) {
    if (!n.children_.empty()) {
      for (const auto& child : n.children_.begin()->second) {
        if (n.read_.getStubName(true) <
            nodes_[child.childNodePos_].read_.getStubName(true)) {
          if (contains(alreadyAdded, n.read_.getStubName(true) +
                                         nodes_[child.childNodePos_]
                                             .read_.getStubName(true))) {
            continue;
          }
          out << nodes_[child.childNodePos_].read_.seqBase_.name_ << "\t"
              << n.read_.seqBase_.name_ << "\t"
              << nodes_[child.childNodePos_].read_.seqBase_.frac_ << "\t"
              << n.read_.seqBase_.frac_;
          alreadyAdded[n.read_.getStubName(true) +
                       nodes_[child.childNodePos_].read_.getStubName(true)] =
              true;
        } else {
          if (contains(alreadyAdded,
                       nodes_[child.childNodePos_].read_.getStubName(true) +
                           n.read_.getStubName(true))) {
            continue;
          }
          out << n.read_.seqBase_.name_ << "\t"
              << nodes_[child.childNodePos_].read_.seqBase_.name_ << "\t"
              << n.read_.seqBase_.frac_ << "\t"
              << nodes_[child.childNodePos_].read_.seqBase_.frac_;
          alreadyAdded[nodes_[child.childNodePos_].read_.getStubName(true) +
                       n.read_.getStubName(true)] = true;
        }
        out << "\t" << child.distance_ << "\t" << child.gapDistance_ << "\t"
            << child.misMatches_ << "\t" << child.identity_ << "\t"
            << convertBoolToString(child.differsByIndels_) << "\t"
            << convertBoolToString(child.differsByMismatches_);
        out << std::endl;
      }
    }
  }
}

void bestDistGraph::printOutAllInfos(std::ostream& out) {
	out << "childName\tparentName\tchildFraq\tparentFraq\t"
			"dis\tgapDis\tmisMatches\tidentity\tdiffbyIndels\tdiffByMismatches\tadded\tgaps\tgappedBases"<< std::endl;
  for (const auto& n : nodes_) {
    for (const auto& childVec : n.children_) {
      for (const auto& child : childVec.second) {
        out << nodes_[child.childNodePos_].read_.seqBase_.name_ << "\t"
            << n.read_.seqBase_.name_ << "\t"
            << nodes_[child.childNodePos_].read_.seqBase_.frac_ << "\t"
            << n.read_.seqBase_.frac_ << "\t" << child.distance_ << "\t"
            << child.gapDistance_ << "\t" << child.misMatches_ << "\t"
            << child.identity_ << "\t"
            << convertBoolToString(child.differsByIndels_) << "\t"
            << convertBoolToString(child.differsByMismatches_) << "\t"
            << convertBoolToString(nodes_[child.childNodePos_].added_) << "\t"
        		<< child.numGaps_ << "\t" << child.numGappedBases_;
        out << std::endl;
      }
    }
  }
}
void bestDistGraph::printOutGV(std::ostream& out) {
  std::unordered_map<std::string, bool> alreadyAdded;
  out << "graph G  { " << std::endl;
  out << "bgcolor =\"#000000\"" << std::endl;
  out << "overlap = false; " << std::endl;
  out << "fixedsize = true; " << std::endl;
  out << "fontcolor = \"#ffffff\"" << std::endl;
  out << "fontsize = 20" << std::endl;
  // out << "node [shape=circle,style=filled,width=.3, height=.3]" << std::endl;
  uint32_t count = 0;
  std::vector<bib::color> colors = bib::evenHuesAll(0.75, 0.45, nodes_.size());
  for (const auto& n : nodes_) {
    out << n.read_.getReadId() << " [shape=circle,style=filled,fixesize =true, "
                                  "color = \"#000000\", fillcolor ="
        << "\"#" << colors[count].hexStr_ << "\""
        << ", width = " << n.read_.seqBase_.frac_ * 50 << "]" << std::endl;
    ++count;
  }
  // uint32_t maxDist = vectorMaximum(getVectorOfMapKeys(bestDists_)) + 1;
  // std::cout <<"maxDist: " << maxDist << std::endl;
  for (const auto& n : nodes_) {
    if (!n.children_.empty()) {
      for (const auto& child : n.children_.begin()->second) {
        if (n.read_.getStubName(true) <
            nodes_[child.childNodePos_].read_.getStubName(true)) {
          if (contains(alreadyAdded, n.read_.getStubName(true) +
                                         nodes_[child.childNodePos_]
                                             .read_.getStubName(true))) {
            continue;
          }
          alreadyAdded[n.read_.getStubName(true) +
                       nodes_[child.childNodePos_].read_.getStubName(true)] =
              true;
        } else {
          if (contains(alreadyAdded,
                       nodes_[child.childNodePos_].read_.getStubName(true) +
                           n.read_.getStubName(true))) {
            continue;
          }
          alreadyAdded[nodes_[child.childNodePos_].read_.getStubName(true) +
                       n.read_.getStubName(true)] = true;
        }
        for (const auto& i : iter::range<uint32_t>(0, child.misMatches_)) {
          out << nodes_[child.childNodePos_].read_.getReadId() << i
              << n.read_.getReadId()
              << " [shape=point,style=filled,color = \"#000000\", fillcolor ="
              << "\"#F8766D\""
              << ", width = " << 0.1 << " , label =\"\"]" << std::endl;
        }
        std::string lastNode = "";
        if (child.misMatches_ != 0) {
          lastNode = nodes_[child.childNodePos_].read_.getReadId() +
                     to_string(0) + n.read_.getReadId();
          out << nodes_[child.childNodePos_].read_.getReadId() << " -- "
              << lastNode << " " << std::endl;
          for (const auto& i : iter::range<uint32_t>(1, child.misMatches_)) {
            out << lastNode << " -- "
                << nodes_[child.childNodePos_].read_.getReadId() << i
                << n.read_.getReadId() << std::endl;
            lastNode = nodes_[child.childNodePos_].read_.getReadId() +
                       to_string(i) + n.read_.getReadId();
          }
          out << lastNode << " -- " << n.read_.getReadId() << " " << std::endl;
        } else {
          out << nodes_[child.childNodePos_].read_.getReadId() << " -- "
              << n.read_.getReadId() << " " << std::endl;
        }
      }
    }
  }
  out << "}" << std::endl;
}
void bestDistGraph::printOutBestMatch(
    std::ostream& out, std::ostream& gvOut, uint32_t nodePos,
    std::unordered_map<std::string, bool>& alreadyAdded, bool best,
    uint32_t misDist, bool recursive,
    const std::unordered_map<std::string, bib::color>& colsForName) {
  auto& n = nodes_[nodePos];
  if (!n.children_.empty()) {
    auto childrenIter = n.children_.begin();
    if (!best) {
      childrenIter = n.children_.find(misDist);
      if (childrenIter == n.children_.end()) {
        return;
      }
    }
    for (const auto& child : childrenIter->second) {
      if (n.read_.getStubName(true) <
          nodes_[child.childNodePos_].read_.getStubName(true)) {
        if (contains(alreadyAdded,
                     n.read_.getStubName(true) +
                         nodes_[child.childNodePos_].read_.getStubName(true))) {
          continue;
        }
        out << nodes_[child.childNodePos_].read_.seqBase_.name_ << "\t"
            << n.read_.seqBase_.name_ << "\t"
            << nodes_[child.childNodePos_].read_.seqBase_.frac_ << "\t"
            << n.read_.seqBase_.frac_;
        alreadyAdded[n.read_.getStubName(true) +
                     nodes_[child.childNodePos_].read_.getStubName(true)] =
            true;
      } else {
        if (contains(alreadyAdded,
                     nodes_[child.childNodePos_].read_.getStubName(true) +
                         n.read_.getStubName(true))) {
          continue;
        }
        out << n.read_.seqBase_.name_ << "\t"
            << nodes_[child.childNodePos_].read_.seqBase_.name_ << "\t"
            << n.read_.seqBase_.frac_ << "\t"
            << nodes_[child.childNodePos_].read_.seqBase_.frac_;
        alreadyAdded[nodes_[child.childNodePos_].read_.getStubName(true) +
                     n.read_.getStubName(true)] = true;
      }
      out << "\t" << child.distance_ << "\t" << child.gapDistance_ << "\t"
          << child.misMatches_ << "\t" << child.identity_ << "\t"
          << convertBoolToString(child.differsByIndels_) << "\t"
          << convertBoolToString(child.differsByMismatches_);
      out << std::endl;
      // std::cout <<"b1" << std::endl;

      std::vector<bib::color> misLineCols = bib::evenHuesInbetweenTwo(
          colsForName.at(nodes_[child.childNodePos_].read_.getReadId()),
          colsForName.at(n.read_.getReadId()), child.misMatches_ + 1);
      // std::cout << "b2" << std::endl;
      // std::cout << "size: " << misLineCols.size() << std::endl;
      // colsForName.at(nodes_[child.childNodePos_].read_.getReadId()).printDescription(std::cout
      // , false);
      // colsForName.at(n.read_.getReadId()).printDescription(std::cout ,
      // false);
      /*uint32_t colCount = 0;
      for (const auto& c : misLineCols) {
        std::cout << "colCount: " << colCount << std::endl;
        c.printDescription(std::cout, false);
        ++colCount;
      }*/
      for (const auto& i : iter::range<uint32_t>(0, child.misMatches_)) {
        gvOut << nodes_[child.childNodePos_].read_.getReadId() << i
              << n.read_.getReadId()
              << " [shape=point,style=filled,color = \"#000000\", fillcolor ="
              << "\"red\""
              << ", width = " << 0.15 << " , label =\"\"]" << std::endl;
      }
      std::string lastNode = "";
      if (child.misMatches_ != 0) {

        lastNode = nodes_[child.childNodePos_].read_.getReadId() +
                   to_string(0) + n.read_.getReadId();
        gvOut << nodes_[child.childNodePos_].read_.getReadId() << " -- "
              << lastNode << " [penwidth=5, color=\"#" << misLineCols[0].hexStr_
              << "\"]" << std::endl;
        for (const auto& i : iter::range<uint32_t>(1, child.misMatches_)) {
          gvOut << lastNode << " -- "
                << nodes_[child.childNodePos_].read_.getReadId() << i
                << n.read_.getReadId() << " [penwidth=5, color=\"#"
                << misLineCols[i].hexStr_ << "\"]" << std::endl;
          lastNode = nodes_[child.childNodePos_].read_.getReadId() +
                     to_string(i) + n.read_.getReadId();
        }
        gvOut << lastNode << " -- " << n.read_.getReadId()
              << " [penwidth=5, color=\"#"
              << misLineCols[child.misMatches_].hexStr_ << "\"]" << std::endl;
      } else {
        gvOut << nodes_[child.childNodePos_].read_.getReadId() << " -- "
              << n.read_.getReadId() << " [penwidth=5, color=\"#"
              << misLineCols[0].hexStr_ << "\"]" << std::endl;
      }
      if (contains(mainClusterNames_, n.read_.getStubName(true)) &&
          !contains(mainClusterNames_,
                    nodes_[child.childNodePos_].read_.getStubName(true))) {
        mainClusterNames_[n.read_.getStubName(true)]
            .emplace_back(nodes_[child.childNodePos_].read_.getStubName(true));
        mainClusterNames_[nodes_[child.childNodePos_].read_.getStubName(true)]
            .emplace_back(n.read_.getStubName(true));
        nodes_[child.childNodePos_].added_ = true;
      }
      if (contains(mainClusterNames_,
                   nodes_[child.childNodePos_].read_.getStubName(true)) &&
          !contains(mainClusterNames_, n.read_.getStubName(true))) {
        mainClusterNames_[n.read_.getStubName(true)]
            .emplace_back(nodes_[child.childNodePos_].read_.getStubName(true));
        mainClusterNames_[nodes_[child.childNodePos_].read_.getStubName(true)]
            .emplace_back(n.read_.getStubName(true));
        n.added_ = true;
      }
      if (recursive) {
        printOutBestMatch(out, gvOut, child.childNodePos_, alreadyAdded, best,
                          misDist, recursive, colsForName);
      }
    }
  }
}
bool bestDistGraph::allHaveBeenAdded() {
  for (const auto& n : nodes_) {
    if (!n.added_) {
      return false;
    }
  }
  return true;
}
void bestDistGraph::createDotBestConnectedFile(
    const std::string& workingDir,
    const std::unordered_map<std::string, bib::color>& colorsForName,
    bool printTextFile) {
  auto colorsToUse = colorsForName;
  if (colorsToUse.empty()) {
    uint32_t count = 0;
    std::vector<bib::color> colors = bib::evenHuesAll(0.75, 0.45, nodes_.size());
    for (const auto& n : nodes_) {
      colorsToUse[n.read_.getReadId()] = colors[count];
      ++count;
    }
  }
  std::ofstream gvOut;
  openTextFile(gvOut, workingDir + otuName_, ".dot", true, false);
  std::stringstream holdThis;
  std::ofstream txtOut;
  if (printTextFile) {
    openTextFile(txtOut, workingDir + otuName_, ".txt", true, false);
  }
  gvOut << "graph G  { " << std::endl;
  gvOut << "bgcolor =\"#000000\"" << std::endl;
  gvOut << "labelloc=\"t\"" << std::endl;
  gvOut << "fontcolor = \"#ffffff\"" << std::endl;
  gvOut << "fontsize = 20" << std::endl;
  gvOut << "label = \"" << otuName_ << "\"" << std::endl;
  gvOut << "fixedsize = true; " << std::endl;
  // uint32_t count = 0;
  // uint32_t colStep = colors.size()/(testDist.nodes_.size()+1);
  for (const auto& n : nodes_) {
    gvOut << n.read_.getReadId() << " [shape=circle,style=filled,fixesize "
                                    "=true, color = \"#000000\", fillcolor ="
          << "\"#" << colorsToUse.at(n.read_.getReadId()).hexStr_ << "\""
          << ", width = " << sqrt(n.read_.seqBase_.frac_ * 50) << "]"
          << std::endl;
    //++count;
  }
  // std::cout <<"mark 1 " << std::endl;
 // std::cout << "nodesSize: " << nodes_.size() << std::endl;
 // std::cout << "parent node: " << parentNode_ << std::endl;

  std::unordered_map<std::string, bool> alreadyAdded;
  nodes_[parentNode_].added_ = true;
  // std::cout <<"mark 2 " << std::endl;
  //std::cout << parentNode_ << std::endl;
  //std::cout << nodes_[parentNode_].read_.getStubName(true) << std::endl;
  mainClusterNames_[nodes_[parentNode_].read_.getStubName(true)] = {};
  // std::cout <<"mark 3 " << std::endl;
  if (printTextFile) {
    printOutBestMatch(txtOut, gvOut, parentNode_, alreadyAdded, true, 1, true,
                      colorsForName);
  } else {
    // std::cout <<"mark 4 " << std::endl;
    printOutBestMatch(holdThis, gvOut, parentNode_, alreadyAdded, true, 1, true,
                      colorsForName);
    // std::cout <<"mark 5 " << std::endl;
  }
  // std::cout <<"mark 6 " << std::endl;
  auto distIter = bestDists_.begin();
  while (!allHaveBeenAdded() && distIter != bestDists_.end()) {
    for (const auto& nodePos : iter::range(nodes_.size())) {
      if (nodes_[nodePos].added_) {
        continue;
      }
      if (printTextFile) {
        printOutBestMatch(txtOut, gvOut, nodePos, alreadyAdded, false,
                          distIter->first, false, colorsForName);
      } else {
        printOutBestMatch(holdThis, gvOut, nodePos, alreadyAdded, false,
                          distIter->first, false, colorsForName);
      }
    }
    // std::cout << distIter->first << std::endl;
    ++distIter;
  }
  gvOut << "}" << std::endl;
}

}  // bib
