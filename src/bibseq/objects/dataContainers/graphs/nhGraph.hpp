#pragma once
//
//  nhGraph.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 12/5/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
#include "bibseq/utils.h"
#include "bibseq/alignment.h"
namespace bibseq {
// template<typname T>
class nhGraph {

 public:
  nhGraph(const std::vector<readObject>& reads, aligner& alignerObj,
          int minOverLap, int maxOverLap);
  struct readFixes {
    readFixes(const readObject& read, int minOverLap, int maxOverLap)
        : _read(read) {
      for (auto i : iter::range(minOverLap, maxOverLap + 1)) {
        _prefixes[i] = _read.seqBase_.seq_.substr(0, i);
        _suffixes[i] =
            _read.seqBase_.seq_.substr(read.seqBase_.seq_.length() - i);
      }
    }
    readObject _read;
    std::unordered_map<int, std::string> _prefixes;
    std::unordered_map<int, std::string> _suffixes;
  };
  struct node {
    node(const readObject& read, int minOverLap, int maxOverLap)
        : _read(read), _readEnds(readFixes(read, minOverLap, maxOverLap)) {}
    readObject _read;
    readFixes _readEnds;
    // forward, for this graph head overlap
    std::vector<std::pair<uint, std::string>> _headConnections;
    // backward, for this graph tail overlap
    std::vector<std::pair<uint, std::string>> _tailConnections;
  };
  struct connection {
    connection(const std::string& value, const std::pair<uint, uint>& firstCon)
        : _value(value), _nodePositions({firstCon}) {}
    std::string _value;
    // first being the tail connection possiton and the second being the head
    // connection position
    std::unordered_map<uint, uint> _nodePositions;
  };
  std::vector<node> _nodes;
  void addNode(const readObject& read, aligner& alignerObj, int minOverLap,
               int maxOverLap);
  void printAdjacencyTable(std::ostream& out, bool printName);
  std::multimap<uint, uint> getAdjacencyMap();
};

}  // namespace bib
#ifndef NOT_HEADER_ONLY
#include "nhGraph.cpp"
#endif
