#pragma once
//
//  otuContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/objects/seqContainers/baseContainer.hpp"

namespace bibseq {

template <typename T>
class otuContainer : public baseContainer<T> {

 public:
  // contructors
  otuContainer() : baseContainer<T>() {}
  otuContainer(const seqInfo& seqBase) : baseContainer<T>(seqBase) {}
  otuContainer(const seqInfo& seqBase, const std::vector<T>& reads)
      : baseContainer<T>(seqBase, reads) {}
  otuContainer(const seqInfo& seqBase, double percentIdentity,
               double gapPercentCutOff)
      : percentIdentity_(percentIdentity),
        gapPercentCutOff_(gapPercentCutOff),
        baseContainer<T>(seqBase) {}
  otuContainer(const seqInfo& seqBase, const std::vector<T>& reads,
               double percentIdentity, double gapPercentCutOff)
      : percentIdentity_(percentIdentity),
        gapPercentCutOff_(gapPercentCutOff),
        baseContainer<T>(seqBase, reads) {}
  // members
  double percentIdentity_ = 0.97;
  double gapPercentCutOff_ = 0.03;

  // functions
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "otuContainer.cpp"
#endif
