#pragma once
//
//  refMapContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/objects/seqContainers/baseContainer.hpp"

namespace bibseq {

template <typename T>
class refMapContainer : public baseContainer<T> {

 public:
  // Constructors
  refMapContainer() : baseContainer<T>() {}
  refMapContainer(const seqInfo& seqBase) : baseContainer<T>(seqBase) {
    initialize();
  }
  refMapContainer(const seqInfo& seqBase, const std::vector<T>& reads)
      : baseContainer<T>(seqBase, reads) {
    initialize();
  }

  void initialize() {
    // no reads mapped yet so initialize counts to zero
    this->seqBase_.cnt_ = 0;
    this->seqBase_.frac_ = 0;
  }
  // members

  // functions
  template<typename READS, typename REFS>
  static std::vector<refMapContainer> createContainers(const std::vector<REFS> & reads){
  	std::vector<refMapContainer> refContainers;
  	for(const auto & read : reads){
  		refContainers.emplace_back(refMapContainer<READS>(read.seqBase_));
  	}
  	return refContainers;
  }
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "refMapContainer.cpp"
#endif
