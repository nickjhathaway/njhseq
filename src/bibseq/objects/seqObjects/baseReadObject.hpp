#pragma once
//
//  baseReadObject.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/31/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
namespace bibseq {

class baseReadObject {

 public:
  baseReadObject() : seqBase_(seqInfo()) {}
  baseReadObject(const seqInfo& seqBase) : seqBase_(seqBase) {}

  seqInfo seqBase_;

  // description
  virtual void printDescription(std::ostream& out, bool deep = false) const;
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "baseReadObject.cpp"
#endif
