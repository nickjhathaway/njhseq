#pragma once
//
//  baseContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
namespace bibseq {

template <typename T>
class baseContainer {

 public:
  // contructors
  baseContainer() : seqBase_(seqInfo()) {}
  baseContainer(const seqInfo& seqBase) : seqBase_(seqBase) {}
  baseContainer(const seqInfo& seqBase, const std::vector<T>& reads)
      : seqBase_(seqBase), reads_(reads) {}
  // members
  seqInfo seqBase_;
  std::vector<T> reads_;
  // functions
  // ading reads
  virtual void addRead(const T& read) {
    reads_.emplace_back(read);
    seqBase_.cnt_ += read.seqBase_.cnt_;
    seqBase_.frac_ += read.seqBase_.frac_;
  }
  virtual void writeReads(std::ofstream & outFile, bool writeBase)const{
  	if(writeBase){
  		seqBase_.outPutFastq(outFile);
  	}
  	for(const auto & read : reads_){
  		read.seqBase_.outPutFastq(outFile);
  	}
  }

  // description
  virtual void printDescription(std::ostream& out, bool deep) const {
    out << "baseContainer{" << std::endl;
    out << "seqBase_:" << std::endl;
    seqBase_.printDescription(out, deep);
    if (deep) {
      out << "reads:{" << std::endl;
      out << "std::vector<T> " << std::endl;
      for (const auto& read : reads_) {
        read.printDescription(out, deep);
      }
    }
    out << "}" << std::endl;
  }
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "baseContainer.cpp"
#endif
