#pragma once
//
//  baseContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 3/7/14.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
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
 virtual ~baseContainer(){}
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "baseContainer.cpp"
#endif
