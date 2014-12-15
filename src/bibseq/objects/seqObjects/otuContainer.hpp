#pragma once
//
//  otuContainer.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 02/28/14.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/baseReadObject.hpp"

namespace bibseq {
/*
template<typename T>
class otuContainer {
public:

        //contructors
        otuContainer(const seqInfo & seqBase): seqBase_(seqBase){}
        otuContainer(const T & read): seqBase_(read.seqBase_), reads_({read}){}

        //members
        seqInfo seqBase_;
        std::vector<T> reads_;

        //functions
        virtual void addRead(const T & read){
                reads_.emaplce_back(read);
                seqBase_.cnt_+= read.seqBase_.cnt_;
                seqBase_.frac_+= read.seqBase_.frac_;
        }


private:
};
*/
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "otuContainer.cpp"
#endif
