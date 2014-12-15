#pragma once
//
//  stdAddition.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/8/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "allSystemIncludes.h"
#include <bibcpp.h>

namespace bibseq {

using bib::count_if;
using bib::for_each;
using bib::sort;
using bib::reverse;
using bib::contains;
using bib::find;
using bib::has;
using bib::in;
using bib::iota;
using bib::shuffle;
using estd::to_string;


template <typename T>
const uint32_t len(const T& c) {
  return static_cast<uint32_t>(c.size());
}



}  // bib



