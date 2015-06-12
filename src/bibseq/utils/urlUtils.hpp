#pragma once
//
//  urlUtils.hpp
//
//  Created by Nick Hathaway on 5/27/15.
//  Copyright (c) 2015 University of Massachusetts Medical School. All rights
// reserved.
//

#include "bibseq/common.h"
#include "bibseq/utils/bitSwaps.hpp"

namespace bibseq {

size_t WriteCallback(char* contents, size_t size, size_t nmemb,
                     std::ostream* stream);
std::string GetURL(const std::string url);
void GetURLStream(const std::string url, std::ostream & out);


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "urlUtils.cpp"
#endif
