#pragma once
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
//
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//

//#error Precompiled header file not found

//C libraries
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <cerrno>
#include <cstddef>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
// Need below to define size_t
#include <cstddef>
#include <cctype>
#include <ctime>


//system files
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>


//containers
#include <map>
#include <vector>
#include <unordered_map>
#include <array>
#include <list>
#include <stack>
#include <queue>
#include <deque>
#include <set>
#include <tuple>
#include <utility> //defines std::pair

//strings
#include <sstream>
#include <string>

//io
#include <iostream> //defines cout, cerr, etc
#include <fstream>
#include <iomanip>

//misc std libraries
#include <algorithm>
#include <iterator>
#include <functional>
#include <limits> //to define numeric_limits<TYPE_HERE>::min(),::max(), etc.
#include <numeric>
#include <random>
#include <regex>
#include <memory> //smart pointers
//std concurrency
#include <thread>
#include <mutex>

//exception handling
#include <exception>
#include <stdexcept>
// The below is needed for run-time type information (RTTI)
// where there is a typeid.
#include <typeinfo>



//non standard libraries
//#include <openssl/md5.h>
#include <curl/curl.h>
#include <cppitertools/range.hpp>

#if __has_include("cppitertools/reverse.hpp")
#include <cppitertools/reverse.hpp>
#else
#include <cppitertools/reversed.hpp>
namespace iter {
/**@todo simple fix for now, should just change all occurrences of reverse to reversed*/
template <typename Container>
iter::impl::Reverser<Container> reverse(Container&& container) {
  return reversed(std::forward<Container>(container));
}
}  // namespace iter
#endif

#include <cppitertools/enumerate.hpp>

//own files
#include "ediannessmacros.h"
