#pragma once
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
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
#include <cppitertools/itertools.hpp>


//own files
#include "ediannessmacros.h"
#if __APPLE__ == 1 && __cpp_lib_shared_timed_mutex < 201402L && __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ <= 101106
#include <sharedMutex.h>
#else
#include <shared_mutex>
#endif


