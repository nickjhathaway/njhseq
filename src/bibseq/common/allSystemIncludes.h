#pragma once
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
#include <openssl/md5.h>
#include <curl/curl.h>
#include <cppitertools/range.hpp>
#include <cppitertools/reverse.hpp>
#include <cppitertools/enumerate.hpp>

//own files
#include "ediannessmacros.h"
