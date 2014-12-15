#pragma once
//
//  Created by Michael Purcaro on 11/7/13.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//

//#error Precompiled header file not found



#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <math.h>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <unordered_map>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <iterator>
#include <cstddef>

#include <functional>

// Some files need this, even if DEBUG or PROG_DEBUG is not defined.
#include <cassert>

// The below is needed for run-time type information (RTTI)
// where there is a typeid.
#include <typeinfo>

#include <cstdio>
#include <stdlib.h>
#include <cstdlib>

// Need below to define size_t
#include <cstddef>
// Need below to define numeric_limits<TYPE_HERE>::min(),::max(), etc.
#include <limits>
#include <numeric>
#include <math.h>

#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <list>
#include <stack>
#include <queue>
#include <deque>
#include <set>
#include <algorithm>
#include <iterator>
// Below defines #pair
#include <utility>

// For toupper()
#include <cctype>

// For isdigit()
#include <ctype.h>

// Below is for file modification times
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <ctime>

#include <sys/types.h>
#include <regex.h>
//
//
#include <openssl/md5.h>
//

#include <exception>
#include <stdexcept>
#include <memory>
#include <thread>
#include <mutex>
#include <curl/curl.h>

#include <cppitertools/range.hpp>
#include <cppitertools/reverse.hpp>
#include <cppitertools/enumerate.hpp>
#include "ediannessmacros.h"
