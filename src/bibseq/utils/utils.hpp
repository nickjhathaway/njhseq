#pragma once
//
//  utils.hpp
//  ampliconCluster
//
//  Created by Nick Hathaway on 8/31/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//

#include "bibseq/common.h"
#include "bibseq/utils/bitSwaps.hpp"
//////durations and date
namespace bibseq {

// will call printDescription on a vector of objects, most of the objects in seqTools have
// such a function

template<typename T>
void printDescriptionVec(const std::vector<T> & vec, std::ostream & out, bool deep){
	for_each(vec, [&](const T & obj){ obj.printDescription(out, deep);});
}

template<typename T>
std::map<T, uint32_t> countVec(const std::vector<T> & vec){
	std::map<T, uint32_t> counts;
	for(const auto & element : vec){
		++counts[element];
	}
	return counts;
}

std::map<uint64_t, uint32_t> printStringLengths(const VecStr& strings,
                                                std::ostream& out = std::cout);

template <typename T>
std::string leftPadNumStr(T num, T highestNumber = 10) {
  std::stringstream ans;
  if (num == 0) {
    ans << std::string(log10(highestNumber), '0');
  } else {
    ans << std::string(((int)log10(highestNumber) - (int)log10(num)), '0');
  }
  ans << num;
  return ans.str();
}

std::string getCurrentDate();

std::string convertBoolToString(bool convert);



/**@b Print out the contents of a map
 *
 * @param theMap The map to print
 * @param delim The delimiter to use
 * @param out The stream to print o
 */
template <typename MAP>
void printOutMapContents(const MAP& theMap,
                         const std::string& delim,
                         std::ostream& out) {
  for (const auto& mValue : theMap) {
    out << mValue.first << delim << mValue.second << "\n";
  }
}


/**@b Get a vector of all the values stored in theMap
 *
 * @param theMap The map to get the values from
 * @return A vector of the mapped values
 */
template <typename MAP>
static std::vector<typename MAP::mapped_type> getVectorOfMapValues(const MAP& theMap) {
  std::vector<typename MAP::mapped_type> ret;
  for (const auto & mValue : theMap) {
    ret.push_back(mValue.second);
  }
  return ret;
}

/**@b Get a vector of all the keys stored in theMap
 *
 * @param theMap The map to get the keys from
 * @return A vector of the key values
 */
template <typename MAP>
static std::vector<typename MAP::key_type> getVectorOfMapKeys(const MAP& theMap) {
  std::vector<typename MAP::key_type> ret;
  for (const auto & mValue : theMap) {
    ret.push_back(mValue.first);
  }
  return ret;
}



// with no header
void printTableOrganized(const std::vector<VecStr>& content, std::ostream& out);
// with header
void printTableOrganized(const std::vector<VecStr>& content,
                         const VecStr& header, std::ostream& out);

template <typename FIR, typename SEC>
VecStr pairToVecStr(const std::pair<FIR, SEC>& p) {
  return VecStr{to_string(p.first), to_string(p.second)};
}


size_t WriteCallback(char* contents, size_t size, size_t nmemb,
                     std::ostream* stream);

std::string GetURL(const std::string url);
void GetURLStream(const std::string url, std::ostream & out);

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || \
    (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif(defined(_AIX) || defined(__TOS__AIX__)) || \
    (defined(__sun__) || defined(__sun) ||       \
     defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS();

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS();
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "utils.cpp"
#endif
