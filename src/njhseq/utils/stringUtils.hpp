#pragma once
//
//  stringUtils.hpp
//
//  Created by Nick Hathaway on 12/24/12.
//
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
#include "njhseq/common/typedefs.hpp"
#include "njhseq/common/allSystemIncludes.h"
#include "utils.hpp"

namespace njhseq {




std::string getStringFromSubstrings(const std::string& seq,
                                    const std::vector<size_t>& positons,
                                    int subStringSize);
void trimStringAtFirstOccurence(std::string& str, const std::string& occurence);

void trimStringsAtFirstOccurence(VecStr& strings, const std::string& occurence);




// remove lower case letters and their corresponding qualities from the
// sequence.

VecStr tokenizeString(const std::string& str, const std::string& delim,
                      bool addEmptyToEnd = false);

std::vector<size_t> findOccurences(const std::string& target,
                                   const std::string& subSeq);

uint32_t countOccurences(const std::string& target, const std::string& subSeq);

uint32_t countRegexOccurrences(const std::string& target, const std::regex & subPat);


void translateStringWithKey(std::string& str, MapStrStr& key);


struct stringSorter {
  void sortStrByLength(VecStr& vec) {
    std::sort(vec.begin(), vec.end(),
              [](const std::string& a, const std::string& b) {
      if (a.size() == b.size()) {
        return a < b;
      }
      return a.size() < b.size();
    });
  }
};

std::string condenseSeqSimple(const std::string& seq);

/////lower case handling functions
std::vector<char> getUpperCaseLetters();
std::vector<char> getLowerCaseLetters();
std::vector<char> getLowerUpperCaseLetters();
std::vector<char> processAlphStrVecChar(const std::string & alphabetStr,
		const std::string & delim);
std::pair<std::vector<char>, std::vector<uint32_t>> processAlphStrVecCharCounts(const std::string & alphabetStr,
		const std::string & delim);
VecStr processAlphStrVecStr(const std::string & alphabetStr, const std::string & delim);
std::vector<char> determineAlph(const std::string & str);

void stringToUpper(std::string& str);
void stringToLower(std::string& str);
std::string stringToUpperReturn(std::string str);
std::string stringToLowerReturn(std::string str);
void changeSubStrToLower(std::string& str, size_t pos, size_t length);
void changeSubStrToLowerToEnd(std::string& str, size_t pos);
void changeSubStrToLowerFromBegining(std::string& str, size_t pos);
void changeSubStrToUpperFromBegining(std::string& str, size_t pos);
void subStrToUpper(std::string & str, uint32_t pos, uint32_t len);
void changeCertainSubStrToLower(std::string& str, const std::string& substring);
void changeStringVectorToLowerCase(VecStr& vec);

bool isIntStr(const std::string& str);
bool isDoubleStr(const std::string& str);
bool isVecOfIntStr(const VecStr& vec);
bool isVecOfDoubleStr(const VecStr& vec);

// conversion between strings and vectors
template <typename T>
std::string vectorToString(const std::vector<T>& vectorToConvert,
                           const std::string& delim = " ") {
  if (vectorToConvert.empty()) {
    return "";
  }
  std::stringstream tempStringStream;
  copy(vectorToConvert.begin(), vectorToConvert.end(),
       std::ostream_iterator<T>(tempStringStream, delim.c_str()));
  std::string returnString = tempStringStream.str().c_str();
  returnString.erase(returnString.size() - (int)delim.size());
  return returnString;
}
template<typename T>
std::string frameVec(const std::vector<T> & vec, const std::string & delim,
		const std::string & frontEnd, const std::string & backEnd){
	return frontEnd + vectorToString(vec, delim) + backEnd;
}



//trimming whitespace
std::size_t findFirstWhitespace(const std::string & str);
bool trimAtFirstWhitespace(std::string & str);
void trimEndWhiteSpace(std::string& str);
std::string trimEndWhiteSpaceReturn(std::string str);
bool allWhiteSpaceStr(const std::string & str);
//trim quotes
std::string stripQuotes( const std::string& str );



template<typename T, typename TOT>
std::string getPercentageString(T partial, TOT total) {
	if (total == 0) {
		return "0";
	}
	std::stringstream ss;
	ss << partial << "(" << std::setprecision(3)
			<< 100 * partial / static_cast<double>(total) << "%)";
	return ss.str();
}


// repeat string
std::string repeatString(const std::string& stringToRepeat, uint32_t n);

template <class T>
void kmerComposition(const T& read, std::map<std::string, int>& kMerMap,
                     int kLength) {
  for (auto i : iter::range(read.seqBase_.seq_.size() - kLength + 1)) {
    kMerMap[read.seqBase_.seq_.substr(i, kLength)] += read.seqBase_.cnt_;
  }
}

std::string removeCharReturn(std::string inputStr, const char& theChar);
void removeChar(std::string& inputStr, const char& theChar);

VecStr checkTwoRotatingStrings(const std::string& str1, const std::string& str2,
                               int allowableMismatches);
int numberOfMismatches(const std::string& str1, const std::string& str2);





void addToStr(std::string & str, char c);
void addToStr(std::string & str, const std::string & otherStr);


VecStr streamToVecStr(std::istream& stream);
template <typename T>
std::vector<T> stringToVector(const std::string& strToConvert) {
  std::vector<T> ans;
  std::stringstream ss;
  ss.clear();
  ss << trimEndWhiteSpaceReturn(strToConvert);
  while (!ss.eof()) {
    T a;
    ss >> a;
    ans.push_back(a);
  }
  return ans;
}


struct stringHasher {
	//http://stackoverflow.com/questions/7968674/unexpected-collision-with-stdhash
	size_t operator()( std::string const& s ) const{
		size_t result = 2166136261U ;
		for ( const auto & c : s ) {
				result = 127 * result + static_cast< unsigned char >( c ) ;
		}
		return result ;
	}
};

/**@b comparer for case insensitivity comparison, can be used in maps
 *
 */
struct strICaseCmp{
	struct charICaseCmp{
		bool operator()(const unsigned char& c1, const unsigned char& c2) const {
			return tolower(c1) < tolower(c2);
		}
	};
	bool operator()(const std::string & str1, const std::string & str2) const {
		return std::lexicographical_compare(str1.begin(), str1.end(), str2.begin(),
				str2.end(), charICaseCmp());
	}
};

uint32_t countBeginChar(const std::string & str);

uint32_t countEndChar(const std::string & str);
}  // namespace njhseq


