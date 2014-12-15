#pragma once
//
//  stringUtils.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 12/24/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//

#include "bibseq/common/typedefs.hpp"
#include "bibseq/common/allSystemIncludes.h"
#include "utils.hpp"

namespace bibseq {
std::string get_cwd();


//finding substrings
bool containsSubString(const std::string& str, const std::string& subString);
bool endsWith(const std::string& a, const std::string& b);
bool beginsWith(const std::string& str, const std::string& target);


std::string getTimeFormat(double timeInSecondsOriginal, bool wordy,
                          int secondsDecimalPlaces = 6);

std::string getStringFromSubstrings(const std::string& seq,
                                    const std::vector<size_t>& positons,
                                    int subStringSize);
void trimStringAtFirstOccurence(std::string& str, const std::string& occurence);

void trimStringsAtFirstOccurence(VecStr& strings, const std::string& occurence);



std::string replaceString(std::string theString,
                          const std::string& replaceSpace = " ",
                          const std::string& newSpace = "_");

// remove lower case letters and their corresponding qualities from the
// sequence.

VecStr tokenizeString(const std::string& str, const std::string& delim,
                      bool addEmptyToEnd = false);

std::vector<size_t> findOccurences(const std::string& target,
                                   const std::string& subSeq);

int countOccurences(const std::string& target, const std::string& subSeq);

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

bool stringContainsAllDigits(const std::string& str);
bool stringContainsAllDigitsDouble(const std::string& str);
bool vectorOfNumberStringsInt(const VecStr& vec);
bool vectorOfNumberStringsDouble(const VecStr& vec);

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



template <typename T, typename TOT>
std::string getPercentageString(T partial, TOT total) {
  if (total == 0) {
    return "0";
  }
  std::stringstream outStream;
  outStream << std::setprecision(3);
  outStream << partial << "("
  		<< 100 * partial / static_cast<double>(total) << "%)";
  return outStream.str();
}

// combiners
std::string combineStrings(const VecStr& strings);

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

std::string intToHex(int i);
uint32_t hexToInt(const std::string& hString);

void rstrip(std::string & str, char c);
std::string rstripReturn(std::string str, char c);

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



}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "stringUtils.cpp"
#endif
