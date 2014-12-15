#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bibseq/objects/helperObjects/gaps.hpp"
#include "bibseq/objects/helperObjects/mismatch.hpp"
#include "bibseq/objects/seqObjects/readObject.hpp"

namespace bibseq {

class mismatchComparison {

 public:
  mismatchComparison()
      : lowestPrimaryQual(0),
        lowestNeighborhoodQual(0),
        kmerPosFreq(0),
        kmerAnywhereFreq(0) {}
  mismatchComparison(int pQual, int lowestNQ, int kPosFreq, int kAnywhereFreq)
      : lowestPrimaryQual(pQual),
        lowestNeighborhoodQual(lowestNQ),
        kmerPosFreq(kPosFreq),
        kmerAnywhereFreq(kAnywhereFreq) {}
  int lowestPrimaryQual;
  int lowestNeighborhoodQual;
  int kmerPosFreq;
  int kmerAnywhereFreq;
  std::string getStringInfo(const std::string& delim = "\t") const {
    std::stringstream out;
    out << lowestPrimaryQual << delim << lowestNeighborhoodQual << delim
        << kmerPosFreq << delim << kmerAnywhereFreq;
    return out.str();
  }
};

class readMatchComparison {

 public:
  readMatchComparison() {}
  // readMatchComparison(int oneBIndels, int twoBIndels, int lIndels, double
  // hScore, const std::map<uint32_t, mismatchComparison> & mInfos ):
  // oneBaseIndels(oneBIndels), twoBaseIndels(twoBIndels), largeIndels(lIndels),
  // homopolymerScore(hScore), mismatchComparisons(mInfos){}
  readMatchComparison(double oneBIndels, double twoBIndels, double lIndels,
                      const std::map<uint32_t, mismatchComparison>& mInfos)
      : oneBaseIndels(oneBIndels),
        twoBaseIndels(twoBIndels),
        largeIndels(lIndels),
        mismatchComparisons(mInfos) {}
  double oneBaseIndels;
  double twoBaseIndels;
  double largeIndels;
  // double homopolymerScore;
  std::map<uint32_t, mismatchComparison> mismatchComparisons;
  VecStr getStringInfos(const std::string& delim = "\t") const {
    std::stringstream firstLine;
    // firstLine<<oneBaseIndels<<delim<<twoBaseIndels<<delim<<largeIndels<<delim<<homopolymerScore;
    firstLine << oneBaseIndels << delim << twoBaseIndels << delim
              << largeIndels;
    VecStr ans;

    if (mismatchComparisons.size() == 0) {
      firstLine << delim << delim << delim << delim;
    } else {
      firstLine << delim
                << (*mismatchComparisons.begin()).second.getStringInfo(delim);
    }
    ans.push_back(firstLine.str());
    int count = 0;
    for (auto& mComIter : mismatchComparisons) {
      if (count != 0) {
        std::stringstream currentLine;
        currentLine << delim << delim << delim << delim;
        currentLine << mComIter.second.getStringInfo(delim);
        ans.push_back(currentLine.str());
      }
      ++count;
    }
    return ans;
  }
};

class readMatch {

 public:
  readMatch(const std::string& otherName, const readObject& read,
            const std::map<uint32_t, gap>& inGaps,
            const std::map<uint32_t, mismatch>& inMismatches,
            const std::map<uint32_t, mismatch>& inLowKmerMismatches,
            double identity, double coverage, double gaps)
      : otherReadName(otherName),
        matchedRead(read),
        gaps(inGaps),
        mismatches(inMismatches),
        lowKmerMismatches(inLowKmerMismatches),
        fractionIdentity(identity),
        queryCoverage(coverage),
        fractionGaps(gaps) {}
  std::string otherReadName;
  readObject matchedRead;
  std::map<uint32_t, gap> gaps;
  std::map<uint32_t, mismatch> mismatches;
  std::map<uint32_t, mismatch> lowKmerMismatches;
  double fractionIdentity;
  double queryCoverage;
  double fractionGaps;
  VecStr getInfoStrings() const;
  VecStr getComparisonStrings(const readObject& readB,
                              bool weighHomopolymers) const;
  size_t getMaxLines() {
    if (mismatches.size() > gaps.size()) {
      return mismatches.size();
    } else {
      return gaps.size();
    }
  }
};

class alignmentProfileInfo {
 public:
  alignmentProfileInfo() {}
  alignmentProfileInfo(const readObject& read) : readB(read) {}
  alignmentProfileInfo(const readObject& read, const readMatch& firstReadMatch)
      : readB(read), bestMatches(std::vector<readMatch>(1, firstReadMatch)) {}
  // convention is readA will be a reference sequence and readB is a real
  // sequence
  readObject readB;
  std::vector<readMatch> bestMatches;
  void addReadMatch(const readMatch& otherBestMatch) {
    bestMatches.push_back(otherBestMatch);
  }

  void outputReferenceComparisonInfo(std::ostream& out);
  void outputMatchComparisons(std::ostream& out, bool weighHomopolymers);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "alignmentInfo.cpp"
#endif
