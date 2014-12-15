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

#include "alignmentInfo.hpp"

namespace bibseq {

VecStr readMatch::getInfoStrings() const {
  std::string firstString = matchedRead.seqBase_.name_ + "\t" +
                            std::to_string(fractionIdentity) + "\t" +
                            std::to_string(queryCoverage) + "\t" +
                            std::to_string(fractionGaps) + "\t" +
                            std::to_string(mismatches.size()) + "\t";

  if (mismatches.empty()) {
    firstString.append("\t\t\t\t\t\t\t\t\t\t\t\t");
  } else {
    firstString.append(mismatches.begin()->second.outputInfoString());
  }
  firstString.append("\t" + std::to_string(gaps.size()) + "\t");
  if (gaps.empty()) {
    firstString.append("\t\t\t\t\t\t");
  } else {
    firstString.append(gaps.begin()->second.outputGapInfoSingleLine());
  }
  VecStr ans(1, firstString);

  for (auto mIter = mismatches.begin(); mIter != mismatches.end(); ++mIter) {
    if (mIter == mismatches.begin()) {
      continue;
    }
    auto currentStr = "\t\t\t\t\t" + mIter->second.outputInfoString();
    ans.push_back(currentStr);
  }
  uint32_t counter = 2;
  for (auto gIter = gaps.begin(); gIter != gaps.end(); ++gIter) {
    if (gIter == gaps.begin()) {
      continue;
    }
    if (counter > ans.size()) {
      auto currentStr = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t" +
                        gIter->second.outputGapInfoSingleLine();

      ans.push_back(currentStr);
    } else {
      ans[counter - 1].append("\t\t" + gIter->second.outputGapInfoSingleLine());
    }
    ++counter;
  }

  return ans;
}

VecStr readMatch::getComparisonStrings(const readObject& readB,
                                       bool weighHomopolymers) const {
  double oneBaseIndels = 0;
  double twoBaseIndels = 0;
  double largeBaseIndels = 0;
  if (weighHomopolymers) {
    for (const auto& gIter : gaps) {
      if (gIter.second.size == 1) {
        oneBaseIndels += gIter.second.homoploymerScore;
      } else if (gIter.second.size == 2) {
        twoBaseIndels += gIter.second.homoploymerScore;
      } else {
        largeBaseIndels += gIter.second.homoploymerScore;
      }
    }
  } else {
    for (const auto& gIter : gaps) {
      if (gIter.second.size == 1) {
        ++oneBaseIndels;
      } else if (gIter.second.size == 2) {
        ++twoBaseIndels;
      } else {
        ++largeBaseIndels;
      }
    }
  }

  std::map<uint32_t, mismatchComparison> mComparisons;
  for (const auto& misIter : mismatches) {
    int lowestPrimaryQual = 99;
    if (misIter.second.refQual < misIter.second.seqQual) {
      lowestPrimaryQual = misIter.second.refQual;
    } else {
      lowestPrimaryQual = misIter.second.seqQual;
    }
    int lowestNeighborhoodQual = 99;
    int lowestReadANeighborhoodQual =
        matchedRead.seqBase_.findLowestNeighborhoodQual(
            misIter.second.refBasePos);
    int lowestReadBNeighborhoodQual =
        readB.seqBase_.findLowestNeighborhoodQual(misIter.second.seqBasePos);

    if (lowestReadANeighborhoodQual < lowestReadBNeighborhoodQual) {
      lowestNeighborhoodQual = lowestReadANeighborhoodQual;
    } else {
      lowestNeighborhoodQual = lowestReadBNeighborhoodQual;
    }
    mComparisons.insert(std::make_pair(
        misIter.first,
        mismatchComparison(lowestPrimaryQual, lowestNeighborhoodQual,
                           misIter.second.kMerFreqByPos,
                           misIter.second.kMerFreq)));
  }
  // readMatchComparison currentComparison=readMatchComparison(oneBaseIndels,
  // twoBaseIndels, largeBaseIndels, homopolymerScore, mComparisons);
  readMatchComparison currentComparison = readMatchComparison(
      oneBaseIndels, twoBaseIndels, largeBaseIndels, mComparisons);
  return currentComparison.getStringInfos();
}

void alignmentProfileInfo::outputReferenceComparisonInfo(std::ostream& out) {
  int count = 0;
  for (const auto& bIter : bestMatches) {
    if (count == 0) {
      out << readB.seqBase_.name_ << "\t" << readB.seqBase_.cnt_ << "\t"
          << readB.seqBase_.frac_ << "\t" << readB.seqBase_.seq_.length()
          << "\t";
    } else {
      out << readB.seqBase_.name_ << "\t\t\t\t";
    }
    ++count;
    VecStr currentInfos = bIter.getInfoStrings();
    for (const auto& sIter : currentInfos) {
      if (sIter == *currentInfos.begin()) {
        out << sIter << std::endl;
      } else {
        out << readB.seqBase_.name_ << "\t"
            << "\t\t\t" << bIter.matchedRead.seqBase_.name_ << sIter
            << std::endl;
      }
    }
  }
}

void alignmentProfileInfo::outputMatchComparisons(std::ostream& out,
                                                  bool weightHomopolyers) {
  int count = 0;
  for (const auto& bIter : bestMatches) {
    if (count == 0) {
      out << readB.seqBase_.name_ << "\t" << readB.seqBase_.cnt_ << "\t"
          << readB.seqBase_.frac_ << "\t" << readB.seqBase_.seq_.length()
          << "\t";
    } else {
      out << readB.seqBase_.name_ << "\t\t\t\t";
    }
    ++count;
    VecStr currentInfos = bIter.getComparisonStrings(readB, weightHomopolyers);
    for (const auto& sIter : currentInfos) {
      if (sIter == *currentInfos.begin()) {
        out << bIter.matchedRead.seqBase_.name_ << "\t" << sIter << std::endl;
      } else {
        out << readB.seqBase_.name_ << "\t"
            << "\t\t\t" << bIter.matchedRead.seqBase_.name_ << sIter
            << std::endl;
      }
    }
  }
}
}  // namespace bib
