#pragma once
//
//  alignmentProfiler.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/3/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/alignment/aligner.hpp"
#include "bibseq/alignment/alignmentInfo.hpp"
#include "bibseq/readVectorManipulation/readVectorOperations.h"

namespace bibseq {

class alignmentProfiler {
 public:
  template <class READS, class REFS>
  static void getAlignmentInformationForReference(
      const std::vector<READS>& reads, const std::vector<REFS>& refSeqs,
      aligner& alignerObj, bool local, const std::string& fileName,
      bool doingKmerChecking, int kLength, bool kmersByPosition,
      bool weighHomopolyers) {
    std::ofstream profile;
    openTextFile(profile, fileName, ".tab.txt", false, false);
    profile << "Read\tclusterSize\tclusterFraction\treadLength\tbestRef\tid"
               "entity\tcoverage\tgaps\tmismatches\ttransiti"
               "onOrTransversion\trefBasePos\trefBase\trefQual\trefLeadQual\tre"
               "fTrailQual\tseqBasePos\tseqBa"
               "se\tseqQual\tseqLeadQual\tseqTrailQual\tkmerFreqByPos\tkmerFreq"
               "\tgaps\trefOrRead\tpos\tgaped"
               "Seq\tsumQual\tsize\thpScore" << std::endl;
    int counter = 1;

    for (auto rIter = reads.begin(); rIter != reads.end(); ++rIter) {
      if (counter % 50 == 0) {
        std::cout << "Currently on " << counter << " of " << reads.size()
                  << std::endl;
      }
      ++counter;
      double bestScore = 0.00;
      std::vector<REFS> bestRefs;

      for (auto refIter = refSeqs.begin(); refIter != refSeqs.end();
           ++refIter) {
        if (refIter->seqBase_.name_ == rIter->seqBase_.name_) {
          continue;
        }
        alignerObj.alignVec(*refIter, *rIter, local);
        if (alignerObj.parts_.score_ == bestScore && alignerObj.parts_.score_ != 0) {
          bestRefs.push_back(*refIter);
        } else if (alignerObj.parts_.score_ > bestScore) {
          bestScore = alignerObj.parts_.score_;
          bestRefs.clear();
          bestRefs.push_back(*refIter);
        }
      }

      auto currentInfo = alignmentProfileInfo(*rIter);

      for (auto bestIter = bestRefs.begin(); bestIter != bestRefs.end();
           ++bestIter) {

        alignerObj.alignVec(*bestIter, *rIter, local);
        alignerObj.profileAlignment(*bestIter, *rIter, kLength, kmersByPosition,
                                    doingKmerChecking, true, true,
                                    weighHomopolyers);
        readMatch currentMatch = readMatch(
            rIter->seqBase_.name_, *bestIter, alignerObj.alignmentGaps_,
            alignerObj.mismatches_, alignerObj.lowKmerMismatches_,
            alignerObj.distances_.percentIdentity_, alignerObj.distances_.queryCoverage_,
            alignerObj.distances_.percentageGaps_);
        currentInfo.bestMatches.push_back(currentMatch);
      }
      currentInfo.outputReferenceComparisonInfo(profile);
    }
  }

  template <class CLUS>
  static void getAlignmentInformationForEveryComparison(
      const std::vector<CLUS>& reads, aligner& alignerObj, bool local,
      const std::string& fileName, bool doingKmerChecking, int kLength,
      bool kmersByPosition, bool weighHomopolyers) {

    std::ofstream profile;
    openTextFile(profile, fileName, ".tab.txt", false, false);
    profile << "Read\tclusterSize\tclusterFraction\treadLength\tbestRef\tid"
               "entity\tcoverage\tgaps\tmismatches\ttransiti"
               "onOrTransversion\trefBasePos\trefBase\trefQual\trefLeadQual\tre"
               "fTrailQual\tseqBasePos\tseqBa"
               "se\tseqQual\tseqLeadQual\tseqTrailQual\tkmerFreqByPos\tkmerFreq"
               "\tgaps\trefOrRead\tpos\tgaped"
               "Seq\tsumQual\tsize\thpScore" << std::endl;
    int counter = 1;

    for (auto rIter = reads.begin(); rIter != reads.end(); ++rIter) {

      if (counter % 50 == 0) {
        std::cout << "Currently on " << counter << " of " << reads.size()
                  << std::endl;
      }
      ++counter;
      alignmentProfileInfo currentInfo = alignmentProfileInfo(*rIter);

      for (auto rIterSecond = rIter; rIterSecond != reads.end();
           ++rIterSecond) {
        if (rIter == rIterSecond) {
          continue;
        }
        alignerObj.alignVec(*rIterSecond, *rIter, local);
        alignerObj.profileAlignment(*rIterSecond, *rIter, kLength,
                                    kmersByPosition, doingKmerChecking, true,
                                    false, weighHomopolyers);
        readMatch currentMatch = readMatch(
            rIter->seqBase_.name_, *rIterSecond, alignerObj.alignmentGaps_,
            alignerObj.mismatches_, alignerObj.lowKmerMismatches_,
            alignerObj.distances_.percentIdentity_, alignerObj.distances_.queryCoverage_,
            alignerObj.distances_.percentageGaps_);
        currentInfo.bestMatches.push_back(currentMatch);
      }
      currentInfo.outputReferenceComparisonInfo(profile);
    }
  }

  static void getAlignmentInformationForReferenceRawReads(
      const std::vector<readObject>& reads,
      const std::vector<readObject>& refSeqs, aligner& alignerObj, bool local,
      const std::string& fileName, int kLength, bool weighHomopolyers);

  template <typename READS, typename REF>
  static void getInfoSingleComparison(const std::vector<READS>& reads,
                                      const REF& refseq, aligner& alignerObj,
                                      bool local, const std::string& fileName,
                                      int kLength, bool weighHomopolyers) {
    std::ofstream profile;
    openTextFile(profile, fileName, ".tab.txt", false, false);
    profile << "Read\tclusterSize\tclusterFraction\treadLength\tclusterName\tid"
               "entity\tcoverage\tgaps\tmismatches\ttransiti"
               "onOrTransversion\trefBasePos\trefBase\trefQual\trefLeadQual\tre"
               "fTrailQual\tseqBasePos\tseqBa"
               "se\tseqQual\tseqLeadQual\tseqTrailQual\tkmerFreqByPos\tkmerFreq"
               "\tgaps\trefOrRead\tpos\tgaped"
               "Seq\tsumQual\tsize\thpScore" << std::endl;
    alignerObj.mismatches_.clear();
    alignerObj.alignmentGaps_.clear();
    alignerObj.lowKmerMismatches_.clear();
    alignerObj.distances_.percentIdentity_ = 1;
    alignerObj.distances_.queryCoverage_ = 1;
    alignerObj.distances_.percentageGaps_ = 0;
    int counter = 1;
    std::multimap<double, std::string, std::greater<double>> infos;
    int readAmount = readVec::getTotalReadCount(reads);
    readMatch currentMatch =
        readMatch(refseq.seqBase_.name_, refseq, alignerObj.alignmentGaps_,
                  alignerObj.mismatches_, alignerObj.lowKmerMismatches_,
                  alignerObj.distances_.percentIdentity_, alignerObj.distances_.queryCoverage_,
                  alignerObj.distances_.percentageGaps_);
    auto exactMatch = alignmentProfileInfo(refseq, currentMatch);
    exactMatch.readB.seqBase_.name_ =
        "exactMatch_" + exactMatch.readB.seqBase_.name_;
    exactMatch.readB.seqBase_.cnt_ = 0;
    for (const auto& read : reads) {
      /*if (read.seqBase_.name_ == refseq.seqBase_.name_) {
        continue;
      }*/
      if (counter % 50 == 0) {
        std::cout << "Currently on " << counter << " of " << reads.size()
                  << std::endl;
      }
      ++counter;
      std::stringstream currentOut;
      alignmentProfileInfo currentInfo = alignmentProfileInfo(read);
      alignerObj.alignVec(refseq, read, local);
      alignerObj.profileAlignment(refseq, read, kLength, true, false, true,
                                  false, weighHomopolyers);
      if (alignerObj.mismatches_.size() == 0 &&
          alignerObj.alignmentGaps_.size() == 0) {
        exactMatch.readB.seqBase_.cnt_ += read.seqBase_.cnt_;
      } else {
        readMatch currentMatch =
            readMatch(read.seqBase_.name_, refseq, alignerObj.alignmentGaps_,
                      alignerObj.mismatches_, alignerObj.lowKmerMismatches_,
                      alignerObj.distances_.percentIdentity_, alignerObj.distances_.queryCoverage_,
                      alignerObj.distances_.percentageGaps_);
        currentInfo.bestMatches.push_back(currentMatch);
      }
      currentInfo.outputReferenceComparisonInfo(currentOut);
      infos.insert({read.seqBase_.cnt_, currentOut.str()});
    }
    exactMatch.readB.setFractionByCount(readAmount);
    exactMatch.outputReferenceComparisonInfo(profile);
    for (const auto& info : infos) {
      profile << info.second;
    }
  }
};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "alignmentProfiler.cpp"
#endif
