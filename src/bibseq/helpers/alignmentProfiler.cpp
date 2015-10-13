#include "alignmentProfiler.hpp"
#include "bibseq/readVectorManipulation/readVectorOperations/massGetters.hpp"

namespace bibseq {

void alignmentProfiler::getAlignmentInformationForReferenceRawReads(
    const std::vector<readObject>& reads,
    const std::vector<readObject>& refSeqs, aligner& alignerObj, bool local,
    const std::string& fileName, int kLength, bool weighHomopolyers) {
  std::ofstream profile;
  openTextFile(profile, fileName, ".tab.txt", false, false);
  profile << "Read\tclusterSize\tclusterFraction\treadLength\tbestRef\tid"
             "entity\tcoverage\tgaps\tmismatches\ttransiti"
             "onOrTransversion\trefBasePos\trefBase\trefQual\trefLeadQual\tre"
             "fTrailQual\tseqBasePos\tseqBa"
             "se\tseqQual\tseqLeadQual\tseqTrailQual\tkmerFreqByPos\tkmerFreq"
             "\tgaps\trefOrRead\tpos\tgaped"
             "Seq\tsumQual\tsize\thpScore" << std::endl;
  std::ofstream comparison;

  openTextFile(comparison, fileName + "_comparison", ".tab.txt", false, false);
  comparison << "Read\tclusterSize\tclusterFraction\tbestMatch\toneBaseIndels\t"
                "twoBaseIndels\tlargeIndels\tmismatchePQual\tmismatchNQual\tkme"
                "rPosFreq\tkmerAnywhereFreq" << std::endl;

  int counter = 1;
  MapStrStr infos;
  MapStrStr comparisonInfos;
  // std::map<std::string, std::stringstream> exactMatcheInfos;
  std::map<std::string, alignmentProfileInfo> exactMatches;
  int readAmount = readVec::getTotalReadCount(reads);
  for (auto rIter = refSeqs.begin(); rIter != refSeqs.end(); ++rIter) {
    readMatch currentMatch =
        readMatch(rIter->seqBase_.name_, *rIter, alignerObj.alignmentGaps_,
                  alignerObj.mismatches_, alignerObj.lowKmerMismatches_,
                  alignerObj.comp_.distances_.percentIdentity_, alignerObj.comp_.distances_.queryCoverage_,
                  alignerObj.comp_.distances_.percentageGaps_);
    exactMatches.insert(std::make_pair(
        rIter->seqBase_.name_, alignmentProfileInfo(*rIter, currentMatch)));
    exactMatches[rIter->seqBase_.name_].readB.seqBase_.name_ = "exactMatch";
    exactMatches[rIter->seqBase_.name_].readB.seqBase_.cnt_--;
  }
  for (std::vector<readObject>::const_iterator rIter = reads.begin();
       rIter != reads.end(); ++rIter) {
    if (counter % 50 == 0) {
      std::cout << "Currently on " << counter << " of " << reads.size()
                << std::endl;
    }
    ++counter;
    double bestScore = 0.00;
    std::vector<readObject> bestRefs;
    std::stringstream currentOut;
    std::stringstream currentComparisonOut;
    for (std::vector<readObject>::const_iterator refIter = refSeqs.begin();
         refIter != refSeqs.end(); ++refIter) {
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

    alignmentProfileInfo currentInfo = alignmentProfileInfo(*rIter);

    for (const auto& best : bestRefs) {
      alignerObj.alignVec(best, *rIter, local);
      alignerObj.profileAlignment(best, *rIter, kLength, true, false, true,
                                  false, weighHomopolyers);
      if (alignerObj.mismatches_.size() == 0 &&
          alignerObj.alignmentGaps_.size() == 0) {
        exactMatches[best.seqBase_.name_].readB.seqBase_.cnt_ +=
            rIter->seqBase_.cnt_ / bestRefs.size();
      } else {
        readMatch currentMatch =
            readMatch(rIter->seqBase_.name_, best, alignerObj.alignmentGaps_,
                      alignerObj.mismatches_, alignerObj.lowKmerMismatches_,
                      alignerObj.comp_.distances_.percentIdentity_, alignerObj.comp_.distances_.queryCoverage_,
                      alignerObj.comp_.distances_.percentageGaps_);
        currentInfo.bestMatches.push_back(currentMatch);
      }
    }
    currentInfo.outputReferenceComparisonInfo(currentOut);
    currentInfo.outputMatchComparisons(currentComparisonOut, weighHomopolyers);
    infos.insert(std::make_pair(rIter->seqBase_.name_, currentOut.str()));
    comparisonInfos.insert(
        std::make_pair(rIter->seqBase_.name_, currentComparisonOut.str()));
  }
  std::map<std::string, alignmentProfileInfo>::iterator exactIter;
  std::cout << exactMatches.size() << std::endl;
  for (exactIter = exactMatches.begin(); exactIter != exactMatches.end();
       ++exactIter) {
    exactIter->second.readB.setFractionByCount(readAmount);
    exactIter->second.outputReferenceComparisonInfo(profile);
    exactIter->second.outputMatchComparisons(comparison, weighHomopolyers);
  }
  for (const auto& cIter : comparisonInfos) {
    comparison << cIter.second;
  }
  MapStrStr::iterator infoIter;
  for (infoIter = infos.begin(); infoIter != infos.end(); ++infoIter) {
    profile << infoIter->second;
  }
}
}  // namespace bib
