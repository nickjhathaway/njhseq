//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "clusterCollapser.hpp"

namespace bibseq {

void clusterCollapser::collapseTandems(std::vector<cluster> &processedReads,
                                       aligner &alignerObj, int runCutOff,
                                       int kLength, bool kMersByPosition,
                                       double freqCutoff, bool local,
                                       bool weighHomopolyer) {

  int counter = 0;
  //auto kMaps = kmerCalculator::indexKmerMpas(processedReads, kLength, runCutOff,
  //                                           runCutOff);
  //auto currentKMaps = alignerObj.getKmerMaps();
  //alignerObj.setKmerMpas(kMaps);
  for (auto &clusIter : iter::reverse(processedReads)) {
    ++counter;
    if (counter % 20 == 0) {
      std::cout << counter << " out of " << processedReads.size() << std::endl;
    }
    for (auto &clusIterSecond : processedReads) {
      if (clusIter.seqBase_.name_ == clusIterSecond.seqBase_.name_) {
        continue;
      }
      alignerObj.alignCache(clusIterSecond, clusIter, local);
      alignerObj.profileAlignment(clusIterSecond, clusIter, kLength,
                                  kMersByPosition, true, true, false,
                                  weighHomopolyer);
      // alignerObj.outPutParameterInfo(std::cout);
      if (alignerObj.checkForTandemRepeatGap() &&
          alignerObj.comp_.hqMismatches_ < 1 &&
          alignerObj.comp_.lqMismatches_ < 1 &&
          clusIterSecond.seqBase_.cnt_ >
              (freqCutoff * clusIter.seqBase_.cnt_)) {
        clusIter.remove = true;
        clusIterSecond.addRead(clusIter);
        break;
      }
    }
  }
  //alignerObj.setKmerMpas(currentKMaps);
  return;
}
}  // namespace bib
