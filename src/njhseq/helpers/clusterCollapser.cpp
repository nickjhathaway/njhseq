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
#include "clusterCollapser.hpp"

namespace njhseq {

void clusterCollapser::collapseTandems(std::vector<cluster> &processedReads,
                                       aligner &alignerObj,
                                       const collapseTandemsPars & pars) {

  int counter = 0;
  //auto kMaps = kmerCalculator::indexKmerMpas(processedReads, kLength, runCutOff,
  //                                           runCutOff);
  //auto currentKMaps = alignerObj.getKmerMaps();
  //alignerObj.setKmerMpas(kMaps);
  for (auto &clusIter : iter::reversed(processedReads)) {
    ++counter;
    if (pars.verbose && counter % 20 == 0) {
      std::cout << counter << " out of " << processedReads.size() << std::endl;
    }
    //assumes this is sorted by counts, so will skip if the second cluster doesn't even meet cluster freq difference multiplier
    for (auto &clusIterSecond : processedReads) {
      if (clusIter.seqBase_.name_ == clusIterSecond.seqBase_.name_ ||
        clusIterSecond.seqBase_.cnt_ < (pars.freqCutoff * clusIter.seqBase_.cnt_) ||
              clusIterSecond.remove) {
        continue;
      }
      alignerObj.alignCacheGlobal(clusIterSecond, clusIter);
      alignerObj.profileAlignment(clusIterSecond, clusIter, true, true, false);
      // alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
      // alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
      // std::cout << "alignerObj.checkForTandemRepeatGap(): " << njh::colorBool(alignerObj.checkForTandemRepeatGap()) << std::endl;

      auto tandems = alignerObj.getTandemRepeatGapsForCurrentAlignment();
      // std::cout << "tandems.size() : " << tandems.size()<< std::endl;
      // std::cout << "alignerObj.comp_.distances_.alignmentGaps_.size(): " << alignerObj.comp_.distances_.alignmentGaps_.size() << std::endl;
      // //alignerObj.checkForTandemRepeatGap()
      // std::cout << "tandems.size() <= pars.allowableTandems : " << njh::colorBool(tandems.size() <= pars.allowableTandems)<< std::endl;
      // std::cout << "alignerObj.comp_.distances_.alignmentGaps_.size() == tandems.size() : " << njh::colorBool(alignerObj.comp_.distances_.alignmentGaps_.size() == tandems.size())<< std::endl;
      // std::cout << "alignerObj.comp_.hqMismatches_ <= pars.allowableMismatches.hqMismatches_ : " << njh::colorBool(alignerObj.comp_.hqMismatches_ <= pars.allowableMismatches.hqMismatches_)<< std::endl;
      // std::cout << "alignerObj.comp_.lqMismatches_ <= pars.allowableMismatches.lqMismatches_ : " << njh::colorBool(alignerObj.comp_.lqMismatches_ <= pars.allowableMismatches.lqMismatches_)<< std::endl;

      if (tandems.size() <= pars.allowableTandems &&
          alignerObj.comp_.distances_.alignmentGaps_.size() == tandems.size() &&
          alignerObj.comp_.hqMismatches_ <= pars.allowableMismatches.hqMismatches_ &&
          alignerObj.comp_.lqMismatches_ <= pars.allowableMismatches.lqMismatches_ &&
          clusIterSecond.seqBase_.cnt_ >= (pars.freqCutoff * clusIter.seqBase_.cnt_)) {
        clusIter.remove = true;
        clusIterSecond.addRead(clusIter);
        break;
      }
    }
  }
}
}  // namespace njh
