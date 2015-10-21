//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//
#include "profiler.hpp"


namespace bibseq {

void profiler::outputHpRunsSimpleMap(
    const std::map<std::string, std::vector<hpRun>>& runs, std::ostream& out) {
  std::vector<hpRun>::const_iterator subIter;
  for (std::map<std::string, std::vector<hpRun>>::const_iterator runIter =
           runs.begin();
       runIter != runs.end(); ++runIter) {
    out << runIter->first << std::endl;
    std::map<int, int> sizeCount;
    std::map<int, int>::iterator sizeCountIter;
    std::map<int, std::vector<int>> sizePositions;
    for (subIter = runIter->second.begin(); subIter != runIter->second.end();
         ++subIter) {
      if (sizeCount.find(subIter->runSize) == sizeCount.end()) {
        sizeCount.insert(std::make_pair(subIter->runSize, subIter->count));
        std::vector<int> tempVec;
        tempVec.push_back(subIter->pos);
        sizePositions.insert(std::make_pair(subIter->runSize, tempVec));
      } else {
        sizeCount[subIter->runSize] += subIter->count;
        sizePositions[subIter->runSize].push_back(subIter->pos);
      }
    }
    for (sizeCountIter = sizeCount.begin(); sizeCountIter != sizeCount.end();
         ++sizeCountIter) {
      out << "\tsize: " << sizeCountIter->first
          << " count: " << sizeCountIter->second << std::endl;
      out << "\t\tpositions: ";
      std::sort(sizePositions[sizeCountIter->first].begin(),
                sizePositions[sizeCountIter->first].end());
      out << vectorToString(sizePositions[sizeCountIter->first], ",")
          << std::endl;
      // out<<"\n\t\t"<<vectorToString(sizePositions[sizeCountIter->first],"\n\t\t")<<std::endl;
    }
  }
}

}  // namespace bib
