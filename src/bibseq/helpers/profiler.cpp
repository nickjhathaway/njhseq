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
