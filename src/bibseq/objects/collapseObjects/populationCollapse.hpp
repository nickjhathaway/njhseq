#pragma once
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
//  populationCollapse.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/collapseObjects/sampleCollapse.hpp"

namespace bibseq {
namespace collapse {
class populationCollapse {

 public:
  populationCollapse(const std::vector<sampleCluster> &inputClusters,
                     const std::string &populationName);

  populationCollapse();
  // members
  // the initial clusters of the sample
  clusterSet input_;
  // collapsed clusters
  clusterSet collapsed_;

  std::string populationName_;

  uint32_t numOfSamps()const;
  // clustering
  void popCluster(collapser &collapserObj,
                  std::map<int, std::vector<double>> iteratorMap,
                  const std::string &sortBy, aligner &alignerObj);
  void popClusterOnId(collapser &collapserObj,
                  std::map<int, std::vector<double>> iteratorMap,
                  const std::string &sortBy, aligner &alignerObj);


  // update the initial infos
  void updateInitialInfos();
  // update the collapsed infos
  void updateCollapsedInfos();

  void renameToOtherPopNames(const std::vector<readObject> &previousPop,
                             aligner &alignerObj);
  void renameClusters();

  void updateInfoWithSampCollapses(const std::map<std::string, sampleCollapse> & sampCollapses);
  void updateInfoWithSampCollapse(const sampleCollapse & sampCollapses);

  //io
  void writeFinal(const std::string &outDirectory, const std::string &outFormat,
                  bool overWrite, bool exitOnFailureToWrite) const;
  void writeFinalInitial(const std::string &outDirectory, const readObjectIOOptions & ioOptions) const;




};
}  // namspace collpase
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "populationCollapse.cpp"
#endif
