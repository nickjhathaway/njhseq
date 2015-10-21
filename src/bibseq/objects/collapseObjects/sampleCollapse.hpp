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
//
//  sampleCollapse.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "collapseCommon.hpp"

namespace bibseq {
namespace collapse {

class sampleCollapse {

 public:
  // constructor,
  // template<typename T>
	sampleCollapse(const std::vector<std::vector<cluster>> &inputClusters,
	               const std::string &sampName, uint32_t sizeCutOff);

  sampleCollapse() {}
  sampleCollapse(const std::string & sampName):sampName_(sampName) {}

  // members
  std::string sampName_;
  // the initial clusters of the sample
  clusterSet input_;
  // the excluded clusters
  clusterSet excluded_;
  // collapsed clusters
  clusterSet collapsed_;
  // functions
  // collapse the input clusters
  void cluster(collapser &collapserObj,
               std::map<int, std::vector<double>> iteratorMap,
               const std::string &sortBy, aligner &alignerObj);

  void clusterOnPerId(collapser &collapserObj,
               std::map<int, std::vector<double>> iteratorMap,
               const std::string &sortBy, aligner &alignerObj);
  // excludes
  void excludeChimeras(bool update);
  void excludeChimeras(bool update, double fracCutOff);
  void excludeFraction(double fractionCutOff, bool update);
  void excludeBySampNum(uint64_t sampsRequired, bool update);
  //
  void renameClusters(const std::string &sortBy);
  // update the exclusioninfos
  void updateExclusionInfos();
  // update the initial infos
  void updateInitialInfos() ;
  // update the collapsed infos
  void updateCollapsedInfos();
  void updateAfterExclustion() ;
  // output reads
  std::vector<sampleCluster> createOutput(bool renameFirst,
                                          const std::string &sortBy);
  // write
	void writeExcluded(const std::string &outDirectory,
			const readObjectIOOptions & ioOptions) const;
	void writeExcludedOriginalClusters(const std::string &outDirectory,
			const readObjectIOOptions & ioOptions) const;
	void writeInitial(const std::string &outDirectory,
			const readObjectIOOptions & ioOptions) const;
	void writeFinal(const std::string &outDirectory,
			const readObjectIOOptions & ioOptions) const;
	void writeFinalOrignalClusters(const std::string &outDirectory,
			const readObjectIOOptions & ioOptions) const;
  // output info
  std::string getSimpleSampInfo(uint32_t clusterId,const std::string &delim = "\t") const;
  static std::string getSimpleSampInfoHeader(const std::string & delim = "\t");
  VecStr getAllInfoVec(bool checkingExpected,
                       const std::string &delim = "\t") const;
  std::map<std::string, std::string, std::greater<std::string>> getAllInfoMap(
      bool checkingExpected, const std::string &delim = "\t", uint32_t maxRepCount = 0) const;





};
}  // namspace collpase
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "sampleCollapse.cpp"
#endif
