#pragma once
//
//  populationCollapse.hpp
//
//  Created by Nicholas Hathaway on 9/13/13.
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


#include "njhseq/objects/collapseObjects/sampleCollapse.hpp"
#include "njhseq/objects/collapseObjects/collapser.hpp"


namespace njhseq {
namespace collapse {
class populationCollapse {

 public:
  populationCollapse(const std::vector<sampleCluster> &inputClusters,
                     const std::string &populationName);

  populationCollapse(const std::string &populationName);

  //for when populationCollapse is initizlized without inputClusters;
  void addInput(const std::vector<sampleCluster> &inputClusters);
  // members
  // the initial clusters of the sample
  clusterSet input_;
  // collapsed clusters
  clusterSet collapsed_;

  std::string populationName_;

  uint32_t numOfSamps()const;
  // clustering
  void popCluster(const collapser &collapserObj,
                  CollapseIterations iteratorMap,
                  const std::string &sortBy, aligner &alignerObj);

  // update the initial infos
  void updateInitialInfos();
  // update the collapsed infos
  void updateCollapsedInfos();

	void renameToOtherPopNames(const std::vector<readObject> &previousPop,
			comparison allowableErrors = comparison());
	void renameToOtherPopNames(const std::vector<readObject> &previousPop,
			aligner & alignerObj, comparison allowableErrors = comparison());
	void addRefMetaToName(const std::vector<readObject> &previousPop,
			comparison allowableErrors = comparison());
	void addRefMetaToName(const std::vector<readObject> &previousPop,
			aligner & alignerObj, comparison allowableErrors = comparison());
	void renameClusters();

  void updateInfoWithSampCollapses(const std::map<std::string, sampleCollapse> & sampCollapses);
  void updateInfoWithSampCollapse(const sampleCollapse & sampCollapses);

	//io
	void writeFinal(const std::string &outDirectory,
			const SeqIOOptions & ioOptions) const;
	void writeFinalInitial(const std::string &outDirectory,
			const SeqIOOptions & ioOptions) const;

	VecStr getPopInfoVec() const;
	static VecStr getPopInfoHeaderVec();

	double getExpectedHeterozygosity() const;


};
}  // namspace collpase
}  // namespace njhseq


