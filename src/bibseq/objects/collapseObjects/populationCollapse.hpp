#pragma once
//
//  populationCollapse.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/collapseObjects/sampleCollapse.hpp"
#include "bibseq/objects/collapseObjects/collapser.hpp"


namespace bibseq {
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
}  // namespace bibseq


