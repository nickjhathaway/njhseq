#pragma once
//
//  sampleCollapse.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/collapseObjects/collapseCommon.h"
#include "bibseq/objects/collapseObjects/collapser.hpp"

namespace bibseq {
namespace collapse {

class sampleCollapse {

public:
	// constructor,
	// template<typename T>
	sampleCollapse(const std::vector<std::vector<cluster>> &inputClusters,
			const std::string &sampName, uint32_t sizeCutOff);

	sampleCollapse() {
	}
	sampleCollapse(const std::string & sampName) :
			sampName_(sampName) {
	}

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
	void cluster(const collapser &collapserObj, CollapseIterations iteratorMap,
			const std::string &sortBy, aligner &alignerObj);
	// excludes
	void excludeChimeras(bool update);
	void excludeChimeras(bool update, double fracCutOff);
	void excludeFraction(double fractionCutOff, bool update);
	void excludeBySampNum(uint32_t sampsRequired, bool update);
	//
	void renameClusters(const std::string &sortBy);
	// update the exclusioninfos
	void updateExclusionInfos();
	// update the initial infos
	void updateInitialInfos();
	// update the collapsed infos
	void updateCollapsedInfos();
	void updateAfterExclustion();
	// output reads
	std::vector<sampleCluster> createOutput(bool renameFirst,
			const std::string &sortBy);
	// write
	void writeExcluded(const std::string &outDirectory,
			const SeqIOOptions & ioOptions) const;
	void writeExcludedOriginalClusters(const std::string &outDirectory,
			const SeqIOOptions & ioOptions) const;
	void writeInitial(const std::string &outDirectory,
			const SeqIOOptions & ioOptions) const;
	void writeFinal(const std::string &outDirectory,
			const SeqIOOptions & ioOptions) const;
	void writeFinalOrignalClusters(const std::string &outDirectory,
			const SeqIOOptions & ioOptions) const;
	// output info
	std::string getSimpleSampInfo(const std::string &delim = "\t") const;
	static std::string getSimpleSampInfoHeader(const std::string & delim = "\t");

	VecStr getSimpleSampInfoVec() const;
	static VecStr getSimpleSampInfoHeaderVec();

	VecStr getAllInfoVec(bool checkingExpected,
			const std::string &delim = "\t") const;
	std::map<std::string, std::string, std::greater<std::string>> getAllInfoMap(
			bool checkingExpected, const std::string &delim = "\t",
			uint32_t maxRepCount = 0) const;

};
}  // namespace collpase
}  // namespace bibseq


