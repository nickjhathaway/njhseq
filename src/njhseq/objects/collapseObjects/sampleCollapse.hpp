#pragma once
//
//  sampleCollapse.hpp
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
#include "njhseq/objects/collapseObjects/collapseCommon.h"
#include "njhseq/objects/collapseObjects/collapser.hpp"

namespace njhseq {
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

	void collapseLowFreqOneOffs(double lowFreqMultiplier, aligner &alignerObj, const collapser &collapserObj);
	// excludes
	void excludeChimeras(bool update);
	void excludeChimeras(bool update, double fracCutOff);
	void excludeChimerasNoReMark(bool update);
	void markChimeras(double fracCutOff);

	void excludeLowFreqOneOffs(bool update, double lowFreqMultiplier, aligner &alignerObj, bool skipChimeras = true);

	//void excludeFractionWithinRep(double fractionCutOff, bool update);
	void excludeFraction(double fractionCutOff, bool update);
	void excludeFractionAnyRep(double fractionCutOff, bool update);
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
}  // namespace njhseq


