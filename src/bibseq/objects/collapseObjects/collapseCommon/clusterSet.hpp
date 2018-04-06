#pragma once
/*
 * clusterSet.hpp
 *
 *  Created on: Aug 7, 2016
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bibseq/objects/collapseObjects/collapseCommon/clusterSetInfo.hpp"
#include "bibseq/IO/SeqIO/SeqIOOptions.hpp"
#include "bibseq/objects/seqObjects/Clusters/sampleCluster.hpp"

#include "bibseq/helpers/profiler.hpp"

namespace bibseq {
namespace collapse {

class clusterSet {
public:
	// constructors
	clusterSet();
	clusterSet(const std::vector<sampleCluster>& clusters);
	// members
	std::vector<sampleCluster> clusters_;
	std::unordered_map<std::string, uint32_t> subClustersPositions_;

	clusterSetInfo info_;

	uint32_t numOfReps() const;

	void setSubClustersPositions();
	void setSetInfo();

	template<typename REF>
	void checkAgainstExpected(const std::vector<REF>& refSeqs,
			aligner& alignerObj, bool local) {
		bool eventBased = true;
		for (auto& clus : clusters_) {

			std::string bestRef = bib::conToStr(
					profiler::compareToRefSingle(refSeqs, clus, alignerObj, local,
							eventBased), ";");

			clus.expectsString = bestRef;

		}
	}


	void writeClusters(std::string filename,
			const SeqIOOptions & ioOptions) const;

	table getReplicateInfo() const;

	double getRMSE() const;
};

} // namespace collapse
} // namespace bibseq

