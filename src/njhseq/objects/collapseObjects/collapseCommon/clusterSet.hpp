#pragma once
/*
 * clusterSet.hpp
 *
 *  Created on: Aug 7, 2016
 *      Author: nick
 */
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
#include "njhseq/objects/collapseObjects/collapseCommon/clusterSetInfo.hpp"
#include "njhseq/IO/SeqIO/SeqIOOptions.hpp"
#include "njhseq/objects/seqObjects/Clusters/sampleCluster.hpp"

#include "njhseq/helpers/profiler.hpp"

namespace njhseq {
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

			std::string bestRef = njh::conToStr(
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
} // namespace njhseq

