#pragma once
/*
 * clusterSet.hpp
 *
 *  Created on: Aug 7, 2016
 *      Author: nick
 */

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
			aligner& alignerObj, bool local, bool weighHomopolyers) {
		bool eventBased = true;
		for (auto& clus : clusters_) {
			std::string bestRef = profiler::getBestRef(refSeqs, clus, alignerObj,
					local, weighHomopolyers, eventBased, true, ",");
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

