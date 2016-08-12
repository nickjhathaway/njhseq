/*
 * clusterSet.cpp
 *
 *  Created on: Aug 7, 2016
 *      Author: nick
 */

#include "clusterSet.hpp"

namespace bibseq {
namespace collapse {

clusterSet::clusterSet() {

}

clusterSet::clusterSet(const std::vector<sampleCluster>& clusters) :
		clusters_(clusters), info_(clusters) {

}

uint32_t clusterSet::numOfReps() const{
	return info_.infos_.size();
}

void clusterSet::setSubClustersPositions() {
	uint32_t pos = 0;
	for (const auto& clus : clusters_) {
		for (const auto& read : clus.reads_) {
			subClustersPositions_[read->getStubName(true)] = pos;
		}
		++pos;
	}
}

void clusterSet::setSetInfo() {
	info_.clear();
	for (const auto& read : clusters_) {
		info_.updateInfo(read.reads_);
	}
}

// writing
void clusterSet::writeClusters(std::string filename,
		const SeqIOOptions & ioOptions) const {
	SeqOutput::write(clusters_, filename, ioOptions);
}

table clusterSet::getReplicateInfo() const {
	VecStr colNames { "clusterName", "repName", "fraction" };
	table ret(colNames);
	for (const auto & clus : clusters_) {
		for (const auto & s : info_.infos_) {
			auto search = clus.sampInfos().find(s.first);
			if (search == clus.sampInfos().end()) {
				ret.content_.emplace_back(toVecStr(clus.seqBase_.name_, s.first, 0));
			} else {
				ret.content_.emplace_back(
						toVecStr(clus.seqBase_.name_, s.first, search->second.fraction_));
			}
		}
	}
	return ret;
}

double clusterSet::getRMSE() const {
	double sumSquares = 0;

	for (const auto & clus : clusters_) {
		auto sampInfosNames = getVectorOfMapKeys(clus.sampInfos());
		std::vector<double> differences;
		for (const auto & pos : iter::range(sampInfosNames.size())) {
			for (const auto & subPos : iter::range(pos)) {
				if (pos == subPos) {
					continue;
				}
				differences.emplace_back(
						clus.sampInfos().at(sampInfosNames[pos]).fraction_
								- clus.sampInfos().at(sampInfosNames[subPos]).fraction_);
			}
		}
		sumSquares += std::pow(vectorMean(differences), 2);
	}
	return std::sqrt(sumSquares / clusters_.size());
}

} // namespace collapse
} // namespace bibseq
