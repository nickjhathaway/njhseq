/*
 * HapsSampleClustering.cpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "HapsSampleClustering.hpp"
#include "njhseq/objects/dataContainers/graphs/UndirWeightedGraph.hpp"
#include "njhseq/concurrency/PairwisePairFactory.hpp"

namespace njhseq {

std::string HapsSampleClustering::measureName(Measure measure){
	switch (measure) {
		case Measure::ccc:
			return "ccc";
		case Measure::rmse:
			return "rmse";
		case Measure::jaccard:
			return "jaccard";
	}
	return "";
}

HapsSampleClustering::Measure HapsSampleClustering::measureFromName(const std::string & name){
	if("ccc" == name){
		return Measure::ccc;
	}
	if("rmse" == name){
		return Measure::rmse;
	}
	if("jaccard" == name){
		return Measure::jaccard;
	}
	std::stringstream ss;
	ss << __PRETTY_FUNCTION__ << ", error " << "unknown measure: " << name << ", options are ccc, rmse or jaccard" << "\n";
	throw std::runtime_error{ss.str()};
}

double HapsSampleClustering::Pars::eps() const{
	return Measure::rmse == measure_ ? cutOff_ : 1 - cutOff_;
}

HapsSampleClustering::Measures::Measures(HapsEncodedMatrix::CCCRMSEResults cccRmse, std::vector<std::vector<double>> jaccard):
		cccRmse_(std::move(cccRmse)), jaccard_(std::move(jaccard)){
}

HapsSampleClustering::Measures HapsSampleClustering::calcMeasures(HapsEncodedMatrix & haps, uint32_t binSize, bool verbose){
	haps.calcHapProbs();
	haps.add_relative_abundance();
	auto jaccard = std::move(haps.genIndexMeasures(binSize, verbose).byHapsTarShared);
	return Measures{haps.calc_ccc_rmse_measures(binSize, verbose), std::move(jaccard)};
}

std::vector<std::vector<double>> HapsSampleClustering::distances(const Measures & measures, Measure measure, uint32_t numThreads){
	if(Measure::rmse == measure){
		return measures.cccRmse_.rmse;
	}
	//CCC runs from -1 to 1 so 1 - CCC runs 0 to 2, 0 being CCC of 1
	auto ret = Measure::ccc == measure ? measures.cccRmse_.ccc : measures.jaccard_;
	PairwisePairFactory pairFactory(ret.size());
	uint32_t pairBatchCount = 100000;
	std::function<void()> transform = [&pairFactory, &pairBatchCount, &ret]() {
		PairwisePairFactory::PairwisePairVec pairs;
		while (pairFactory.setNextPairs(pairs, pairBatchCount)) {
			for (const auto & pair : pairs.pairs_) {
				ret[pair.row_][pair.col_] = -1 * (ret[pair.row_][pair.col_] - 1);
				ret[pair.col_][pair.row_] = ret[pair.row_][pair.col_];
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(transform, numThreads);
	for (uint32_t pos = 0; pos < ret.size(); ++pos) {
		ret[pos][pos] = 0;
	}
	return ret;
}

HapsSampleClustering::Results HapsSampleClustering::cluster(const HapsEncodedMatrix & haps, const Measures & measures, const Pars & pars){
	const auto distanceMatrix = distances(measures, pars.measure_, pars.numThreads_);
	const auto lociCoveragePerSample = haps.getTargetCoveragePerSample();
	const uint32_t sampleCount = distanceMatrix.size();
	Results ret;
	ret.lowCoverage_ = std::vector<bool>(sampleCount, false);
	for(uint32_t pos = 0; pos < sampleCount; ++pos){
		ret.lowCoverage_[pos] = lociCoveragePerSample.at(haps.sampNamesVec_[pos]) < pars.minLociCoverage_;
	}
	njhUndirWeightedGraph<double, uint32_t>::dbscanPars dbscanPars(pars.eps(), pars.minNeighbors_);
	njhUndirWeightedGraph<double, uint32_t> distGraph;
	for (uint32_t pos = 0; pos < sampleCount; ++pos) {
		distGraph.addNode(estd::to_string(pos), pos);
	}
	{
		struct PairDist {
			uint32_t row_;
			uint32_t col_;
			double dist_;
		};
		std::vector<PairDist> belowEps;
		std::mutex belowEpsMut;
		PairwisePairFactory pairFactory(sampleCount);
		uint32_t pairBatchCount = 100000;
		std::function<void()> gatherBelowEps = [&]() {
			PairwisePairFactory::PairwisePairVec pairs;
			std::vector<PairDist> currentBelowEps;
			while (pairFactory.setNextPairs(pairs, pairBatchCount)) {
				for (const auto & pair : pairs.pairs_) {
					if (ret.lowCoverage_[pair.row_] || ret.lowCoverage_[pair.col_]) {
						continue;
					}
					auto dist = distanceMatrix[pair.row_][pair.col_];
					if (dist <= dbscanPars.eps_) {
						currentBelowEps.emplace_back(PairDist{static_cast<uint32_t>(pair.row_), static_cast<uint32_t>(pair.col_), dist});
					}
				}
			}
			std::lock_guard<std::mutex> lock(belowEpsMut);
			addOtherVec(belowEps, currentBelowEps);
		};
		njh::concurrent::runVoidFunctionThreaded(gatherBelowEps, pars.numThreads_);
		//edges go in by pair so the clustering doesn't depend on thread timing
		njh::sort(belowEps, [](const PairDist & first, const PairDist & second){
			return first.row_ == second.row_ ? first.col_ < second.col_ : first.row_ < second.row_;
		});
		for (const auto & below : belowEps) {
			distGraph.addEdge(estd::to_string(below.row_), estd::to_string(below.col_), below.dist_);
		}
	}
	distGraph.dbscan(dbscanPars);
	if (dbscanPars.minEpNeighbors_ > 2) {
		distGraph.add_small_groups_in_off_nodes(dbscanPars);
	}
	ret.numberOfGroups_ = distGraph.numberOfGroups_;
	ret.groups_ = std::vector<uint32_t>(sampleCount, noGroup);
	for (const auto & node : distGraph.nodes_) {
		ret.groups_[node->value_] = node->group_;
		ret.groupIndexes_[node->group_].emplace_back(node->value_);
	}
	for (const auto & group : ret.groupIndexes_) {
		if (noGroup == group.first) {
			continue;
		}
		std::vector<double> cccsWithinGroup;
		std::vector<double> rmsesWithinGroup;
		std::vector<double> jaccardsWithinGroup;
		PairwisePairFactory pfac(group.second.size());
		PairwisePairFactory::PairwisePair pair;
		while (pfac.setNextPair(pair)) {
			cccsWithinGroup.emplace_back(measures.cccRmse_.ccc[group.second[pair.col_]][group.second[pair.row_]]);
			rmsesWithinGroup.emplace_back(measures.cccRmse_.rmse[group.second[pair.col_]][group.second[pair.row_]]);
			jaccardsWithinGroup.emplace_back(measures.jaccard_[group.second[pair.col_]][group.second[pair.row_]]);
		}
		ret.cccStats_[group.first] = getStatsOnVec(cccsWithinGroup);
		ret.rmseStats_[group.first] = getStatsOnVec(rmsesWithinGroup);
		ret.jaccardStats_[group.first] = getStatsOnVec(jaccardsWithinGroup);
	}
	for (const auto & e : distGraph.edges_) {
		if (e->on_) {
			auto node1 = e->nodeToNode_.begin()->second.lock();
			auto node2 = e->nodeToNode_.rbegin()->second.lock();
			ret.edges_.emplace_back(Edge{node1->value_, node2->value_, node1->group_});
		}
	}
	return ret;
}

std::vector<uint32_t> HapsSampleClustering::averageLinkageOrder(const std::vector<uint32_t> & items,
		const std::vector<std::vector<double>> & distances, double missingDistance){
	const uint32_t n = items.size();
	if (n < 3) {
		return items;
	}
	std::vector<double> dist(static_cast<uint64_t>(n) * n, 0);
	for (uint32_t i = 0; i < n; ++i) {
		for (uint32_t j = i + 1; j < n; ++j) {
			auto d = distances[items[i]][items[j]];
			if (std::isnan(d)) {
				d = missingDistance;
			}
			dist[static_cast<uint64_t>(i) * n + j] = d;
			dist[static_cast<uint64_t>(j) * n + i] = d;
		}
	}
	std::vector<uint32_t> size(n, 1);
	std::vector<bool> active(n, true);
	std::vector<uint32_t> nodeOf(n);
	for (uint32_t i = 0; i < n; ++i) {
		nodeOf[i] = i;
	}
	std::vector<uint32_t> lefts;
	std::vector<uint32_t> rights;
	std::vector<uint32_t> chain;
	uint32_t remaining = n;
	uint32_t firstActive = 0;
	while (remaining > 1) {
		if (chain.empty()) {
			while (!active[firstActive]) {
				++firstActive;
			}
			chain.emplace_back(firstActive);
		}
		const uint32_t a = chain.back();
		const int64_t prev = chain.size() > 1 ? chain[chain.size() - 2] : -1;
		int64_t best = -1;
		double bestDist = std::numeric_limits<double>::infinity();
		for (uint32_t k = 0; k < n; ++k) {
			if (!active[k] || k == a) {
				continue;
			}
			const auto d = dist[static_cast<uint64_t>(a) * n + k];
			//ties go to the previous link so the chain can't cycle
			if (d < bestDist || (d == bestDist && k == prev)) {
				bestDist = d;
				best = k;
			}
		}
		if (best == prev) {
			chain.pop_back();
			chain.pop_back();
			const auto prevPos = static_cast<uint32_t>(prev);
			for (uint32_t k = 0; k < n; ++k) {
				if (!active[k] || k == a || k == prevPos) {
					continue;
				}
				const auto merged = (size[a] * dist[static_cast<uint64_t>(a) * n + k] + size[prevPos] * dist[static_cast<uint64_t>(prevPos) * n + k]) / (size[a] + size[prevPos]);
				dist[static_cast<uint64_t>(a) * n + k] = merged;
				dist[static_cast<uint64_t>(k) * n + a] = merged;
			}
			lefts.emplace_back(nodeOf[a]);
			rights.emplace_back(nodeOf[prevPos]);
			nodeOf[a] = n + lefts.size() - 1;
			size[a] += size[prevPos];
			active[prevPos] = false;
			--remaining;
		} else {
			chain.emplace_back(static_cast<uint32_t>(best));
		}
	}
	std::vector<uint32_t> ret;
	ret.reserve(n);
	std::vector<uint32_t> stack{static_cast<uint32_t>(n + lefts.size() - 1)};
	while (!stack.empty()) {
		const auto node = stack.back();
		stack.pop_back();
		if (node < n) {
			ret.emplace_back(items[node]);
		} else {
			stack.emplace_back(rights[node - n]);
			stack.emplace_back(lefts[node - n]);
		}
	}
	return ret;
}

}  // namespace njhseq
