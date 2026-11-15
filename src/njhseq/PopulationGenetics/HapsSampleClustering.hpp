#pragma once

/*
 * HapsSampleClustering.hpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "njhseq/PopulationGenetics/HapsEncodedMatrix.hpp"

namespace njhseq {

/**@brief Between sample CCC, RMSE and Jaccard on the targets two samples share, and DBSCAN clustering of samples on one of them
 *
 */
class HapsSampleClustering {
public:
	enum class Measure {
		ccc,
		rmse,
		jaccard
	};

	/**@brief ccc, rmse or jaccard
	 *
	 */
	[[nodiscard]] static std::string measureName(Measure measure);
	/**@brief from ccc, rmse or jaccard, throws on anything else
	 *
	 */
	[[nodiscard]] static Measure measureFromName(const std::string & name);

	struct Pars {
		Measure measure_ = Measure::ccc;
		double cutOff_ = 0.95; /**< samples are neighbors with CCC or Jaccard at or above this, or RMSE at or below it*/
		uint32_t minNeighbors_ = 5; /**< the DBSCAN minimum neighbors*/
		double minLociCoverage_ = 0.90; /**< samples with data for less than this fraction of the targets aren't clustered*/
		uint32_t numThreads_ = 1;

		/**@brief the DBSCAN distance, 1 - cutOff_ for CCC and Jaccard, cutOff_ for RMSE
		 *
		 */
		[[nodiscard]] double eps() const;
	};

	struct Measures {
		Measures(HapsEncodedMatrix::CCCRMSEResults cccRmse, std::vector<std::vector<double>> jaccard);
		HapsEncodedMatrix::CCCRMSEResults cccRmse_;
		std::vector<std::vector<double>> jaccard_; /**< HapsEncodedMatrix::IndexResults::byHapsTarShared, NaN for samples sharing no targets*/
	};

	/**@brief calc relative abundances on haps, then CCC, RMSE and Jaccard in one pass (HapsEncodedMatrix::calcSharedTargetMeasures())
	 *
	 * @param haps the encoded haplotypes
	 * @param binSize the number of sample pairs each thread takes at a time
	 * @param verbose print progress
	 */
	[[nodiscard]] static Measures calcMeasures(HapsEncodedMatrix & haps, uint32_t binSize = 1000, bool verbose = false);

	/**@brief the distance matrix DBSCAN uses, 1 - CCC, RMSE or 1 - Jaccard, the diagonal 0
	 *
	 */
	[[nodiscard]] static std::vector<std::vector<double>> distances(const Measures & measures, Measure measure, uint32_t numThreads = 1);

	struct Edge {
		uint32_t sample1_;
		uint32_t sample2_;
		uint32_t group_; /**< the group of sample1_*/
	};

	static constexpr uint32_t noGroup = std::numeric_limits<uint32_t>::max();

	struct Results {
		std::vector<uint32_t> groups_; /**< group by sample index, noGroup when it's in none*/
		std::vector<bool> lowCoverage_; /**< by sample index, below the min loci coverage so not clustered*/
		std::map<uint32_t, std::vector<uint32_t>> groupIndexes_; /**< sample indexes by group, noGroup holds the rest including low coverage*/
		std::map<uint32_t, std::map<std::string, double>> cccStats_; /**< stats (getStatsOnVec) on the pairs within each group, noGroup left out*/
		std::map<uint32_t, std::map<std::string, double>> rmseStats_;
		std::map<uint32_t, std::map<std::string, double>> jaccardStats_;
		std::vector<Edge> edges_; /**< the graph edges still on after clustering*/
		uint32_t numberOfGroups_{0};
	};

	/**@brief DBSCAN on the pars measure, then small groups of off nodes (when min neighbors is above 2)
	 *
	 * @param haps the encoded haplotypes measures came from
	 * @param measures from calcMeasures()
	 * @param pars the measure, cut off and DBSCAN settings
	 */
	[[nodiscard]] static Results cluster(const HapsEncodedMatrix & haps, const Measures & measures, const Pars & pars);

	/**@brief order of items from average linkage clustering, built with the nearest neighbor chain, items in the same cluster next to each other
	 *
	 * @param items indexes into distances
	 * @param distances the full distance matrix
	 * @param missingDistance used in place of NaN distances
	 * @return items reordered
	 */
	[[nodiscard]] static std::vector<uint32_t> averageLinkageOrder(const std::vector<uint32_t> & items,
			const std::vector<std::vector<double>> & distances, double missingDistance);
};

}  // namespace njhseq
