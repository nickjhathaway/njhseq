#pragma once

/*
 * kMeans.hpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "njhseq/common.h"

namespace njhseq {

/**@brief k-means clustering of numeric points with squared Euclidean distance, and farthest point picking for starting centers
 *
 * Every tie goes to the lowest index so the results are always the same; seqServer's utils.js has the same functions (njhKMeans*) for the viewers.
 */
class KMeans {
public:
	using Point = std::vector<double>;

	struct Pars {
		uint32_t maxIterations_ = 50;
		/**@brief move each center to the closest point not already taken by an earlier center, so centers are always actual (distinct) points
		 *
		 */
		bool snapToPoints_ = false;
	};

	struct Results {
		std::vector<Point> centers_;
		std::vector<uint32_t> centerPoints_; /**< with snapToPoints_, the point each center is*/
		std::vector<uint32_t> assignments_; /**< the center of each point*/
		uint32_t iterations_{0}; /**< rounds of assigning points*/
		bool converged_{false}; /**< the assignments stopped changing before maxIterations_*/
	};

	/**@brief squared Euclidean distance
	 *
	 */
	[[nodiscard]] static double dist2(const Point & first, const Point & second);

	/**@brief the mean of the points, in order
	 *
	 */
	[[nodiscard]] static Point mean(const std::vector<Point> & points);

	/**@brief k points spread out: the first is the point farthest from start, each next one is the point farthest from those already picked
	 *
	 * With k as the number of points this orders all of them, each as far as it can be from the ones before it.
	 *
	 * @param points the points
	 * @param k how many to pick, at most the number of points
	 * @param start where to measure the first pick from, e.g. the mean
	 * @return the positions of the picked points, in pick order
	 */
	[[nodiscard]] static std::vector<uint32_t> farthestPoints(const std::vector<Point> & points, uint32_t k, const Point & start);

	/**@brief Lloyd's k-means from starting centers
	 *
	 * Each round assigns points to their closest center, stops when no assignment changes, otherwise moves each center to the mean of its points
	 * (a center with no points stays put), snapped to a point when Pars::snapToPoints_ is set.
	 *
	 * @param points the points
	 * @param startPoints the positions of the points to start the centers at
	 * @param pars rounds and snapping
	 */
	[[nodiscard]] static Results run(const std::vector<Point> & points, const std::vector<uint32_t> & startPoints, const Pars & pars);

	/**@brief run() with the default Pars
	 *
	 */
	[[nodiscard]] static Results run(const std::vector<Point> & points, const std::vector<uint32_t> & startPoints);
};

}  // namespace njhseq
