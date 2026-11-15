/*
 * kMeans.cpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "kMeans.hpp"

namespace njhseq {

double KMeans::dist2(const Point & first, const Point & second) {
	double ret = 0;
	for (uint32_t dim = 0; dim < first.size(); ++dim) {
		const double diff = first[dim] - second[dim];
		ret += diff * diff;
	}
	return ret;
}

KMeans::Point KMeans::mean(const std::vector<Point> & points) {
	if (points.empty()) {
		return {};
	}
	Point ret(points.front().size(), 0);
	for (const auto & point : points) {
		for (uint32_t dim = 0; dim < ret.size(); ++dim) {
			ret[dim] += point[dim];
		}
	}
	for (auto & val : ret) {
		val /= points.size();
	}
	return ret;
}

std::vector<uint32_t> KMeans::farthestPoints(const std::vector<Point> & points, uint32_t k, const Point & start) {
	if (k > points.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "can't pick " << k << " from " << points.size() << " points" << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::vector<uint32_t> ret;
	std::vector<double> minDists(points.size());
	for (uint32_t pos = 0; pos < points.size(); ++pos) {
		minDists[pos] = dist2(points[pos], start);
	}
	while (ret.size() < k) {
		uint32_t next = 0;
		for (uint32_t pos = 1; pos < minDists.size(); ++pos) {
			if (minDists[pos] > minDists[next]) {
				next = pos;
			}
		}
		ret.emplace_back(next);
		for (uint32_t pos = 0; pos < points.size(); ++pos) {
			const double dist = dist2(points[pos], points[next]);
			if (1 == ret.size() || dist < minDists[pos]) {
				minDists[pos] = dist;
			}
		}
		//picked points are out of the running
		for (const auto & picked : ret) {
			minDists[picked] = -1;
		}
	}
	return ret;
}

KMeans::Results KMeans::run(const std::vector<Point> & points, const std::vector<uint32_t> & startPoints, const Pars & pars) {
	Results ret;
	if (points.empty() || startPoints.empty()) {
		return ret;
	}
	for (const auto & startPoint : startPoints) {
		ret.centers_.emplace_back(points[startPoint]);
	}
	ret.centerPoints_ = startPoints;
	ret.assignments_ = std::vector<uint32_t>(points.size(), std::numeric_limits<uint32_t>::max());
	const uint32_t dims = points.front().size();
	for (uint32_t iteration = 0; iteration < pars.maxIterations_; ++iteration) {
		++ret.iterations_;
		bool changed = false;
		for (uint32_t pos = 0; pos < points.size(); ++pos) {
			uint32_t best = 0;
			double bestDist = dist2(points[pos], ret.centers_[0]);
			for (uint32_t centerPos = 1; centerPos < ret.centers_.size(); ++centerPos) {
				const double dist = dist2(points[pos], ret.centers_[centerPos]);
				if (dist < bestDist) {
					bestDist = dist;
					best = centerPos;
				}
			}
			if (ret.assignments_[pos] != best) {
				ret.assignments_[pos] = best;
				changed = true;
			}
		}
		if (!changed) {
			ret.converged_ = true;
			break;
		}
		std::vector<Point> sums(ret.centers_.size(), Point(dims, 0));
		std::vector<uint32_t> counts(ret.centers_.size(), 0);
		for (uint32_t pos = 0; pos < points.size(); ++pos) {
			for (uint32_t dim = 0; dim < dims; ++dim) {
				sums[ret.assignments_[pos]][dim] += points[pos][dim];
			}
			++counts[ret.assignments_[pos]];
		}
		std::vector<bool> taken(points.size(), false);
		for (uint32_t centerPos = 0; centerPos < ret.centers_.size(); ++centerPos) {
			if (counts[centerPos] > 0) {
				for (uint32_t dim = 0; dim < dims; ++dim) {
					ret.centers_[centerPos][dim] = sums[centerPos][dim] / counts[centerPos];
				}
			}
			if (pars.snapToPoints_) {
				uint32_t best = std::numeric_limits<uint32_t>::max();
				double bestDist = std::numeric_limits<double>::max();
				for (uint32_t pos = 0; pos < points.size(); ++pos) {
					if (taken[pos]) {
						continue;
					}
					const double dist = dist2(points[pos], ret.centers_[centerPos]);
					if (dist < bestDist) {
						bestDist = dist;
						best = pos;
					}
				}
				taken[best] = true;
				ret.centerPoints_[centerPos] = best;
				ret.centers_[centerPos] = points[best];
			}
		}
	}
	return ret;
}

KMeans::Results KMeans::run(const std::vector<Point> & points, const std::vector<uint32_t> & startPoints) {
	return run(points, startPoints, Pars{});
}

}  // namespace njhseq
