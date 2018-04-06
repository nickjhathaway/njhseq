/*
 * GenomicRegionCounter.cpp
 *
 *  Created on: Mar 31, 2017
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

#include "GenomicRegionCounter.hpp"

namespace bibseq {

GenomicRegionCounter::GenomicRegionCount::GenomicRegionCount(
		const GenomicRegion & region) :
		region_(region) {
}
GenomicRegionCounter::GenomicRegionCount::GenomicRegionCount(
		const GenomicRegion & region, uint32_t count) :
		region_(region), count_(count) {
}

bool GenomicRegionCounter::GenomicRegionCount::sameRegion(
		const GenomicRegion& region) const {
	return region_.sameRegion(region);
}
void GenomicRegionCounter::GenomicRegionCount::increaseCount(uint32_t count) {
	count_ += count;
}

bool GenomicRegionCounter::hasRegion(const GenomicRegion & region) const {
	auto search = counts_.find(region.createUidFromCoords());
	if (counts_.end() != search) {
		return true;
	}
	return false;
}

void GenomicRegionCounter::increaseCount(const GenomicRegion & region,
		uint32_t count) {
	if (hasRegion(region)) {
		counts_.at(region.createUidFromCoords()).increaseCount(count);
	} else {
		GenomicRegionCount counter(region, count);
		counts_.emplace(region.createUidFromCoords(), counter);
	}
}

std::vector<GenomicRegion> GenomicRegionCounter::getRegionsLargestOnTop() const {
	auto uids = getVectorOfMapKeys(counts_);
	bib::sort(uids,
			[this](const std::string & key1, const std::string & key2) {return counts_.at(key1).count_ > counts_.at(key2).count_;});
	std::vector<GenomicRegion> ret;
	for (const auto & key : uids) {
		ret.emplace_back(counts_.at(key).region_);
	}
	return ret;
}

void GenomicRegionCounter::increaseCounts(
		const std::vector<BamTools::BamAlignment> & bAlns,
		const BamTools::RefVector & refData) {
	for (const auto & bAln : bAlns) {
		GenomicRegion gRegion(bAln, refData);
		gRegion.setUidWtihCoords();
		increaseCount(gRegion, 1);
	}
}

}  // namespace bibseq

