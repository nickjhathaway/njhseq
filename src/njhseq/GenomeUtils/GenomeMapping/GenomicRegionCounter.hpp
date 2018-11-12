#pragma once
/*
 * GenomicRegionCounter.hpp
 *
 *  Created on: Mar 31, 2017
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
#include "njhseq/utils.h"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"

namespace njhseq {
class GenomicRegionCounter {
public:

	class GenomicRegionCount {
	public:
		GenomicRegionCount(const GenomicRegion & region);
		GenomicRegionCount(const GenomicRegion & region, uint32_t count);
		GenomicRegion region_;
		uint32_t count_ = 1;

		bool sameRegion(const GenomicRegion& region) const;
		void increaseCount(uint32_t count);
	};

	std::unordered_map<std::string, GenomicRegionCount> counts_;

	bool hasRegion(const GenomicRegion & region) const;

	void increaseCount(const GenomicRegion & region, uint32_t count = 1);

	std::vector<GenomicRegion> getRegionsLargestOnTop() const;

	void increaseCounts(const std::vector<BamTools::BamAlignment> & bAlns,
			const BamTools::RefVector & refData);

	std::set<std::string> getIntersectingGffIds(const bfs::path & gffFnp, const VecStr & features = {"gene"})const ;

	static GenomicRegionCounter countRegionsInBam(const bfs::path & bamFnp);
};





}  // namespace njhseq

