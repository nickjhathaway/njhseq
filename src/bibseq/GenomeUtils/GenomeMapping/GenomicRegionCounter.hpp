#pragma once
/*
 * GenomicRegionCounter.hpp
 *
 *  Created on: Mar 31, 2017
 *      Author: nick
 */

#include "bibseq/utils.h"
#include "bibseq/objects/BioDataObject/GenomicRegion.hpp"

namespace bibseq {
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

};

}  // namespace bibseq

