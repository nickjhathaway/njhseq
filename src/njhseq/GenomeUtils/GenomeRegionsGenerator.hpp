#pragma once

/*
 * GenomeRegionsGenerator.hpp
 *
 *  Created on: May 26, 2019
 *      Author: nicholashathaway
 */



#include "njhseq/common.h"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"
#include "njhseq/GenomeUtils/GenomeUtils.hpp"

namespace njhseq {


class GenomeRegionsGenerator {
public:
	struct Params {
		uint32_t step_{100};
		uint32_t window_{100};
	};

	GenomeRegionsGenerator(const std::vector<GenomicRegion> & regions,
			const Params & pars);

	std::vector<GenomicRegion> regions_;
	Params pars_;

	uint32_t currentRegPos_{0};
	uint32_t currentPositionInReg_{0};
	std::mutex mut_;

	bool genRegionLockFree(GenomicRegion & reg);
	bool genRegion(GenomicRegion & reg);
	std::vector<GenomicRegion> genRegions(uint32_t numberOfRegions);

};





}  // namespace njhseq


