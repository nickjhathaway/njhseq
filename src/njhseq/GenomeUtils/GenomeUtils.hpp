#pragma once

/*
 * GenomeUtils.hpp
 *
 *  Created on: May 25, 2019
 *      Author: nicholashathaway
 */

#include <api/BamReader.h>
#include "njhseq/common.h"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"


namespace njhseq {

std::vector<GenomicRegion> genChromosomeGenRegions(const BamTools::RefVector & refData);

std::vector<GenomicRegion> genChromosomeGenRegions(const std::unordered_map<std::string, uint32_t> & refData);

}  // namespace njhseq



