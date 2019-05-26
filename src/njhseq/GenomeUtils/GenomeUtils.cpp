/*
 * GenomeUtils.cpp
 *
 *  Created on: May 25, 2019
 *      Author: nicholashathaway
 */



#include "GenomeUtils.hpp"
namespace njhseq {

std::vector<GenomicRegion> genChromosomeGenRegions(const BamTools::RefVector & refData){
	std::vector<GenomicRegion> ret;
	for(const auto & ref : refData){
		ret.emplace_back(GenomicRegion(ref.RefName,ref.RefName, 0, ref.RefLength, false));
	}
	return ret;
}

std::vector<GenomicRegion> genChromosomeGenRegions(const std::unordered_map<std::string, uint32_t> & refData){
	std::vector<GenomicRegion> ret;
	for(const auto & ref : refData){
		ret.emplace_back(GenomicRegion(ref.first,ref.first, 0, ref.second, false));
	}
	return ret;
}


}  // namespace njhseq
