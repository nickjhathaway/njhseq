/*
 * GenomeRegionsGenerator.cpp
 *
 *  Created on: May 26, 2019
 *      Author: nicholashathaway
 */


#include "GenomeRegionsGenerator.hpp"

namespace njhseq {


GenomeRegionsGenerator::GenomeRegionsGenerator(const std::vector<GenomicRegion> & regions,
		const Params & pars) :
		regions_(regions), pars_(pars) {

}



bool GenomeRegionsGenerator::genRegionLockFree(GenomicRegion & reg){
	if(currentRegPos_ < regions_.size()){
		reg = GenomicRegion("",
				regions_[currentRegPos_].chrom_,
				regions_[currentRegPos_].start_ + currentPositionInReg_,
				std::min(regions_[currentRegPos_].start_ + currentPositionInReg_ + pars_.window_, regions_[currentRegPos_].end_),
				regions_[currentRegPos_].reverseSrand_);
		currentPositionInReg_ += pars_.step_;
		if(regions_[currentRegPos_].start_ + currentPositionInReg_ >= regions_[currentRegPos_].end_){
			++currentRegPos_;
			currentPositionInReg_ = 0;
		}
		return true;
	}
	return false;
}

bool GenomeRegionsGenerator::genRegion(GenomicRegion & reg){
	std::lock_guard<std::mutex> lock(mut_);
	return genRegionLockFree(reg);
}

std::vector<GenomicRegion> GenomeRegionsGenerator::genRegions(uint32_t numberOfRegions){
	std::vector<GenomicRegion> ret;
	std::lock_guard<std::mutex> lock(mut_);
	GenomicRegion reg;
	while(genRegionLockFree(reg)){
		ret.emplace_back(reg);
	}
	return ret;
}



}  // namespace njhseq
