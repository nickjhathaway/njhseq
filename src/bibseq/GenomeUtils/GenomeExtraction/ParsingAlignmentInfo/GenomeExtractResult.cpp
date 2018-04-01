/*
 * GenomeExtractResult.cpp
 *
 *  Created on: May 3, 2017
 *      Author: nick
 */



#include "GenomeExtractResult.hpp"


namespace bibseq {



GenomeExtractResult::GenomeExtractResult(const std::shared_ptr<AlignmentResults> & ext,
		const std::shared_ptr<AlignmentResults> & lig) :
		ext_(ext), lig_(lig) {

}


void GenomeExtractResult::setRegion() {
	if (ext_->gRegion_.chrom_ != lig_->gRegion_.chrom_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error extention chrom, "
				<< ext_->gRegion_.chrom_ << "doesn't equal ligation chrom "
				<< lig_->gRegion_.chrom_ << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (ext_->gRegion_.reverseSrand_ == lig_->gRegion_.reverseSrand_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error extention and ligation are on the same strand, should be mapping to opposite strands"
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	if (ext_->gRegion_.reverseSrand_) {
		if (ext_->gRegion_.start_ < lig_->gRegion_.start_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error if extention is mapping to the reverse strand, it's start, "
					<< ext_->gRegion_.start_
					<< ", should be greater than ligation start, "
					<< lig_->gRegion_.start_ << "\n";
			throw std::runtime_error { ss.str() };
		}
	}else{
		if (ext_->gRegion_.start_ > lig_->gRegion_.start_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error if extention is mapping to the plus strand, it's start, "
					<< ext_->gRegion_.start_
					<< ", should be less than than ligation start, "
					<< lig_->gRegion_.start_ << "\n";
			throw std::runtime_error { ss.str() };
		}
	}

	size_t start = ext_->gRegion_.start_;
	size_t end = lig_->gRegion_.end_;
	size_t innerStart = ext_->gRegion_.end_;
	size_t innereEnd = lig_->gRegion_.start_;
	if(ext_->gRegion_.reverseSrand_){
		start = lig_->gRegion_.start_;
		end = ext_->gRegion_.end_;
		innerStart = lig_->gRegion_.end_;
		innereEnd = ext_->gRegion_.start_;
	}

	gRegion_ = std::make_shared<GenomicRegion>(ext_->gRegion_.uid_ + "-" + lig_->gRegion_.uid_, ext_->gRegion_.chrom_, start, end, ext_->gRegion_.reverseSrand_);
	gRegionInner_ = std::make_shared<GenomicRegion>(ext_->gRegion_.uid_ + "-" + lig_->gRegion_.uid_, ext_->gRegion_.chrom_, innerStart, innereEnd, ext_->gRegion_.reverseSrand_);
}


std::vector<GenomeExtractResult> getPossibleGenomeExtracts(const std::vector<std::shared_ptr<AlignmentResults>> & alnResultsExt,
		const std::vector<std::shared_ptr<AlignmentResults>> & alnResultsLig, const size_t insertSizeCutOff){
	std::vector<GenomeExtractResult> ret;
	//same chrom, opposite strands, less than the insert size
	for (const auto & ext : alnResultsExt) {
		for (const auto & lig : alnResultsLig) {
			//need to be on the same chromosome
			//need to be on opposite strands (should both should be in 5'->3' direction
			//and they shouldn't overlap
			if (ext->gRegion_.chrom_ == lig->gRegion_.chrom_
					&& ext->gRegion_.reverseSrand_ != lig->gRegion_.reverseSrand_
					&& !ext->gRegion_.overlaps(lig->gRegion_)
					&& ext->gRegion_.start_ != lig->gRegion_.end_
					&& ext->gRegion_.end_ != lig->gRegion_.start_ ) {

				if(ext->gRegion_.reverseSrand_){
					if(ext->gRegion_.start_ > lig->gRegion_.start_){
						GenomeExtractResult extraction(ext, lig);
						extraction.setRegion();
						if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
							ret.emplace_back(extraction);
						}
					}
				}else{
					if(ext->gRegion_.start_ < lig->gRegion_.start_){
						GenomeExtractResult extraction(ext, lig);
						extraction.setRegion();
						if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
							ret.emplace_back(extraction);
						}
					}
				}
			}
		}
	}
	return ret;
}


}  // namespace bibseq


