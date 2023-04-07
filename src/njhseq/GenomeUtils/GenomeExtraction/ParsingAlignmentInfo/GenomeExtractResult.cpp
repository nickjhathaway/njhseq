/*
 * GenomeExtractResult.cpp
 *
 *  Created on: May 3, 2017
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


#include "GenomeExtractResult.hpp"

#include <utility>


namespace njhseq {


GenomeExtractResult::GenomeExtractResult(const std::shared_ptr<AlignmentResults> &ext,
                                         const std::shared_ptr<AlignmentResults> &lig) : GenomeExtractResult(
    *ext, *lig) {

}

GenomeExtractResult::GenomeExtractResult(const std::shared_ptr<ReAlignedSeq> &ext,
                                         const std::shared_ptr<ReAlignedSeq> &lig) : GenomeExtractResult(
    *ext, *lig) {

}

GenomeExtractResult::GenomeExtractResult(const AlignmentResults &ext,
                                         const AlignmentResults &lig) : GenomeExtractResult(
    ext.gRegion_, ext.comp_, lig.gRegion_, lig.comp_) {

}

GenomeExtractResult::GenomeExtractResult(const ReAlignedSeq &ext,
                                         const ReAlignedSeq &lig) : GenomeExtractResult(
    ext.gRegion_, ext.comp_, lig.gRegion_, lig.comp_) {

}


GenomeExtractResult::GenomeExtractResult(GenomicRegion ext,
                                         comparison extComp,
                                         GenomicRegion lig,
                                         comparison ligComp) :
    extRegion_(std::move(ext)), extComp_(std::move(extComp)),
    ligRegion_(std::move(lig)), ligComp_(std::move(ligComp)) {

}




void GenomeExtractResult::setRegion() {
	if (extRegion_.chrom_ != ligRegion_.chrom_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error extension chrom, "
				<< extRegion_.chrom_ << "doesn't equal ligation chrom "
				<< ligRegion_.chrom_ << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (extRegion_.reverseSrand_ == ligRegion_.reverseSrand_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error extension and ligation are on the same strand, should be mapping to opposite strands"
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	if (extRegion_.reverseSrand_) {
		if (extRegion_.start_ < ligRegion_.start_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error if extension is mapping to the reverse strand, it's start, "
					<< extRegion_.start_
					<< ", should be greater than ligation start, "
					<< ligRegion_.start_ << "\n";
			throw std::runtime_error { ss.str() };
		}
	}else{
		if (extRegion_.start_ > ligRegion_.start_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error if extension is mapping to the plus strand, it's start, "
					<< extRegion_.start_
					<< ", should be less than than ligation start, "
					<< ligRegion_.start_ << "\n";
			throw std::runtime_error { ss.str() };
		}
	}

	size_t start = extRegion_.start_;
	size_t end = ligRegion_.end_;
	size_t innerStart = extRegion_.end_;
	size_t innereEnd = ligRegion_.start_;
	if(extRegion_.reverseSrand_){
		start = ligRegion_.start_;
		end = extRegion_.end_;
		innerStart = ligRegion_.end_;
		innereEnd = extRegion_.start_;
	}

	gRegion_ = std::make_shared<GenomicRegion>(extRegion_.uid_ + "-" + ligRegion_.uid_, extRegion_.chrom_, start, end, extRegion_.reverseSrand_);
	gRegionInner_ = std::make_shared<GenomicRegion>(extRegion_.uid_ + "-" + ligRegion_.uid_, extRegion_.chrom_, innerStart, innereEnd, extRegion_.reverseSrand_);
}





}  // namespace njhseq


