#pragma once
/*
 * GenomeExtractResult.hpp
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
#include "njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/AlignmentResults.hpp"
#include "njhseq/BamToolsUtils/ReAlignedSeq.hpp"

namespace njhseq {


class GenomeExtractResult {
public:
	GenomeExtractResult(GenomicRegion extRegion,
      comparison extComp,
			GenomicRegion ligRegion,
      comparison ligComp);
  GenomeExtractResult(const std::shared_ptr<AlignmentResults> & ext,
                      const std::shared_ptr<AlignmentResults> & lig);
  GenomeExtractResult(const std::shared_ptr<ReAlignedSeq> & ext,
                      const std::shared_ptr<ReAlignedSeq> & lig);

  GenomeExtractResult(const AlignmentResults & ext,
                      const AlignmentResults & lig);
  GenomeExtractResult(const ReAlignedSeq & ext,
                      const ReAlignedSeq & lig);
  GenomicRegion extRegion_;
  comparison extComp_;

  GenomicRegion ligRegion_;
  comparison ligComp_;


  std::shared_ptr<GenomicRegion> gRegion_ ;//= nullptr; it's silly but this is due to eclipse for now, should be able to revert once cdt 9.3 is released june 28 2017
	std::shared_ptr<GenomicRegion> gRegionInner_ ;

	void setRegion();

};


template<typename ALN_RESULTS>
std::vector<GenomeExtractResult> getPossibleGenomeExtracts(const std::vector<ALN_RESULTS> & alnResultsExt,
		const std::vector<ALN_RESULTS> & alnResultsLig, const size_t insertSizeCutOff = std::numeric_limits<size_t>::max()){
  std::vector<GenomeExtractResult> ret;
  //same chrom, opposite strands, less than the insert size
  for (const auto & ext : alnResultsExt) {
    for (const auto & lig : alnResultsLig) {
      //need to be on the same chromosome
      //need to be on opposite strands (should both should be in 5'->3' direction
      //and they shouldn't overlap
      if (getRef(ext).gRegion_.chrom_ == getRef(lig).gRegion_.chrom_
          && getRef(ext).gRegion_.reverseSrand_ != getRef(lig).gRegion_.reverseSrand_
          && !getRef(ext).gRegion_.overlaps(getRef(lig).gRegion_)
          && getRef(ext).gRegion_.start_ != getRef(lig).gRegion_.end_
          && getRef(ext).gRegion_.end_ != getRef(lig).gRegion_.start_ ) {

        if(getRef(ext).gRegion_.reverseSrand_){
          if(getRef(ext).gRegion_.start_ > getRef(lig).gRegion_.start_){
            GenomeExtractResult extraction(ext, lig);
            extraction.setRegion();
            if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
              ret.emplace_back(extraction);
            }
          }
        }else{
          if(getRef(ext).gRegion_.start_ < getRef(lig).gRegion_.start_){
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

}  // namespace njhseq





