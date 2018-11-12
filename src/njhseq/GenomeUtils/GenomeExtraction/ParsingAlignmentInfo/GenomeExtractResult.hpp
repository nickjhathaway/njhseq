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

namespace njhseq {


class GenomeExtractResult {
public:
	GenomeExtractResult(const std::shared_ptr<AlignmentResults> & ext,
			const std::shared_ptr<AlignmentResults> & lig);
	std::shared_ptr<AlignmentResults> ext_;
	std::shared_ptr<AlignmentResults> lig_;

	std::shared_ptr<GenomicRegion> gRegion_ ;//= nullptr; it's silly but this is due to eclipse for now, should be able to revert once cdt 9.3 is released june 28 2017
	std::shared_ptr<GenomicRegion> gRegionInner_ ;

	void setRegion();

};


std::vector<GenomeExtractResult> getPossibleGenomeExtracts(const std::vector<std::shared_ptr<AlignmentResults>> & alnResultsExt,
		const std::vector<std::shared_ptr<AlignmentResults>> & alnResultsLig, const size_t insertSizeCutOff = std::numeric_limits<size_t>::max());





}  // namespace njhseq





