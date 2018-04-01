#pragma once
/*
 * GenomeExtractResult.hpp
 *
 *  Created on: May 3, 2017
 *      Author: nick
 */

#include "bibseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/AlignmentResults.hpp"

namespace bibseq {


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





}  // namespace bibseq





