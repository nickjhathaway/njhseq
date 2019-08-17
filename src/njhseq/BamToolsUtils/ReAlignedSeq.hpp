#pragma once

/*
 * ReAlignedSeq.hpp
 *
 *  Created on: Aug 16, 2019
 *      Author: nicholashathaway
 */


#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"
#include "njhseq/alignment/aligner.h"



namespace njhseq {

class ReAlignedSeq{
public:

	BamTools::BamAlignment bAln_;

	GenomicRegion gRegion_;

	seqInfo querySeq_;
	seqInfo refSeq_;
	seqInfo alnQuerySeq_;
	seqInfo alnRefSeq_;


	comparison comp_;

	struct genRealignmentPars{
		bool adjustForSoftClipping = true;
		uint32_t extendAmount = 5;
	};


	static ReAlignedSeq genRealignment(const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData,
			aligner & alignerObj,
			const std::unordered_map<std::string, uint32_t> & chromLengths,
			TwoBit::TwoBitFile & tReader,
		const genRealignmentPars & pars);
};





}  // namespace njhseq



