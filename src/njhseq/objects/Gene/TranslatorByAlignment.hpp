#pragma once

/*
 * TranslatorByAlignment.hpp
 *
 *  Created on: Aug 14, 2019
 *      Author: nicholashathaway
 */



#include "njhseq/objects/BioDataObject.h"
#include "njhseq/GenomeUtils.h"
#include "njhseq/objects/Gene/GeneFromGffs.hpp"
#include "njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo.h"
#include "njhseq/BamToolsUtils.h"



namespace njhseq {

class TranslatorByAlignment{
public:

	struct TranslatorByAlignmentPars{
		TranslatorByAlignmentPars();
		bfs::path gffFnp_ = "";
		BioCmdsUtils::LastZPars lzPars_;
		bool useLastz_ = false;
		std::string additionalBowtieArguments_;
		bfs::path workingDirtory_;

		bool keepTemporaryFiles_ = false;
	};

	struct TranslateSeqRes {
		seqInfo translation_;
		seqInfo queryAlnTranslation_;
		seqInfo refAlnTranslation_;

		comparison comp_;
	};

	struct TranslatorByAlignmentResult{
		std::set<std::string> geneIds_;
		std::unordered_map<std::string, std::unordered_map<std::string, TranslateSeqRes>> translations_;

		std::unordered_map<std::string, std::vector<std::shared_ptr<AlignmentResults>>> seqAlns_;
	};


	static std::unordered_map<std::string, TranslateSeqRes> translateBasedOnAlignment(
			const BamTools::BamAlignment & bAln,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			TwoBit::TwoBitFile & tReader,
			aligner & alignerObj,
			const BamTools::RefVector & refData);


	TranslatorByAlignmentPars pars_;

	TranslatorByAlignment(const TranslatorByAlignmentPars & pars);

	TranslatorByAlignmentResult run(const SeqIOOptions & seqOpts);

};

}  // namespace njhseq



