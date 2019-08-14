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


	struct VariantsInfo {
		std::map<uint32_t, std::map<char, uint32_t>> allBases;

		std::map<uint32_t, std::map<char, uint32_t>> snps;
		std::map<uint32_t, std::map<std::string,uint32_t>> insertions;
		std::map<uint32_t, std::map<std::string,uint32_t>> deletions;

		std::map<uint32_t, std::map<char, uint32_t>> snpsFinal;
		std::map<uint32_t, std::map<std::string,uint32_t>> insertionsFinal;
		std::map<uint32_t, std::map<std::string,uint32_t>> deletionsFinal;

		std::set<uint32_t> variablePositons_;
	};

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

		//by transcript name
		std::unordered_map<std::string, VariantsInfo> proteinVariants_;
		std::unordered_map<std::string, std::string> proteinForTranscript_;
		//by chromosome
		std::unordered_map<std::string, VariantsInfo> seqVariants_;
		std::unordered_map<std::string, std::unordered_map<uint32_t, char>> baseForPosition_;
	};

	struct RunPars {
		uint32_t occurrenceCutOff = 2;
		double lowVariantCutOff = 0.005;
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

	TranslatorByAlignmentResult run(const SeqIOOptions & seqOpts,
			const std::unordered_map<std::string, uint32_t> & sampCountsForHaps,
			const RunPars & rPars);

};

}  // namespace njhseq



