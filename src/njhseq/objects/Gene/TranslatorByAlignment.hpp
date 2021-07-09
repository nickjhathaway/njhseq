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
#include "njhseq/BamToolsUtils.h"
#include "njhseq/BamToolsUtils/ReAlignedSeq.hpp"
#include "njhseq/programUtils/seqSetUp.hpp"



namespace njhseq {

class TranslatorByAlignment{
public:

	struct Codon{
		Codon(char aa, std::tuple<char, char, char> bases): aa_(aa), bases_(bases){

		}
		char aa_;
		std::tuple<char, char, char> bases_;
	};

	struct AAInfo {
		AAInfo(std::string transcriptName,
				uint32_t zeroBasedPos,
				Codon cod,
				bool knownMut) : transcriptName_(transcriptName),
				 zeroBasedPos_(zeroBasedPos), cod_(cod), knownMut_(
						knownMut) {
		}
		std::string transcriptName_;
		uint32_t zeroBasedPos_;
		Codon cod_;
		bool knownMut_;
	};

	struct RunPars {
		uint32_t occurrenceCutOff = 2;
		double lowVariantCutOff = 0.005;
		ReAlignedSeq::genRealignmentPars realnPars;

	};

	struct VariantsInfo {

		VariantsInfo(const Bed3RecordCore & region, const seqInfo & refSeq);
		Bed3RecordCore region_;
		seqInfo seqBase_;


		std::map<uint32_t, std::map<char, uint32_t>> allBases;
		std::map<uint32_t, uint32_t> depthPerPosition;
		std::map<uint32_t, std::unordered_set<std::string>> samplesPerPosition;

		std::map<uint32_t, std::map<char, uint32_t>> snps;
		std::map<uint32_t, std::map<std::string,uint32_t>> insertions;
		std::map<uint32_t, std::map<std::string,uint32_t>> deletions;

		std::map<uint32_t, std::map<char, uint32_t>> snpsFinal;
		std::map<uint32_t, std::map<std::string,uint32_t>> insertionsFinal;
		std::map<uint32_t, std::map<std::string,uint32_t>> deletionsFinal;

		std::set<uint32_t> variablePositons_;

		Bed3RecordCore getVariableRegion();


		void addVariantInfo(const std::string & alignedRefSeq,
				const std::string & alignedQuerySeq,
				uint32_t querySeqCount,
				const std::unordered_set<std::string> & samples,
				const comparison & comp,
				uint32_t offSetStart);


		void setFinals(const RunPars & rPars);
		//void setFinals(const RunPars & rPars, uint32_t totalPopCount);

		char getBaseForGenomicRegionNoCheck(const uint32_t pos) const;

		char getBaseForGenomicRegion(const uint32_t pos) const;

		void writeVCF(const OutOptions &vcfOutOpts) const;

		void writeSNPTable(const OutOptions &snpTabOutOpts) const;

		uint32_t getFinalNumberOfSegratingSites() const;

		/**@brief write out info for positions, will throw if no info present
		 *
		 * The order written, name, position, ref base, query base, count, frequency, allele Depth, samples depth
		 *
		 * @param out the stream to write to
		 * @param name a name to output with the info
		 * @param snpPositions the positions
		 * @param oneBased whether to write out the positions as one base or zero based positions
		 */
		void writeOutSNPsInfo(std::ostream & out,
				const std::string & name,
				const std::set<uint32_t> & snpPositions,
				bool oneBased = false);
		/**@brief write out info determined snps
		 *
		 * The order written, name, position, ref base, query base, count, frequency, allele Depth, samples depth
		 *
		 * @param out the stream to write to
		 * @param name a name to output with the info
		 * @param oneBased whether to write out the positions as one base or zero based positions
		 */
		void writeOutSNPsFinalInfo(std::ostream & out,
				const std::string & name,
				bool oneBased = false);
		/**@brief write out info all positions
		 *
		 * The order written, name, position, ref base, query base, count, frequency, allele Depth, samples depth
		 *
		 * @param out the stream to write to
		 * @param name a name to output with the info
		 * @param oneBased whether to write out the positions as one base or zero based positions
		 */
		void writeOutSNPsAllInfo(std::ostream & out,
				const std::string & name,
				bool oneBased = false);

		static VecStr SNPHeaderAminoAcid(){
			return VecStr{"transcript","position(1-based)","refAA","AA","count","freq","alleleDepth","samples"};
		}
		static VecStr SNPHeaderGenomic(){
			return VecStr{"transcript","position(0-based)","ref","base","count","freq","alleleDepth","samples"};
		}
	};

	struct TranslatorByAlignmentPars{
		TranslatorByAlignmentPars();
		bfs::path gffFnp_ = "";
		BioCmdsUtils::LastZPars lzPars_;
		bool useLastz_ = false;
		std::string additionalBowtieArguments_;
		bfs::path workingDirtory_;

		bool keepTemporaryFiles_ = false;

		bfs::path knownAminoAcidMutationsFnp_; //! at least 2 columns, transcriptid, aaposition/

		void setOptions(seqSetUp & setUp);
	};

	struct TranslateSeqRes {
		std::string transcriptName_;

		seqInfo cDna_;
		seqInfo translation_;
		seqInfo queryAlnTranslation_;
		seqInfo refAlnTranslation_;

		std::tuple<GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo> firstAminoInfo_;
		std::tuple<GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo> lastAminoInfo_;


		comparison comp_;

		Codon getCodonForAARefPos(uint32_t aaRefPos) const;
	};



	struct TranslatorByAlignmentResult{
		std::set<std::string> geneIds_;
		std::unordered_map<std::string, std::unordered_map<std::string, TranslateSeqRes>> translations_;

		std::unordered_map<std::string, std::vector<ReAlignedSeq>> seqAlns_;

		//by transcript name
		std::unordered_map<std::string, VariantsInfo> proteinVariants_;
		std::unordered_map<std::string, std::string> proteinForTranscript_;
		//by chromosome
		std::unordered_map<std::string, VariantsInfo> seqVariants_;
		std::unordered_map<std::string, std::unordered_map<uint32_t, char>> baseForPosition_;

		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>>  transcriptInfosForGene_;
		std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>  translationInfoForTranscirpt_;


		VecStr seqsUnableToBeMapped_;

		void writeSeqLocations(std::ostream & out) const;
	};




	static std::unordered_map<std::string, TranslateSeqRes> translateBasedOnAlignment(
			const BamTools::BamAlignment & bAln,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			TwoBit::TwoBitFile & tReader,
			aligner & alignerObj,
			const BamTools::RefVector & refData);

	static std::unordered_map<std::string, TranslateSeqRes> translateBasedOnAlignment(
			const ReAlignedSeq & realigned,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			aligner & alignerObj);


	TranslatorByAlignmentPars pars_;

	std::unordered_map<std::string, std::vector<uint32_t>> knownAminoAcidPositions_;
	std::unordered_map<std::string, std::unordered_map<uint32_t, MetaDataInName>> metaDataAssociatedWithAminoacidPosition_;

	TranslatorByAlignment(const TranslatorByAlignmentPars & pars);

	TranslatorByAlignmentResult run(const SeqIOOptions & seqOpts,
			const std::unordered_map<std::string, std::unordered_set<std::string>> & sampCountsForHaps,
			const RunPars & rPars);


	static std::unordered_map<std::string, std::set<uint32_t>> readInAAPositions(const bfs::path & knownAminoAcidChangesFnp);

};

}  // namespace njhseq



