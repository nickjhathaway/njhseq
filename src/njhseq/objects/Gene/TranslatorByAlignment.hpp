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

class VCFOutput {
public:
	std::string vcfFormatVersion_{"VCFv4.0"};
	struct FormatOrInfoEntry {
		FormatOrInfoEntry() = default;
		FormatOrInfoEntry(		std::string id,
		std::string number,
		std::string type,
		std::string description): id_(std::move(id)),
			number_(std::move(number)),
			type_(std::move(type)),
			description_(std::move(description)) {

		}
		std::string id_;
		std::string number_;
		std::string type_;
		std::string description_;
	};

	std::vector<FormatOrInfoEntry> formatEntries_;
	std::vector<FormatOrInfoEntry> infoEntries_;

	class VCFRecord {
	public:
		VCFRecord() = default;
		std::string chrom_;/**<chromosome name */
		uint32_t pos_{std::numeric_limits<uint32_t>::max()};/**<position, 1-based to conform with VCF standards */
		std::string id_{"."};/**<id for this variant */
		std::string ref_;/**<reference sequence for this variant */
		VecStr alts_;/**<all the alternative alleles */
		uint32_t qual_{40};/**<the quality for this variant */
		std::string filter_{"PASS"};/**<whether or not this variant passes QC checks */

		MetaDataInName info_;/**<the info entry for this variant */
		std::map<std::string, MetaDataInName> sampleFormatInfos_;/**<the sample info for this variant */
	};

	std::vector<VCFRecord> records;

	void sortRecords();

	/**
	 * \brief write out the header and the CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO fields
	 * \param out output stream
	 */
	void writeOutFixedOnly(std::ostream & out);

	/**
	 * \brief write out the header and the CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO fields and then also FORMAT column and the sample columns after that
	 * \param out the output stream to write to
	 */
	void writeOutFixedAndSampleMeta(std::ostream & out);
};


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

		RunPars();
		uint32_t occurrenceCutOff = 2;
		double lowVariantCutOff = 0.005;
		ReAlignedSeq::genRealignmentPars realnPars;

	};

	struct VariantsInfo {

		VariantsInfo(const Bed3RecordCore & region, const seqInfo & refSeq);
		Bed3RecordCore region_;//!< the reference region being compared to
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

		VCFOutput createVCFOutputFixed() const;

		VCFOutput writeVCF(const OutOptions &vcfOutOpts) const;
		VCFOutput writeVCF(std::ostream & out) const;

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
		void writeOutSNPsInfo(
				std::ostream & out,
				const std::string & name,
				const std::set<uint32_t> & snpPositions,
				bool oneBased = false);
		void writeOutSNPsInfo(
				const OutOptions & outOpts,
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
		void writeOutSNPsFinalInfo(
				std::ostream & out,
				const std::string & name,
				bool oneBased = false);
		void writeOutSNPsFinalInfo(
				const OutOptions & outOpts,
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
		void writeOutSNPsAllInfo(
				std::ostream & out,
				const std::string & name,
				bool oneBased = false);
		void writeOutSNPsAllInfo(
				const OutOptions & outOpts,
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
		bool writeOutGeneInfos_ = true;

		bfs::path knownAminoAcidMutationsFnp_; //! at least 2 columns, transcriptid, aaposition/

		uint32_t aaExpand_ = 10;
		bool useFullProtein_ = false;

		uint32_t allowableStopCodons_ {1};

		void setOptions(seqSetUp & setUp, bool requireGenome = false);
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

		Bed6RecordCore genBedRec() const;
	};



	struct TranslatorByAlignmentResult{
		std::set<std::string> geneIds_;
		std::unordered_map<std::string, std::unordered_map<std::string, TranslateSeqRes>> translations_;
		std::unordered_map<std::string, std::unordered_map<std::string, TranslateSeqRes>> filteredOffTranslations_;

		std::unordered_map<std::string, std::vector<ReAlignedSeq>> seqAlns_;

		//by transcript name
		std::unordered_map<std::string, VariantsInfo> proteinVariants_;
		//by chromosome
		std::unordered_map<std::string, VariantsInfo> seqVariants_;

		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>>  transcriptInfosForGene_;
		std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>  translationInfoForTranscirpt_;

		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> fullAATypedWithCodonInfo_;
		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> variantAATypedWithCodonInfo_;

		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> translated_fullAATypedWithCodonInfo_;
		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> translated_variantAATypedWithCodonInfo_;

		std::map<std::string, std::string> genSeqSNPTypedStr() const;


		std::map<std::string, std::string> genAATypedStr() const;
		std::map<std::string, std::string> genAATypedStrOnlyKnowns() const;
		std::map<std::string, std::string> genAATypedStrOnlyPopVariant() const;

		std::map<std::string, std::string> translated_genAATypedStr() const;
		std::map<std::string, std::string> translated_genAATypedStrOnlyKnowns() const;
		std::map<std::string, std::string> translated_genAATypedStrOnlyPopVariant() const;

		VecStr seqsUnableToBeMapped_;
		VecStr seqsTranslationFiltered_;

		std::set<std::string> getAllSeqNames() const;

		void writeSeqLocations(std::ostream & out) const;
		void writeSeqLocationsTranslation(std::ostream & out) const;

		void writeOutSeqAlnIndvVars(const OutOptions & outopts) const;
		void writeOutTranslatedIndvVars(const OutOptions & outOpts, const std::unordered_map<std::string, std::set<uint32_t>> & knownAminoAcidPositions1Based = {}) const;

		void writeOutAATypedInfo(const OutOptions &outOpts) const;
	};




//	std::unordered_map<std::string, TranslateSeqRes> translateBasedOnAlignment(
//			const BamTools::BamAlignment & bAln,
//			const GeneFromGffs & currentGene,
//			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
//			TwoBit::TwoBitFile & tReader,
//			aligner & alignerObj,
//			const BamTools::RefVector & refData);

	std::unordered_map<std::string, TranslateSeqRes> translateBasedOnAlignment(
			const ReAlignedSeq & realigned,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			aligner & alignerObj);


	TranslatorByAlignmentPars pars_;

	std::unordered_map<std::string, std::set<uint32_t>> knownAminoAcidPositions_;
	std::unordered_map<std::string, std::unordered_map<uint32_t, MetaDataInName>> metaDataAssociatedWithAminoacidPosition_;

	TranslatorByAlignment(const TranslatorByAlignmentPars & pars);

	std::set<uint32_t> getAllInterestingAAPosZeroBased(const std::string & transcript, const TranslatorByAlignmentResult & results);

	TranslatorByAlignmentResult run(const SeqIOOptions & seqOpts,
			const std::unordered_map<std::string, std::unordered_set<std::string>> & sampCountsForHaps,
			const RunPars & rPars);


	static std::unordered_map<std::string, std::set<uint32_t>> readInAAPositions(const bfs::path & knownAminoAcidChangesFnp);

};

}  // namespace njhseq



