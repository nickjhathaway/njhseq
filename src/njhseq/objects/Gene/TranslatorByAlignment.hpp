#pragma once

/*
 * TranslatorByAlignment.hpp
 *
 *  Created on: Aug 14, 2019
 *      Author: nicholashathaway
 */





#include <utility>

#include "njhseq/objects/BioDataObject.h"
#include "njhseq/GenomeUtils.h"
#include "njhseq/objects/Gene/GeneFromGffs.hpp"
#include "njhseq/BamToolsUtils.h"
#include "njhseq/BamToolsUtils/ReAlignedSeq.hpp"
#include "njhseq/programUtils/seqSetUp.hpp"
#include "njhseq/objects/Gene/VCFOutput.hpp"



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
		       bool knownMut,
		       Codon refCod) : transcriptName_(std::move(transcriptName)),
		                       zeroBasedPos_(zeroBasedPos), cod_(std::move(cod)), refCod_(std::move(refCod)),
		                       knownMut_(
			                       knownMut) {
		}
		std::string transcriptName_;
		uint32_t zeroBasedPos_;
		Codon cod_;
		Codon refCod_;
		bool knownMut_;
	};

	struct getComplexPositionsPars{
		uint32_t withinDist = 15;
		double fractionOfCoveredSamples = 0.95;/**<fraction of samples that must cover a position for it to be linked into a complex variant */
		double fractionOfSamplesCoveringSpan = 0.5;/**<fraction of samples that must cover a complex variant's full span to report it */
	};

	struct RunPars {

		RunPars();
		uint32_t occurrenceCutOff = 2;/**<The cutoff for the number of samples to appear in */
		double lowVariantCutOff = 0.005;/**<The allele frequency cut off(not to be confused with sample prevalence) */
		uint32_t totalReadDepthCutOff = 0;/**<The total read depth summed across all samples the haplotype appears in */

		// uint32_t withinDistComplexVariant = 15;/**<the distance within to collapse complex variants */
		getComplexPositionsPars complexVarPars;
		ReAlignedSeq::genRealignmentPars realnPars;
		uint32_t ploidy{2};/**<The set ploidy to use for vcf output */
	};

	struct VariantsInfo {

		VariantsInfo(const Bed3RecordCore & region, seqInfo  refSeq);
		Bed3RecordCore region_;//!< the reference region being compared to
		seqInfo seqBase_;

		std::set<uint32_t> alwaysReportLocations_{};//!< always report these positions no matter the count
		std::unordered_map<uint32_t, std::vector<std::string>> forcedAltCalls_{};//!< always report these alts if not detected, will default to zero

		std::map<uint32_t, std::map<char, uint32_t>> allBases;//!< the allele depth per base at position
		std::map<uint32_t, uint32_t> depthPerPosition;//!< the allele depth per position
		std::map<uint32_t, std::unordered_map<std::string, uint32_t>> samplesPerPosition;//!< the samples depth per position
		struct posAlleleCountSamples{
			posAlleleCountSamples() = default;
			uint32_t alleleCount_ = 0;

			std::unordered_map<std::string, uint32_t> sampleReadCnts_;
			[[nodiscard]] uint32_t getTotalReadDepth() const;
			[[nodiscard]] Json::Value toJson() const;
		};

		std::map<uint32_t, std::map<char, posAlleleCountSamples>> snps;
		std::map<uint32_t, std::map<std::string, posAlleleCountSamples>> insertions;
		std::map<uint32_t, std::map<std::string, posAlleleCountSamples>> deletions;
		std::map<uint32_t, posAlleleCountSamples> refForSnps;

		std::map<uint32_t, std::map<char, posAlleleCountSamples>> snpsFinal;
		std::map<uint32_t, std::map<std::string, posAlleleCountSamples>> insertionsFinal;
		std::map<uint32_t, std::map<std::string, posAlleleCountSamples>> deletionsFinal;

		//for complex need to have both the reference and the alt
		//this will hold not only the complex variants but also the snps, indels so there can be one final collection to carry them all
		std::map<uint32_t, std::map<std::string, std::map<std::string, posAlleleCountSamples>>> complex;
		std::map<uint32_t, std::map<std::string, posAlleleCountSamples>> complexRefCount;

		//after re-applying filters after determining complex variants
		std::map<uint32_t, std::map<std::string, std::map<std::string, posAlleleCountSamples>>> complexFinal;


		std::set<uint32_t> variablePositons_;

		Bed3RecordCore getVariableRegion();
		//void setFinals(const RunPars & rPars, uint32_t totalPopCount);
		struct PosStartSize {
			PosStartSize() = default;

			PosStartSize(uint32_t start, uint32_t size): start_(start), size_(size) {
			}

			uint32_t start_{std::numeric_limits<uint32_t>::max()};
			uint32_t size_{std::numeric_limits<uint32_t>::max()};
			uint32_t count_{1};

			[[nodiscard]] uint32_t end() const {
				return start_ + size_;
			}

			[[nodiscard]] bool overlaps(const PosStartSize& other) const {
				return other.start_ < end() && start_ < other.end();
			}

			[[nodiscard]] bool withinDist(const PosStartSize& other, uint32_t dist) const {
				//overlap handles if there is any overlap between the two,
				//if not either the other start is within the distance to the start (other region is after this region) or
				//the end of the other region is within the start of this region (other region is before this region)
				return overlaps(other) || uAbsdiff(other.start_, end()) < dist || uAbsdiff(other.end(), start_) < dist;
			}

			static void sortComplexPositions(std::vector<PosStartSize>& complexPositions) {
				njh::sort(complexPositions, [](const PosStartSize& a, const PosStartSize& b) {
					if (a.start_ == b.start_) {
						return a.end() < b.end();
					}
					return a.start_ < b.start_;
				});
			}
		};

		std::vector<PosStartSize> getComplexPositions(const getComplexPositionsPars & pars);

		struct SpanCoverage {
			uint32_t depth_{0};/**< the haplotype count covering the whole span */
			std::set<std::string> samples_;/**< the samples with a haplotype covering the whole span */
		};
		/**
		 * \brief the coverage of a complex variant from the haplotypes that cover its whole span, both ref and alts (including alts filtered from complexFinal)
		 * \param pos the start of the complex variant
		 * \param ref the reference sequence of the complex variant
		 * \return the depth and samples covering the span
		 */
		[[nodiscard]] SpanCoverage getComplexSpanCoverage(uint32_t pos, const std::string & ref) const;


		void addVariantInfo(const std::string& alignedRefSeq,
		                    const std::string& alignedQuerySeq,
		                    uint32_t querySeqCount,
		                    const std::unordered_map<std::string, uint32_t>& sampleReadCnts,
		                    const comparison& comp,
		                    uint32_t offSetStart);

		void addComplexVariantInfo(const std::string& alignedRefSeq,
		                           const std::string& alignedQuerySeq,
		                           uint32_t querySeqCount,
		                           const std::unordered_map<std::string, uint32_t>& sampleReadCnts,
		                           const comparison& comp,
		                           uint32_t offSetStart,
		                           const std::vector<PosStartSize> & complexPositions);


		void setFinals(const RunPars& rPars);
		void setComplexFinals(const RunPars& rPars, const std::vector<PosStartSize> & complexPositions);


		[[nodiscard]] char getBaseForGenomicRegionNoCheck(const uint32_t pos) const;

		[[nodiscard]] char getBaseForGenomicRegion(const uint32_t pos) const;

		[[nodiscard]] VCFOutput createVCFOutputFixed() const;

		[[nodiscard]] VCFOutput createVCFOutputComplexFixed(const std::vector<PosStartSize>& complexPositions) const;

		[[nodiscard]] VCFOutput createVCFOutputComplexFixedWithSampleInfo(uint32_t ploidy, const std::string& chrom,
		                                                                  uint32_t chromLength,
		                                                                  const std::vector<PosStartSize>& complexPositions)
		const;


		[[nodiscard]] VCFOutput writeVCF(const OutOptions &vcfOutOpts) const;
		VCFOutput writeVCF(std::ostream & out) const;

		void writeSNPTable(const OutOptions &snpTabOutOpts) const;

		[[nodiscard]] uint32_t getFinalNumberOfSegregatingSites() const;
		[[nodiscard]] uint32_t getFinalNumberOfSegregatingSites(const std::set<std::string> & for_samples) const;

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

		/**
		 * @brief Write out the final snp info
		 * @param outOpts output options to write to
		 * @param name a name to output with the info
		 * @param oneBased whether to write out the positioons as one base or zero based positions
		 */
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

		[[nodiscard]] std::set<std::string> getAllSamples() const;
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

	  bool correct_small_homopolymer_errors_ = false;
	  uint32_t min_hp_run_size_to_correct_ = 8;
	  uint32_t max_hp_run_gap_size_ = 3;

		uint32_t allowableStopCodons_ {1};

		/**
		 * @brief Set the default options for translation by alignment
		 * @param setUp the sequence option setter object
		 * @param requireGenome whether or not the genome is required
		 */
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

		std::map<std::string, std::vector<AAInfo>> fullAATypedWithCodonInfo_;
		std::map<std::string, std::vector<AAInfo>> variantAATypedWithCodonInfo_;

		std::map<std::string, std::vector<AAInfo>> translated_fullAATypedWithCodonInfo_;
		std::map<std::string, std::vector<AAInfo>> translated_variantAATypedWithCodonInfo_;

		[[nodiscard]] std::map<std::string, std::string> genSeqSNPTypedStr() const;


		[[nodiscard]] std::map<std::string, std::string> genAATypedStr() const;
		[[nodiscard]] std::map<std::string, std::string> genAATypedStrOnlyKnowns() const;
		[[nodiscard]] std::map<std::string, std::string> genAATypedStrOnlyPopVariant() const;

		[[nodiscard]] std::map<std::string, std::string> translated_genAATypedStr() const;
		[[nodiscard]] std::map<std::string, std::string> translated_genAATypedStrOnlyKnowns() const;
		[[nodiscard]] std::map<std::string, std::string> translated_genAATypedStrOnlyPopVariant() const;

		VecStr seqsUnableToBeMapped_;
		VecStr seqsTranslationFiltered_;/**< seqs aligned to a coding gene that no gene could translate, sorted, set by setSeqsTranslationFiltered() */

		/**
		 * \brief set seqsTranslationFiltered_ to the seqs in it or in filteredOffTranslations_ that have no translation in translations_, run after all translations are added
		 */
		void setSeqsTranslationFiltered();

		[[nodiscard]] std::set<std::string> getAllSeqNames() const;

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

	static std::unordered_map<std::string, TranslateSeqRes> translateBasedOnAlignment(
			const ReAlignedSeq & realigned,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			aligner & alignerObj,
			const TranslatorByAlignmentPars & pars
			);


	TranslatorByAlignmentPars pars_;

	std::unordered_map<std::string, std::set<uint32_t>> knownAminoAcidPositions_;
	std::unordered_map<std::string, std::unordered_map<uint32_t, MetaDataInName>> metaDataAssociatedWithAminoacidPosition_;

	TranslatorByAlignment(TranslatorByAlignmentPars  pars);

	std::set<uint32_t> getAllInterestingAAPosZeroBased(const std::string & transcript, const TranslatorByAlignmentResult & results);

	/**
	 * \brief align sequences to the genome with an gff file and translate and get variant calls for protein and sequence
	 * \param seqOpts input sequence option to read seqs from
	 * \param sampReadCountsForHaps the sample read counts per seq in input file
	 * \param rPars the parameters for this run, has variant freq/occurrence cut offs
	 * \return the alignment/translation results
	 */
	TranslatorByAlignmentResult run(const SeqIOOptions & seqOpts,
	                                const std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> & sampReadCountsForHaps,
	                                const RunPars & rPars);

	/**
	* \brief align sequences to a speicific genomic region and translate and get variant calls for protein and sequence
	 * \param seqOpts input sequence option to read seqs from
	 * \param sampReadCountsForHaps the sample read counts per seq in input file
	 * \param refSeqRegion the genomic region to force alignment to
	 * \param rPars the parameters for this run, has variant freq/occurrence cut offs
	 * \return the alignment/translation results
	 */
	TranslatorByAlignmentResult run(const SeqIOOptions & seqOpts,
	                                const std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> & sampReadCountsForHaps,
	                                const GenomicRegion & refSeqRegion,
	                                const RunPars & rPars);

	template<typename T>
	TranslatorByAlignmentResult run(
		const std::vector<T>& seqs,
		//const std::vector<seqInfo> & seqs,
		const std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> & sampReadCountsForHaps,
		const GenomicRegion& refSeqRegion,
		const RunPars& rPars);

	static std::unordered_map<std::string, std::set<uint32_t>> readInAAPositions(const bfs::path & knownAminoAcidChangesFnp);


	struct GetGenomicLocationsForAminoAcidPositionsPars {
		bfs::path proteinMutantTypingFnp = "";
		bfs::path gffFnp = "";
		bfs::path twoBitFnp = "";
		bool zeroBased = false;
		bool collapsePerId = false;
		std::string addMetaField;
		OutOptions outOpts{bfs::path(""), ".bed"};
		bool doNotWrite = false;
	};

	struct GetGenomicLocationsForAminoAcidPositionsRet {
		std::vector<Bed6RecordCore> genomicLocs; //! zero-based positioning
		std::vector<Bed6RecordCore> transcriptLocs; //! zero-based positioning
	};

	static GetGenomicLocationsForAminoAcidPositionsRet getGenomicLocationsForAminoAcidPositions(
		const GetGenomicLocationsForAminoAcidPositionsPars& pars);

	/**
	 * \brief the min and max genomic positions covered by alignments on a chromosome
	 */
	struct MinMaxPos{
		MinMaxPos()= default;
		size_t minPos_{std::numeric_limits<uint32_t>::max()};
		size_t maxPos_{0};
	};

	/**
	 * \brief the known locations and the coding genes used by the run() functions
	 */
	struct RunSetUp {
		GetGenomicLocationsForAminoAcidPositionsRet locs;
		std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;/**< genes with no non-coding transcripts */
		uint64_t proteinMaxLen{0};
	};

	/**
	 * \brief the aligners used by the run() functions
	 */
	struct RunAligners {
		explicit RunAligners(uint64_t maxLen);
		aligner alignObj;
		aligner alignObjAdjusted;/**< gives more weight to matches to help with large deletions or insertions */
		aligner protein_alignObj;

		/**
		 * \brief get the aligner to realign a sequence with
		 * \param averageLen the average length of the input sequences
		 * \param len the length of the sequence being realigned
		 * \return alignObjAdjusted if len differs from averageLen by more than 100, otherwise alignObj
		 */
		aligner & forLength(uint32_t averageLen, uint32_t len);
	};

	/**
	 * \brief set ret.geneIds_, the transcript infos and the protein variant holders for the genes in ids, and write out the gene infos if requested
	 * \param ids the gene ids overlapping the input sequences
	 * \param seqMaxLen the max input sequence length, the minimum for RunSetUp::proteinMaxLen
	 * \param tReader reader for the genome
	 * \param twoBitFnp the path to the genome 2bit file
	 * \param ret the results to set up
	 * \return the known locations and the coding genes
	 */
	RunSetUp setUpGenesForRun(const std::set<std::string> & ids,
	                          uint64_t seqMaxLen,
	                          TwoBit::TwoBitFile & tReader,
	                          const std::string & twoBitFnp,
	                          TranslatorByAlignmentResult & ret) const;

	/**
	 * \brief get the realignment parameters for a sequence, extending the realignment by how much len differs from averageLen
	 * \param realnPars the starting realignment parameters
	 * \param averageLen the average length of the input sequences
	 * \param len the length of the sequence being realigned
	 * \return the realignment parameters to use
	 */
	static ReAlignedSeq::genRealignmentPars genRealignParsForLength(const ReAlignedSeq::genRealignmentPars & realnPars,
	                                                                uint32_t averageLen,
	                                                                uint32_t len);

	/**
	 * \brief translate an alignment with each gene in geneIds, adding to ret.translations_, ret.filteredOffTranslations_ or ret.seqsTranslationFiltered_, call ret.setSeqsTranslationFiltered() after all alignments
	 * \param seqName the name of the aligned sequence
	 * \param realigned the alignment to the genome
	 * \param geneIds the genes to translate with, must all be in setUp.genes
	 * \param setUp the run set up
	 * \param protein_alignObj aligner for the protein alignment
	 * \param ret the results to add to
	 */
	void addTranslationsForAlignment(const std::string & seqName,
	                                 const ReAlignedSeq & realigned,
	                                 const VecStr & geneIds,
	                                 const RunSetUp & setUp,
	                                 aligner & protein_alignObj,
	                                 TranslatorByAlignmentResult & ret) const;

	/**
	 * \brief call the sequence and protein variants and type the sequences, run after all alignments and translations are in ret
	 * \param minMaxPositionsPerChrom the regions covered by the alignments
	 * \param tReader reader for the genome
	 * \param locs the known locations to always report
	 * \param sampReadCountsForHaps the sample read counts per seq
	 * \param rPars the parameters for this run
	 * \param watch timer to log laps to
	 * \param ret the results to add to
	 */
	void callVariantsForRun(const std::unordered_map<std::string, MinMaxPos> & minMaxPositionsPerChrom,
	                        TwoBit::TwoBitFile & tReader,
	                        const GetGenomicLocationsForAminoAcidPositionsRet & locs,
	                        const std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> & sampReadCountsForHaps,
	                        const RunPars & rPars,
	                        njh::stopWatch & watch,
	                        TranslatorByAlignmentResult & ret);

	/**
	 * \brief call the protein variants and type the sequences by amino acid, run after all translations are in ret, called by callVariantsForRun()
	 * \param locs the known locations to always report
	 * \param sampReadCountsForHaps the sample read counts per seq
	 * \param rPars the parameters for this run
	 * \param watch timer to log laps to
	 * \param ret the results to add to
	 */
	void callProteinVariantsForRun(const GetGenomicLocationsForAminoAcidPositionsRet & locs,
	                               const std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> & sampReadCountsForHaps,
	                               const RunPars & rPars,
	                               njh::stopWatch & watch,
	                               TranslatorByAlignmentResult & ret);

};


template<typename T>
TranslatorByAlignment::TranslatorByAlignmentResult TranslatorByAlignment::run(
	const std::vector<T>& seqs,
	const std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> & sampReadCountsForHaps,
	const GenomicRegion& refSeqRegion,
	const RunPars& rPars) {
	njh::stopWatch watch;
	watch.setLapName("initial set up for variant calling");
	TranslatorByAlignmentResult ret;

	//get some length stats on input seqs
	uint64_t seqMaxLen = readVec::getMaxLength(seqs);
	seqMaxLen = seqMaxLen + rPars.realnPars.extendAmount * 2;
	auto averageLen = static_cast<uint32_t>(readVec::getAvgLength(seqs));

	{
		watch.startNewLap("set up for aligning");
		auto gprefix = bfs::path(pars_.lzPars_.genomeFnp).replace_extension("");
		auto twoBitFnp = gprefix.string() + ".2bit";
		TwoBit::TwoBitFile tReader(twoBitFnp);
		//get overlaping geneinfo
		auto ids = getFeatureIdsFromOverlappingRegions({refSeqRegion}, pars_.gffFnp_);
		const auto setUp = setUpGenesForRun(ids, seqMaxLen, tReader, twoBitFnp, ret);
		RunAligners aligners(setUp.proteinMaxLen);
		auto chromLengths = tReader.getSeqLens();

		//translate with every coding gene overlapping the region
		VecStr geneIdsToTranslate;
		for(const auto & g : ret.geneIds_){
			if(njh::in(g, setUp.genes)){
				geneIdsToTranslate.emplace_back(g);
			}
		}

		std::unordered_map<std::string, MinMaxPos> minMaxPositionsPerChrom;
		watch.startNewLap("aligning");
		for (const auto & seq : seqs) {
			const auto & seqName = getSeqBase(seq).name_;
			const auto seqLen = getSeqBase(seq).seq_.size();
			auto results = ReAlignedSeq::genRealignment(seq, refSeqRegion,
				aligners.forLength(averageLen, seqLen),
				chromLengths, tReader,
				genRealignParsForLength(rPars.realnPars, averageLen, seqLen));
			if(results.alnRefSeq_.seq_.front() != '-' && results.alnRefSeq_.seq_.back() != '-') {
				minMaxPositionsPerChrom[results.gRegion_.chrom_].minPos_ = std::min(minMaxPositionsPerChrom[results.gRegion_.chrom_].minPos_,results.gRegion_.start_);
				minMaxPositionsPerChrom[results.gRegion_.chrom_].maxPos_ = std::max(minMaxPositionsPerChrom[results.gRegion_.chrom_].maxPos_,results.gRegion_.end_);
				ret.seqAlns_[seqName].emplace_back(results);
				addTranslationsForAlignment(seqName, results, geneIdsToTranslate, setUp, aligners.protein_alignObj, ret);
			} else {
				ret.seqsUnableToBeMapped_.emplace_back(seqName);
			}
		}
		callVariantsForRun(minMaxPositionsPerChrom, tReader, setUp.locs, sampReadCountsForHaps, rPars, watch, ret);
	}

	OutOptions timeLogOpts(njh::files::make_path(pars_.workingDirtory_, "timeLog.txt"));
	timeLogOpts.overWriteFile_ = true;
	OutputStream timeLogOut(timeLogOpts);
	watch.logLapTimes(timeLogOut, true, 6, true);
	return ret;
}


}  // namespace njhseq



