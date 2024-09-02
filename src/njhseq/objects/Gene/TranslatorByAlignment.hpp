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
		double fractionOfCoveredSamples = 0.95;
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
				return (other.end() > start_ && other.end() <= end()) || (other.start_ >= start_ && other.start_ < end());
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
		                           std::vector<PosStartSize> complexPositions);


		void setFinals(const RunPars& rPars);
		void setComplexFinals(const RunPars& rPars, std::vector<PosStartSize> complexPositions);


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

		[[nodiscard]] uint32_t getFinalNumberOfSegratingSites() const;

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
		VecStr seqsTranslationFiltered_;

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

	// TranslatorByAlignmentResult run(
	// 	const std::vector<seqInfo> & seqs,
	// 	const std::unordered_map<std::string, std::unordered_set<std::string>> & sampCountsForHaps,
	// 	const GenomicRegion & refSeqRegion,
	// 	const RunPars & rPars);

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

	struct GetGenomicLocationsForAminoAcidPositionsRet{
		std::vector<Bed6RecordCore> genomicLocs;
		std::vector<Bed6RecordCore> transcriptLocs;
	};

	static GetGenomicLocationsForAminoAcidPositionsRet getGenomicLocationsForAminoAcidPositions(const GetGenomicLocationsForAminoAcidPositionsPars & pars);

};


template<typename T>
TranslatorByAlignment::TranslatorByAlignmentResult TranslatorByAlignment::run(
	const std::vector<T>& seqs,
	//const std::vector<seqInfo> & seqs,
	const std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> & sampReadCountsForHaps,
	const GenomicRegion& refSeqRegion,
	const RunPars& rPars) {
	njh::stopWatch watch;
	watch.setLapName("initial set up for variant calling");
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	TranslatorByAlignmentResult ret;

	//get some length stats on input seqs
	uint64_t seqMaxLen = readVec::getMaxLength(seqs);
	seqMaxLen = seqMaxLen + rPars.realnPars.extendAmount * 2;
	auto averageLen = static_cast<uint32_t>(readVec::getAvgLength(seqs));

	//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
	{
		watch.startNewLap("set up for aligning");
		auto gprefix = bfs::path(pars_.lzPars_.genomeFnp).replace_extension("");
		auto twoBitFnp = gprefix.string() + ".2bit";
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		TwoBit::TwoBitFile tReader(twoBitFnp);
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		//get overlaping geneinfo
		auto ids = getFeatureIdsFromOverlappingRegions({refSeqRegion}, pars_.gffFnp_);
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		ret.geneIds_ = ids;
		// get gene information
		auto geneInfoDir = njh::files::make_path(pars_.workingDirtory_, "geneInfos");
		if(pars_.keepTemporaryFiles_ || pars_.writeOutGeneInfos_){
			njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});
		}
		GetGenomicLocationsForAminoAcidPositionsRet locs;
		if(!pars_.knownAminoAcidMutationsFnp_.empty()) {
			//add known to bed files
			GetGenomicLocationsForAminoAcidPositionsPars parsForBedFileGen;
			parsForBedFileGen.gffFnp = pars_.gffFnp_;
			parsForBedFileGen.twoBitFnp = twoBitFnp;
			parsForBedFileGen.doNotWrite = true;
			parsForBedFileGen.proteinMutantTypingFnp = pars_.knownAminoAcidMutationsFnp_;
			locs = getGenomicLocationsForAminoAcidPositions(parsForBedFileGen);
		}
		OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		std::unordered_map<std::string, VecStr> idToTranscriptName;
		std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> rawGenes = GeneFromGffs::getGenesFromGffForIds(pars_.gffFnp_, ids);
		// std::cout << "rawGenes.size(): " << rawGenes.size() << std::endl;

		std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;

		for(const auto & gene : rawGenes){
			bool failFilter = false;
			for(const auto & transcript : gene.second->mRNAs_){
				if(njh::in(njh::strToLowerRet(transcript->type_), VecStr{"rrna", "trna", "snorna","snrna","ncrna"}) ){
					////std::cout << __FILE__ << " " << __LINE__ << std::endl;
					transcript->writeGffRecord(std::cout);
					failFilter = true;
					break;
				}
			}
			//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
			if(!failFilter){
				 //std::cout << __FILE__ << " " << __LINE__ << std::endl;
				genes[gene.first] = gene.second;
			}
		}
		// std::cout << "genes.size(): " << genes.size() << std::endl;
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		//std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>> geneTranscriptInfos;
		uint64_t proteinMaxLen = seqMaxLen;
		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(const auto & gene : genes){
			genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
			for(const auto & transcript : gene.second->mRNAs_){
				idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
			}
			ret.transcriptInfosForGene_[gene.first] = gene.second->generateGeneSeqInfo(tReader, false);
			for(const auto & transcriptInfo : ret.transcriptInfosForGene_[gene.first]){
				ret.translationInfoForTranscirpt_[transcriptInfo.first] = transcriptInfo.second;
				if(pars_.useFullProtein_) {
					readVec::getMaxLength(transcriptInfo.second->protein_, proteinMaxLen);
				}
	//			ret.proteinForTranscript_[transcriptInfo.first] = transcriptInfo.second->protein_.seq_;
				ret.proteinVariants_.emplace(transcriptInfo.first,
										VariantsInfo{Bed3RecordCore(transcriptInfo.first, 0, transcriptInfo.second->protein_.seq_.size()), transcriptInfo.second->protein_});
			}
			if(pars_.keepTemporaryFiles_ || pars_.writeOutGeneInfos_){
				gene.second->writeOutGeneInfo(tReader, outOpts);
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		aligner alignObj(proteinMaxLen, gapScoringParameters(6,1,0,0,0,0), substituteMatrix::createDegenScoreMatrixLessN(2,-2));
		aligner protein_alignObj(proteinMaxLen, gapScoringParameters(6,1,0,0,0,0), substituteMatrix(2,-2));

		//aligner alignObjSeq(seqMaxLen + rPars.realnPars.extendAmount * 2, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2));
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		auto chromLengths = tReader.getSeqLens();
		//std::cout << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << std::endl;
		struct MinMaxPos{
			MinMaxPos()= default;
			size_t minPos_{std::numeric_limits<uint32_t>::max()};
			size_t maxPos_{0};
		};
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

		std::unordered_map<std::string, MinMaxPos> minMaxPositionsPerChrom;
		aligner alignObjAdjusted(proteinMaxLen, gapScoringParameters(6,1,0,0,0,0), substituteMatrix::createDegenScoreMatrixLessN(10, -2));
		watch.startNewLap("aligning");
		for (const auto & seq : seqs) {
			auto reAlignParsCopy = rPars.realnPars;
			if(uAbsdiff(averageLen, getSeqBase(seq).seq_.size()) > reAlignParsCopy.extendAmount){
				reAlignParsCopy.extendAmount = reAlignParsCopy.extendAmount + uAbsdiff(averageLen, getSeqBase(seq).seq_.size());
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::shared_ptr<ReAlignedSeq> initialResults;
			if(uAbsdiff(averageLen, getSeqBase(seq).seq_.size()) > 100){
				initialResults = std::make_shared<ReAlignedSeq>(ReAlignedSeq::genRealignment(seq, refSeqRegion, alignObjAdjusted, chromLengths, tReader, rPars.realnPars));
			} else {
				initialResults = std::make_shared<ReAlignedSeq>(ReAlignedSeq::genRealignment(seq, refSeqRegion, alignObj, chromLengths, tReader, rPars.realnPars));
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if(initialResults->alnRefSeq_.seq_.front() != '-' && initialResults->alnRefSeq_.seq_.back() != '-') {
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				minMaxPositionsPerChrom[initialResults->gRegion_.chrom_].minPos_ = std::min(minMaxPositionsPerChrom[initialResults->gRegion_.chrom_].minPos_,initialResults->gRegion_.start_);
				minMaxPositionsPerChrom[initialResults->gRegion_.chrom_].maxPos_ = std::max(minMaxPositionsPerChrom[initialResults->gRegion_.chrom_].maxPos_,initialResults->gRegion_.end_);

				ret.seqAlns_[getSeqBase(seq).name_].emplace_back(*initialResults);

				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				for (const auto & g : ret.geneIds_) {
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					const auto & currentGene = njh::mapAt(genes, g);
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					const auto & currentGeneInfo = njh::mapAt(ret.transcriptInfosForGene_, g);
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;

					try {
						std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes> translations;
//            std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
						translations = translateBasedOnAlignment(*initialResults, *currentGene, currentGeneInfo, protein_alignObj, pars_);
						//////std::cout << __FILE__ << " " << __LINE__ << std::endl;
//            std::cout << "translations.size(): " << translations.size() << std::endl;
						for(const auto & trans : translations){
							auto queryTransStart = trans.second.queryAlnTranslation_.seq_.find_first_not_of('-');
							if('-' == trans.second.refAlnTranslation_.seq_[queryTransStart] || countOccurences(trans.second.queryAlnTranslation_.seq_, "*") > pars_.allowableStopCodons_){
								//probably should do a more intensive check here fo
								ret.filteredOffTranslations_[getSeqBase(seq).name_].emplace(trans);
//                std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
							} else {
								ret.translations_[getSeqBase(seq).name_].emplace(trans);
//                std::cout << __PRETTY_FUNCTION__  << " " << __LINE__ << std::endl;
							}
						}
					} catch (std::exception & e) {
//						std::cerr << e.what() << std::endl;
						//Generally if an exception happen while attempting to translate alignment is mangled
						if(!njh::in(getSeqBase(seq).name_, ret.translations_)){
							ret.seqsTranslationFiltered_.emplace_back(getSeqBase(seq).name_);
						}
					}
				}
			} else {
				ret.seqsUnableToBeMapped_.emplace_back(getSeqBase(seq).name_);
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		 //////std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//list the seqs no tran
		for(const auto & filteredOff : ret.filteredOffTranslations_){
			if(!njh::in(filteredOff.first, ret.translations_)){
				ret.seqsTranslationFiltered_.emplace_back(filteredOff.first);
			}
		}

		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap("set up for variant calling");
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//index snps
		for(const auto & positons : minMaxPositionsPerChrom){
			Bed3RecordCore chromRegion(
					positons.first,
					positons.second.minPos_,
					positons.second.maxPos_);
			auto refSeq = GenomicRegion(chromRegion).extractSeq(tReader);
			ret.seqVariants_.emplace(chromRegion.chrom_, VariantsInfo(chromRegion, refSeq));
	//		for(uint32_t seqPos = 0; seqPos < refSeq.seq_.size(); ++seqPos){
	//			ret.baseForPosition_[chromRegion.chrom_][chromRegion.chromStart_ + seqPos] = refSeq.seq_[seqPos];
	//		}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap("seq variant calling");
		for(const auto & seqName : ret.seqAlns_){
			for(const auto & aln : seqName.second){
				uint32_t popCount = njh::mapAt(sampReadCountsForHaps, aln.querySeq_.name_).size();
				ret.seqVariants_.at(aln.gRegion_.chrom_).addVariantInfo(
						aln.alnRefSeq_.seq_,
						aln.alnQuerySeq_.seq_,
						popCount,
						njh::mapAt(sampReadCountsForHaps, aln.querySeq_.name_),
						aln.comp_,
						aln.gRegion_.start_);
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for (auto& varPerChrom: ret.seqVariants_) {
			for(const auto & knowLocs : locs.genomicLocs) {
				if(knowLocs.chrom_ == varPerChrom.first) {
					auto coveredPositions = getVectorOfMapKeys(varPerChrom.second.allBases);
					auto minLocs = vectorMinimum(coveredPositions);
					auto maxLocs = vectorMaximum(coveredPositions);
					for(const auto forcePositon : iter::range(knowLocs.chromStart_, knowLocs.chromEnd_)) {
						if(forcePositon >= minLocs && forcePositon <= maxLocs) {
							varPerChrom.second.alwaysReportLocations_.emplace(forcePositon);
						}
					}
				}
			}
		}
		//set finals for the snps
		std::unordered_map<std::string, std::vector<VariantsInfo::PosStartSize>> complexPositionsPerChrom;
		for (auto& varPerChrom: ret.seqVariants_) {
			varPerChrom.second.setFinals(rPars);
			{
				auto complex_positions = varPerChrom.second.getComplexPositions(rPars.complexVarPars);
				complexPositionsPerChrom[varPerChrom.first] = complex_positions;
				// std::cout << __FILE__ << " " << __LINE__ << std::endl;
				// std::cout << varPerChrom.first << std::endl;
				// for (const auto& pos: complex_positions) {
				// 	std::cout << "\t" << pos.start_ << "\t" << pos.size_ << "\t" << pos.count_ << std::endl;
				// }

			}

		}
		for(const auto & seqName : ret.seqAlns_){
			for(const auto & aln : seqName.second){
				uint32_t popCount = njh::mapAt(sampReadCountsForHaps, aln.querySeq_.name_).size();
				ret.seqVariants_.at(aln.gRegion_.chrom_).addComplexVariantInfo(
						aln.alnRefSeq_.seq_,
						aln.alnQuerySeq_.seq_,
						popCount,
						njh::mapAt(sampReadCountsForHaps, aln.querySeq_.name_),
						aln.comp_,
						aln.gRegion_.start_,
						complexPositionsPerChrom.at(aln.gRegion_.chrom_));
			}
		}

		for (auto& varPerChrom: ret.seqVariants_) {
			varPerChrom.second.setComplexFinals(rPars, complexPositionsPerChrom.at(varPerChrom.first));
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//index amino acid changes per transcript
		watch.startNewLap("protein variant calling");
		for(const auto & seqName : ret.translations_){
			for(const auto & transcript : seqName.second){
				if(countOccurences(transcript.second.queryAlnTranslation_.seq_, "*") > 1){
					//should log which ones have messed up translations
				}else{
					auto popCount = njh::mapAt(sampReadCountsForHaps, seqName.first).size();
					ret.proteinVariants_.at(transcript.first).addVariantInfo(
							transcript.second.refAlnTranslation_.seq_,
							transcript.second.queryAlnTranslation_.seq_,
							popCount,
							njh::mapAt(sampReadCountsForHaps, seqName.first),
							transcript.second.comp_,
							0);
				}
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for (auto& varPerTrans : ret.proteinVariants_) {
			for(const auto & knowLocs : locs.transcriptLocs) {
				if(knowLocs.chrom_ == varPerTrans.first) {
					auto coveredPositions = getVectorOfMapKeys(varPerTrans.second.allBases);
					auto minLocs = vectorMinimum(coveredPositions);
					auto maxLocs = vectorMaximum(coveredPositions);
					if(knowLocs.chromStart_ >= minLocs && knowLocs.chromStart_ <= maxLocs) {
						varPerTrans.second.alwaysReportLocations_.emplace(knowLocs.chromStart_);
						MetaDataInName meta(knowLocs.extraFields_[0]);
						if(meta.containsMeta("KnownAlts") && !meta.getMeta("KnownAlts").empty()) {
							auto alts = tokenizeString(meta.getMeta("KnownAlts"), ",");
							varPerTrans.second.forcedAltCalls_[knowLocs.chromStart_] = alts;
						}
					}
				}
			}
		}
		for(auto & varPerTrans : ret.proteinVariants_){
			varPerTrans.second.setFinals(rPars);
		}
		std::unordered_map<std::string, std::unordered_map<uint32_t, std::tuple<GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo,GeneSeqInfo::GenePosInfo>>> allCodonInfoByAAPos;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(const auto & seqName : ret.translations_) {
			for(const auto & transcript : seqName.second) {
				if(!njh::in(transcript.first, allCodonInfoByAAPos)) {
					allCodonInfoByAAPos[transcript.first] = ret.translationInfoForTranscirpt_.at(transcript.first)->getInfosByAAPos();
				}
			}
		}
		//type sequences for significant protein variation including codon info
		watch.startNewLap("type sequences for significant protein variation including codon info");
		for(auto & varPerTrans : ret.proteinVariants_){

			std::set<uint32_t> knownMutationsLocationsZeroBased;
			for(const auto pos : knownAminoAcidPositions_[varPerTrans.first]){
				knownMutationsLocationsZeroBased.emplace(pos - 1);
			}
			//this includes both known mutations as well as locations with significant variation
			std::set<uint32_t> allLocations = getAllInterestingAAPosZeroBased(varPerTrans.first, ret);
			for (auto & seqName : ret.translations_) {
				if (njh::in(varPerTrans.first, seqName.second)) {
					for (const auto & loc : allLocations) {
						if(loc < std::get<0>(seqName.second[varPerTrans.first].firstAminoInfo_).aaPos_ || loc > std::get<0>(seqName.second[varPerTrans.first].lastAminoInfo_).aaPos_){
							//location is not within the aligned translation
							continue;
						}
						auto codon = seqName.second[varPerTrans.first].getCodonForAARefPos(loc);
						const auto & refCodonInfo = allCodonInfoByAAPos[varPerTrans.first].at(loc);
						Codon refCodon(std::get<0>(refCodonInfo).aa_,
						std::make_tuple(std::get<0>(refCodonInfo).base_,std::get<1>(refCodonInfo).base_, std::get<2>(refCodonInfo).base_));

						ret.fullAATypedWithCodonInfo_[seqName.first].emplace_back(
								TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
										njh::in(loc, knownMutationsLocationsZeroBased), refCodon));
						ret.translated_fullAATypedWithCodonInfo_[njh::pasteAsStr(seqName.first, "[transcript=", varPerTrans.first, "]")].emplace_back(
								TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
										njh::in(loc, knownMutationsLocationsZeroBased), refCodon));
						if(njh::in(loc, varPerTrans.second.snpsFinal)){
							ret.variantAATypedWithCodonInfo_[seqName.first].emplace_back(
									TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
											njh::in(loc, knownMutationsLocationsZeroBased), refCodon));
							ret.translated_variantAATypedWithCodonInfo_[njh::pasteAsStr(seqName.first, "[transcript=", varPerTrans.first, "]")].emplace_back(
									TranslatorByAlignment::AAInfo(varPerTrans.first, loc, codon,
											njh::in(loc, knownMutationsLocationsZeroBased), refCodon));
						}
					}
				}
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> translated_fullAATypedWithCodonInfo_;
		std::map<std::string, std::vector<TranslatorByAlignment::AAInfo>> translated_variantAATypedWithCodonInfo_;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

		//sort
		for(auto & seqName : ret.fullAATypedWithCodonInfo_){
			njh::sort(seqName.second, [](const TranslatorByAlignment::AAInfo & info1, const TranslatorByAlignment::AAInfo & info2){
				if(info1.transcriptName_ == info2.transcriptName_){
					return info1.zeroBasedPos_ < info2.zeroBasedPos_;
				}else{
					return info1.transcriptName_ < info2.transcriptName_;
				}
			});
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	OutputStream timeLogOut(njh::files::make_path(pars_.workingDirtory_, "timeLog.txt"));
	watch.logLapTimes(timeLogOut, true, 6, true);
	return ret;
}


}  // namespace njhseq



