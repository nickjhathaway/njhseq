#pragma once
/*
 * bamExtractUtils.hpp
 *
 *  Created on: Feb 2, 2017
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
#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"
#include "njhseq/objects/seqObjects/Paired/PairedRead.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers/readChecker.hpp"


namespace njhseq {



class BamExtractor {
public:


	struct BamExtractSeqsResults {
		std::vector<PairedRead> pairs_;
		std::vector<seqInfo> singlets_;
	};

	struct ExtractCounts {
		uint32_t pairedReads_ = 0;
		uint32_t pairedReadsMateUnmapped_ = 0;//! one mate is mapped and the other mate is unmapped
		uint32_t pairedReadsMateUnmappedFailedSoftClip_ = 0;//! one mate is mapped and the other mate is unmapped but the mapped mate failed soft-clip filter
		uint32_t pairedReadsMateFailedSoftClip_ = 0;//! one mate mapped fully but the other mate failed a soft clip filter
		uint32_t pairedReadsBothFailedSoftClip_ = 0;//! one mate mapped fully but the other mate failed a soft clip filter
		uint32_t unpaiedReads_ = 0;
		uint32_t orphans_ = 0; /**< reads that are paired but their mates weren't found */
		uint32_t orphansFiltered_ = 0; /**< reads that are paired but their mates weren't found and eventually are filtered off */

		uint32_t orphansUnmapped_ = 0; /**< reads that are paired but their mates weren't found are unmapped */
		uint32_t pairsUnMapped_ = 0;
		uint32_t unpairedUnMapped_ = 0;
		uint32_t unpairedFailedSoftClip_ = 0;//! single but failed a soft clip filter


//		uint32_t mateFilteredOff_ = 0; //! mate has been filtered
//		uint32_t mateFilteredOffUnmapped_ = 0;//! mate was unmapped

		uint32_t bothMatesFilteredOff_ = 0; //! both mates has been filtered
		uint32_t singlesFilteredOff_ = 0; //! single has been filtered


		uint32_t mateFilteredOff_ = 0; //! mate has been filtered
		uint32_t mateFilteredOffUnmapped_ = 0;//! mate was unmapped

		uint32_t discordant_ = 0;
		uint32_t inverse_ = 0;


		uint32_t getTotal();
		void log(std::ostream & out, const bfs::path & bamFnp) ;
	};

	struct BamExtractSeqsResultsAlns : public ExtractCounts{
		std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment>> pairs_;
		std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment>> pairsMateUnmapped_;
		std::vector<BamTools::BamAlignment> singlets_;

		std::vector<BamTools::BamAlignment> thrownAwayMates_;

		std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment>> pairsUnmapped_;
		std::vector<BamTools::BamAlignment> singletsUnmapped_;

		std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment>> inversePairs_;
		std::vector<std::pair<BamTools::BamAlignment, BamTools::BamAlignment>> discordantPairs_;

	};

	struct ExtractedFilesOpts :  public ExtractCounts{
		SeqIOOptions inPairs_;
		SeqIOOptions inPairsMateUnmapped_;
		SeqIOOptions inThrownAwayUnmappedMate_;
		SeqIOOptions inUnpaired_;
		SeqIOOptions inPairsUnMapped_;
		SeqIOOptions inUnpairedUnMapped_;

		SeqIOOptions inFilteredPairs_;
		SeqIOOptions inFilteredSingles_;

		SeqIOOptions inInverse_;
		SeqIOOptions inDiscordant_;
		void removeAllInFiles();

	};

	struct ExtractedFilesWithStichingOpts : public ExtractedFilesOpts{
		SeqIOOptions notCombinedPairs_;
		SeqIOOptions stitched_;
	};

	BamExtractor(bool verbose = false);


	BamExtractSeqsResults extractReadsFromBamRegion(
			const bfs::path & bamFnp,
			const GenomicRegion & region,
			double percInRegion);

	void writeExtractReadsFromBamRegion(
			const bfs::path & bamFnp,
			const GenomicRegion & region,
			double percInRegion,
			const OutOptions & outOpts);

//	ExtractedFilesOpts writeExtractReadsFromBamRegionHandelOrientation(
//			const SeqIOOptions & inOutOpts,
//			GenomicRegion & region,
//			double percInRegion,
//			bool originalOrientation,
//			bool throwAwayUnmappedMate);


	BamExtractSeqsResultsAlns extractReadsFromBamRegionAlns(
			const bfs::path & bamFnp,
			const GenomicRegion & region,
			double percInRegion);

	ExtractedFilesWithStichingOpts writeExtractReadsFromBamRegionStitch(
			const bfs::path & bamFnp, const GenomicRegion & region,
			double percInRegion, const OutOptions & outOpts,
			const std::string & extraFlashCms = "");

	void writeExtractReadsFromBam(const bfs::path & bamFnp,
			const OutOptions & outOpts);

	void writeExtractReadsFromBamOnlyMapped(const bfs::path & bamFnp,
			const OutOptions & outOpts);



	/**@brief extract concordant mapping pairs and singles, no unmapped or discordant alignments extracted
	 *
	 * @param bamFnp the file path to the bam file to extract from
	 * @return a struct with the alignments extracted
	 */
	//BamExtractSeqsResultsAlns extractReadsFromBamAlns(const bfs::path & bamFnp);

	ExtractedFilesOpts extractReadsFromBamToSameOrientationContigs(
			const SeqIOOptions & opts,
			bool throwAwayUnmmpaedMates);

	ExtractedFilesOpts extractReadsFromBamWrite(const SeqIOOptions & opts, bool referenceOrientation, bool throwAwayUnmmpaedMates);
	ExtractedFilesOpts extractReadsFromBamWriteAsSingles(const SeqIOOptions & opts, bool referenceOrientation = false);


	struct extractReadsWtihCrossRegionMappingPars {
		double percInRegion_ = 0.5;
		bool originalOrientation_ = false;
		bool throwAwayUnmappedMate_ = false;
		bool tryToFindOrphansMate_ = false;
		bool keepMarkedDuplicate_ = false;
		uint32_t minAlnMapSize_ = 35;
		Json::Value toJson() const;

	};
	ExtractedFilesOpts extractReadsWtihCrossRegionMapping(
			const SeqIOOptions & inOutOpts,
			const std::vector<GenomicRegion> & regions,
			const extractReadsWtihCrossRegionMappingPars & extractPars);

	ExtractedFilesOpts extractReadsWtihCrossRegionMappingAsSingles(
				const SeqIOOptions & inOutOpts,
				const std::vector<GenomicRegion> & regions, double percInRegion,
				bool originalOrientation);


	ExtractedFilesOpts writeUnMappedSeqs(const SeqIOOptions & opts);
	struct writeUnMappedSeqsAndSmallAlnsWithFiltersPars {
		uint32_t maxAlnSize_ { 35 };
		uint32_t minSotClip_ { 40 };
		uint32_t minQuerySize_ { 36 };
		bool tryToFindOrphanMate_ { false };

		bool filterOffLowEntropyShortAlnsRecruits_{true};
		double filterOffLowEntropyShortAlnsRecruitsCutOff_{1.5};
		uint32_t entropyKlen_{2};

		Json::Value toJson() const;
	};

	ExtractedFilesOpts writeUnMappedSeqsAndSmallAlnsWithFilters(const SeqIOOptions & opts,
			const ReadCheckerQualCheck & qualChecker,
			const writeUnMappedSeqsAndSmallAlnsWithFiltersPars & alnSizeFilt,
			const std::vector<GenomicRegion> & skipRegions);


	bool verbose_ = false;
	bool debug_ = false;

	uint32_t insertLengthCutOff_ = 1000;

};




}  // namespace njhseq



