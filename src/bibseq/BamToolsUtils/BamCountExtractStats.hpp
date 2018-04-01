#pragma once
/*
 * BamCountExtractStats.hpp
 *
 *  Created on: Feb 29, 2016
 *      Author: nick
 */

#include "bibseq/common.h"
#include "bibseq/IO.h"

namespace bibseq {

class BamCountExtractStats {
public:
	uint32_t readsMapped = 0;
	uint32_t readsUsed = 0;
	uint32_t readsNotMapped = 0;
	uint32_t readsBellowMappingQuality = 0;
	uint32_t readsBellowMinLen = 0;
	uint32_t totalReads = 0;
	uint64_t highQualityBases = 0;
	uint64_t lowQaulityBases = 0;
	uint64_t largeInsertSize = 0;
	uint64_t disCordantMapping = 0;
	uint64_t pairedReadsTotal = 0;
	uint64_t pairedReadsUsed = 0;
	uint64_t overLappingPairedReads = 0;
	uint64_t failedToFindMate = 0;


	void addStats(const BamCountExtractStats & otherStats);

	void createBaseFilterFile(const OutOptions & fileOpts,
			uint32_t qualCutOff) const;



	void createReadFilterFile(const OutOptions & fileOpts, uint32_t minLen,
			uint32_t mappingQual, uint32_t insertSizeCutOff) const;

	void createPairedReadInfoFile(const OutOptions & fileOpts) const;

	static void createBaseFilterFileMultiple(
			const std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
			const OutOptions & fileOpts, uint32_t qualCutOff);

	static void createPairedReadInfoFileMultiple(
			const std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
			const OutOptions & fileOpts);

	static void createReadFilterFileMultiple(
			const std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
			const OutOptions & fileOpts, uint32_t minLen, uint32_t mappingQual,
			uint32_t insertSizeCutOff);
};

}  // namespace bibseq

