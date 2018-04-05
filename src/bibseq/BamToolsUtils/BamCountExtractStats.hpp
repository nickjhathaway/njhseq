#pragma once
/*
 * BamCountExtractStats.hpp
 *
 *  Created on: Feb 29, 2016
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
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

