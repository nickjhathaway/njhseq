/*
 * BamCountExtractStats.cpp
 *
 *  Created on: Feb 29, 2016
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

#include "BamCountExtractStats.hpp"

namespace njhseq {

void BamCountExtractStats::addStats(const BamCountExtractStats & otherStats) {
	highQualityBases += otherStats.highQualityBases;
	lowQaulityBases += otherStats.lowQaulityBases;
	readsBellowMappingQuality += otherStats.readsBellowMappingQuality;
	readsBellowMinLen += otherStats.readsBellowMinLen;
	readsMapped += otherStats.readsMapped;
	readsNotMapped += otherStats.readsNotMapped;
	readsUsed += otherStats.readsUsed;
	totalReads += otherStats.totalReads;
	disCordantMapping += otherStats.disCordantMapping;
	largeInsertSize += otherStats.largeInsertSize;
	pairedReadsTotal += otherStats.pairedReadsTotal;
	pairedReadsUsed += otherStats.pairedReadsUsed;
	overLappingPairedReads += otherStats.overLappingPairedReads;
	failedToFindMate += otherStats.failedToFindMate;
}

void BamCountExtractStats::createBaseFilterFile(const OutOptions & fileOpts,
		uint32_t qualCutOff) const {
	std::ofstream baseFilterInfoFile;
	openTextFile(baseFilterInfoFile, fileOpts);
	baseFilterInfoFile << "totalBases\tbases>" << qualCutOff << "\tbase<="
			<< qualCutOff << std::endl;
	baseFilterInfoFile << highQualityBases + lowQaulityBases << "\t"
			<< getPercentageString(highQualityBases,
					highQualityBases + lowQaulityBases) << "\t"
			<< getPercentageString(lowQaulityBases,
					highQualityBases + lowQaulityBases) << std::endl;
}

void BamCountExtractStats::createReadFilterFile(const OutOptions & fileOpts, uint32_t minLen,
		uint32_t mappingQual, uint32_t insertSizeCutOff) const {
	std::ofstream filterInfoFile;
	openTextFile(filterInfoFile, fileOpts);
	filterInfoFile << "totalReads\treadsUsed\tpairedReads\tpairedReadsUsed\tlen<" << minLen
			<< "\tMappingQuality<" << mappingQual << "\treadsNotMapped"
			<< "\tlargeInsertSize>" << insertSizeCutOff
			<< "\tfailedToFindMate\tdiscordantMapping"
			<< std::endl;
	filterInfoFile << totalReads << "\t"
			<< getPercentageString(readsUsed, totalReads) << "\t"
			<< getPercentageString(pairedReadsTotal, totalReads) << "\t"
			<< getPercentageString(pairedReadsUsed, totalReads) << "\t"
			<< getPercentageString(readsBellowMinLen, totalReads) << "\t"
			<< getPercentageString(readsBellowMappingQuality, totalReads) << "\t"
			<< getPercentageString(readsNotMapped, totalReads) << "\t"
			<< getPercentageString(largeInsertSize, totalReads) << "\t"
			<< getPercentageString(failedToFindMate, totalReads) << "\t"
			<< getPercentageString(disCordantMapping, totalReads) << std::endl;
}

void BamCountExtractStats::createPairedReadInfoFile(const OutOptions & fileOpts) const{
	std::ofstream filterInfoFile;
	openTextFile(filterInfoFile, fileOpts);
	filterInfoFile << "totalReads\treadsPairedTotal\treadsPairedUsed\treadsOverlapping\tfailedToFindMate"<< std::endl;
	filterInfoFile << totalReads << "\t"
			<< getPercentageString(pairedReadsTotal, totalReads) << "\t"
			<< getPercentageString(pairedReadsUsed, pairedReadsTotal) << "\t"
			<< getPercentageString(overLappingPairedReads, pairedReadsUsed) << "\t"
			<< getPercentageString(failedToFindMate, pairedReadsTotal)
			<< std::endl;
}


void BamCountExtractStats::createPairedReadInfoFileMultiple(
			const std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
			const OutOptions & fileOpts){
	BamCountExtractStats allStats;
	for (const auto & stat : refExtracts) {
		if (0 == stat.second.readsUsed) {
			continue;
		}
		allStats.addStats(stat.second);
	}
	std::ofstream filterInfoFile;
	//OutOptions(setUp.pars_.directoryName_ + "baseFilterInfo.tab.txt")
	openTextFile(filterInfoFile, fileOpts);
	filterInfoFile << "refName\ttotalReads\treadsPairedTotal\treadsPairedUsed\treadsOverlapping\tfailedToFindMate"<< std::endl;

	auto statKeys = getVectorOfMapKeys(refExtracts);
	njh::sort(statKeys);
	for (const auto & statKey : statKeys) {
		const auto & stat = refExtracts.at(statKey);
		if (0 == stat.readsUsed) {
			continue;
		}
		filterInfoFile << statKey << "\t"
		    << stat.totalReads << "\t"
				<< getPercentageString(stat.pairedReadsTotal, stat.totalReads) << "\t"
				<< getPercentageString(stat.pairedReadsUsed,stat.pairedReadsTotal) << "\t"
				<< getPercentageString(stat.overLappingPairedReads, stat.pairedReadsUsed) << "\t"
				<< getPercentageString(stat.failedToFindMate, stat.pairedReadsTotal)
				<< std::endl;
	}
	filterInfoFile << "all" << "\t"
	    << allStats.totalReads << "\t"
			<< getPercentageString(allStats.pairedReadsTotal, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.pairedReadsUsed,allStats.pairedReadsTotal) << "\t"
			<< getPercentageString(allStats.overLappingPairedReads, allStats.pairedReadsUsed) << '\t'
			<< getPercentageString(allStats.failedToFindMate, allStats.pairedReadsTotal)
			<< std::endl;

}
void BamCountExtractStats::createBaseFilterFileMultiple(
		const std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
		const OutOptions & fileOpts, uint32_t qualCutOff) {
	BamCountExtractStats allStats;
	for (const auto & stat : refExtracts) {
		if (0 == stat.second.readsUsed) {
			continue;
		}
		allStats.addStats(stat.second);
	}
	std::ofstream baseFilterInfoFile;
	//OutOptions(setUp.pars_.directoryName_ + "baseFilterInfo.tab.txt")
	openTextFile(baseFilterInfoFile, fileOpts);
	baseFilterInfoFile << "refName\ttotalBases\tbases>" << qualCutOff
			<< "\tbase<=" << qualCutOff << std::endl;
	auto statKeys = getVectorOfMapKeys(refExtracts);
	njh::sort(statKeys);
	for (const auto & statKey : statKeys) {
		const auto & stat = refExtracts.at(statKey);
		if (0 == stat.readsUsed) {
			continue;
		}
		baseFilterInfoFile << statKey << "\t"
				<< stat.highQualityBases + stat.lowQaulityBases << "\t"
				<< getPercentageString(stat.highQualityBases,
						stat.highQualityBases + stat.lowQaulityBases)
				<< "\t"
				<< getPercentageString(stat.lowQaulityBases,
						stat.highQualityBases + stat.lowQaulityBases)
				<< std::endl;
	}
	baseFilterInfoFile << "all" << "\t"
			<< allStats.highQualityBases + allStats.lowQaulityBases << "\t"
			<< getPercentageString(allStats.highQualityBases,
					allStats.highQualityBases + allStats.lowQaulityBases) << "\t"
			<< getPercentageString(allStats.lowQaulityBases,
					allStats.highQualityBases + allStats.lowQaulityBases) << std::endl;
}

void BamCountExtractStats::createReadFilterFileMultiple(
		const std::unordered_map<std::string, BamCountExtractStats>& refExtracts,
		const OutOptions & fileOpts,
		uint32_t minLen, uint32_t mappingQual,
		uint32_t insertSizeCutOff) {
	BamCountExtractStats allStats;

	for (const auto & stat : refExtracts) {
		if (0 == stat.second.readsUsed) {
			continue;
		}
		allStats.addStats(stat.second);
	}
	std::ofstream filterInfoFile;
	openTextFile(filterInfoFile, fileOpts);
	filterInfoFile << "refName\ttotalReads\treadsUsed\tpairedReads\tpairedReadsUsed\tlen<" << minLen
			<< "\tMappingQuality<" << mappingQual << "\treadsNotMapped"
			<< "\tlargeInsertSize>" << insertSizeCutOff
			<< "\tfailedToFindMate\tdiscordantMapping"
			<< std::endl;
	auto statKeys = getVectorOfMapKeys(refExtracts);
	njh::sort(statKeys);
	for (const auto & statKey : statKeys) {
		const auto & stat = refExtracts.at(statKey);
		if (0 == stat.readsUsed) {
			continue;
		}
		filterInfoFile << statKey << "\t" << stat.totalReads << "\t"
				<< getPercentageString(stat.readsUsed, stat.totalReads) << "\t"
				<< getPercentageString(stat.pairedReadsTotal, stat.totalReads) << "\t"
				<< getPercentageString(stat.pairedReadsUsed, stat.totalReads) << "\t"
				<< getPercentageString(stat.readsBellowMinLen, stat.totalReads) << "\t"
				<< getPercentageString(stat.readsBellowMappingQuality, stat.totalReads) << "\t"
				<<	0 << "\t"
				<< getPercentageString(stat.largeInsertSize, stat.totalReads) << "\t"
				<< getPercentageString(stat.failedToFindMate, stat.totalReads) << "\t"
				<< getPercentageString(stat.disCordantMapping, stat.totalReads) << std::endl;

	}
	filterInfoFile << "all" << "\t" << allStats.totalReads << "\t"
			<< getPercentageString(allStats.readsUsed, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.pairedReadsTotal, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.pairedReadsUsed, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.readsBellowMinLen, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.readsBellowMappingQuality, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.readsNotMapped, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.largeInsertSize, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.failedToFindMate, allStats.totalReads) << "\t"
			<< getPercentageString(allStats.disCordantMapping, allStats.totalReads) << std::endl;
}

}  // namespace njhseq

