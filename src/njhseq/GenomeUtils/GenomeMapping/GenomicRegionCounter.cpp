/*
 * GenomicRegionCounter.cpp
 *
 *  Created on: Mar 31, 2017
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

#include "GenomicRegionCounter.hpp"
#include "njhseq/objects/BioDataObject/BioDataFileIO.hpp"
#include "njhseq/objects/BioDataObject/GFFCore.hpp"
#include "njhseq/BamToolsUtils.h"



namespace njhseq {

GenomicRegionCounter::GenomicRegionCount::GenomicRegionCount(
		const GenomicRegion & region) :
		region_(region) {
}
GenomicRegionCounter::GenomicRegionCount::GenomicRegionCount(
		const GenomicRegion & region, uint32_t count) :
		region_(region), count_(count) {
}

bool GenomicRegionCounter::GenomicRegionCount::sameRegion(
		const GenomicRegion& region) const {
	return region_.sameRegion(region);
}
void GenomicRegionCounter::GenomicRegionCount::increaseCount(uint32_t count) {
	count_ += count;
}

bool GenomicRegionCounter::hasRegion(const GenomicRegion & region) const {
	auto search = counts_.find(region.createUidFromCoords());
	if (counts_.end() != search) {
		return true;
	}
	return false;
}

void GenomicRegionCounter::increaseCount(const GenomicRegion & region,
		uint32_t count) {
	if (hasRegion(region)) {
		counts_.at(region.createUidFromCoords()).increaseCount(count);
	} else {
		GenomicRegionCount counter(region, count);
		counts_.emplace(region.createUidFromCoords(), counter);
	}
}

std::vector<GenomicRegion> GenomicRegionCounter::getRegionsLargestOnTop() const {
	auto uids = getVectorOfMapKeys(counts_);
	njh::sort(uids,
			[this](const std::string & key1, const std::string & key2) {return counts_.at(key1).count_ > counts_.at(key2).count_;});
	std::vector<GenomicRegion> ret;
	for (const auto & key : uids) {
		ret.emplace_back(counts_.at(key).region_);
	}
	return ret;
}

void GenomicRegionCounter::increaseCounts(
		const std::vector<BamTools::BamAlignment> & bAlns,
		const BamTools::RefVector & refData) {
	for (const auto & bAln : bAlns) {
		GenomicRegion gRegion(bAln, refData);
		gRegion.setUidWtihCoords();
		increaseCount(gRegion, 1);
	}
}


std::set<std::string> GenomicRegionCounter::getIntersectingGffIds(const bfs::path & gffFnp, const VecStr & features)const {
	std::set<std::string> idsFromData;
	BioDataFileIO<GFFCore> reader { IoOptions(InOptions(gffFnp)) };
	reader.openIn();
	// uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	while (nullptr != gRecord) {
		if (njh::in(gRecord->type_, features) ) {
			for (const auto & gCount : counts_) {
				if (gCount.second.region_.overlaps(*gRecord)) {
					idsFromData.emplace(gRecord->getIDAttr());
					break;
				}
			}
		}
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (njh::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			njh::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gRecord = reader.readNextRecord();
		// ++count;
	}
	return idsFromData;
}



GenomicRegionCounter GenomicRegionCounter::countRegionsInBam(const bfs::path & bamFnp){
		GenomicRegionCounter gCounter;

		BamTools::BamReader bReader;
		bReader.Open(bamFnp.string());
		checkBamOpenThrow(bReader, bamFnp);

		BamTools::BamAlignment bAln;
		auto refData = bReader.GetReferenceData();
		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped()) {
				gCounter.increaseCount(GenomicRegion(bAln, refData), 1);
			}
		}
		return gCounter;
	}

}  // namespace njhseq

