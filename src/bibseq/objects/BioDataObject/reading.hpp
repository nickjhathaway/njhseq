#pragma once
//
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
/*
 * reading.hpp
 *
 *  Created on: Mar 1, 2016
 *      Author: nick
 */

#include "bibseq/objects/BioDataObject/BedRecordCore.hpp"
#include "bibseq/objects/BioDataObject/GFFCore.hpp"
#include "bibseq/objects/BioDataObject/RefSeqGeneRecord.hpp"
#include "bibseq/objects/BioDataObject/RepeatMaskerRecord.hpp"
#include "bibseq/objects/BioDataObject/GenomicRegion.hpp"

#include "bibseq/objects/BioDataObject/swisProt.hpp"
#include "bibseq/objects/BioDataObject/TandemRepeatFinderRecord.hpp"
#include "bibseq/objects/BioDataObject/BioDataFileIO.hpp"


namespace bibseq {
std::vector<std::shared_ptr<GFFCore>> getGFFs(const bfs::path & filename);

std::vector<std::shared_ptr<Bed6RecordCore>> getBeds(const bfs::path & filename);

std::vector<std::shared_ptr<Bed3RecordCore>> getBed3s(const bfs::path & filename);

std::vector<std::shared_ptr<RefSeqGeneRecord>> getRefSeqGeneRecords(const bfs::path & filename);

std::vector<std::shared_ptr<RepeatMaskerRecord>> getRepeatMaskerRecords(const bfs::path & filename);


std::vector<std::shared_ptr<TandemRepeatFinderRecord>> getTandemRepeatFinderRecords(const bfs::path & filename);


void checkPositionSortedBedThrow(const bfs::path & bedFnp,
		const std::string & funcName);

void checkPositionSortedNoOverlapsBedThrow(const bfs::path & bedFnp,
		const std::string & funcName);

struct intersectBedLocsWtihGffRecordsPars {
	intersectBedLocsWtihGffRecordsPars();
	intersectBedLocsWtihGffRecordsPars(const bfs::path & gffFnp);
	intersectBedLocsWtihGffRecordsPars(const bfs::path & gffFnp,
			const VecStr & extraAttributes, const VecStr & selectFeatures);
	bfs::path gffFnp_;
	VecStr extraAttributes_;
	VecStr selectFeatures_;

	Json::Value toJson() const;
};

//Json::Value intersectBedLocsWtihGffRecords(
//		const std::vector<std::shared_ptr<Bed3RecordCore>> & beds,
//		const intersectBedLocsWtihGffRecordsPars & pars);

template <typename BEDREC>
Json::Value intersectBedLocsWtihGffRecords(
		std::vector<BEDREC> & beds,
		const intersectBedLocsWtihGffRecordsPars & pars){
	Json::Value ret;
	for (auto & inputRegion : beds) {
		getRef(inputRegion).extraFields_.emplace_back("");
	}
	std::unordered_map<std::string, std::vector<uint32_t>> bedsByChrome;

	BioDataFileIO<GFFCore> reader { IoOptions(InOptions(pars.gffFnp_)) };
	reader.openIn();
	uint32_t count = 0;
	std::string line = "";
	std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
	for (const auto & bPos : iter::range(beds.size())) {
		bedsByChrome[getRef(beds[bPos]).chrom_].emplace_back(bPos);
	}
	while (nullptr != gRecord) {
		if (pars.selectFeatures_.empty() || bib::in(gRecord->type_, pars.selectFeatures_)) {
			auto gRegion = GenomicRegion(*gRecord);
			for (auto & inputRegionPos : bedsByChrome[gRegion.chrom_]) {

				if (GenomicRegion(getRef(beds[inputRegionPos])).overlaps(gRegion)) {
					if("" != getRef(beds[inputRegionPos]).extraFields_.back()){
						getRef(beds[inputRegionPos]).extraFields_.back().append(",");
					}
					getRef(beds[inputRegionPos]).extraFields_.back().append("[");
					getRef(beds[inputRegionPos]).extraFields_.back().append(
							"ID=" + gRecord->getAttr("ID") + ";");
					if(pars.selectFeatures_.empty() || 1 != pars.selectFeatures_.size()){
						getRef(beds[inputRegionPos]).extraFields_.back().append("feature=" + gRecord->type_ + ";");
					}
					if (!ret.isMember(gRecord->getAttr("ID"))) {
						ret[gRecord->getAttr("ID")] = gRecord->toJson();
					}
					for (const auto & attr : pars.extraAttributes_) {
						if (gRecord->hasAttr(attr)) {
							getRef(beds[inputRegionPos]).extraFields_.back().append(
									attr + "=" + gRecord->getAttr(attr) + ";");
						} else {
							getRef(beds[inputRegionPos]).extraFields_.back().append(
									attr + "=" + "NA" + ";");
						}
					}
					getRef(beds[inputRegionPos]).extraFields_.back().append("]");
				}
			}
		}
		bool end = false;
		while ('#' == reader.inFile_->peek()) {
			if (bib::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
				end = true;
				break;
			}
			bib::files::crossPlatGetline(*reader.inFile_, line);
		}
		if (end) {
			break;
		}
		gRecord = reader.readNextRecord();
		++count;
	}
	return ret;
}



}  // namespace bibseq

