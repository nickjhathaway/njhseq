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
#include "bibseq/objects/BioDataObject/swisProt.hpp"
#include "bibseq/objects/BioDataObject/TandemRepeatFinderRecord.hpp"


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





}  // namespace bibseq

