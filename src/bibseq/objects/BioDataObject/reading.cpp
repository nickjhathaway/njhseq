//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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
 * reading.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: nick
 */

#include "reading.hpp"

namespace bibseq {

std::vector<std::shared_ptr<GFFCore>> getGFFs(const std::string & filename) {
	std::vector<std::shared_ptr<GFFCore>> ret;
	DataFileReader reader(filename);
	uint32_t count = 0;
	std::shared_ptr<GFFCore> gff = reader.readNextRecord<GFFCore>();
	while (gff) {
		ret.emplace_back(gff);
		if (bib::files::nextLineBeginsWith(reader.file_, "##FASTA")) {
			break;
		}
		gff = reader.readNextRecord<GFFCore>();
		++count;
	}
	return ret;
}

std::vector<std::shared_ptr<BedRecordCore>> getBeds(
		const std::string & filename) {
	DataFileReader bedReader(filename);
	std::vector<std::shared_ptr<BedRecordCore>> beds;
	std::shared_ptr<BedRecordCore> bed =
			bedReader.readNextRecord<BedRecordCore>();
	while (bed) {
		beds.emplace_back(bed);
		bed = bedReader.readNextRecord<BedRecordCore>();
	}
	return beds;
}

}  // namespace bibseq

