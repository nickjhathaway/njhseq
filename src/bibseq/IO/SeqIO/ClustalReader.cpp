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
#include "ClustalReader.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {







std::vector<readObject> ClustalReader::readClustal(std::string filename, bool processed) {
	std::vector<readObject> ret;
	table inTab(filename);
	std::vector<std::pair<std::string, readObject>> readMap;
	std::vector<std::pair<std::string, readObject>>::iterator readMapIter;
	for (const auto& row : inTab.content_) {
		if (row.size() != 2 || row[0][0] == '*') {
		} else {
			bool foundMatch = false;
			for (readMapIter = readMap.begin(); readMapIter != readMap.end();
					++readMapIter) {
				if (readMapIter->first == row[0]) {
					foundMatch = true;
					readMapIter->second.seqBase_.append(row[1]);
					break;
				}
			}
			if (!foundMatch) {
				readMap.emplace_back(row[0], readObject(seqInfo(row[0], row[1])));
			}
		}
	}
	for (readMapIter = readMap.begin(); readMapIter != readMap.end();
			readMapIter++) {
		ret.emplace_back(readMapIter->second);
	}
	return ret;
}

std::vector<readObject> ClustalReader::readClustalNg(std::string filename, bool processed) {
	std::vector<readObject> ret;
	ret = readClustal(filename, processed);
	readVec::removeGapsFromReads(ret);
	return ret;
}









}  // namespace bib
