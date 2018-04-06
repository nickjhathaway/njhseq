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

#include "ClustalReader.hpp"

#include <bibcpp/bashUtils.h>


#include "bibseq/readVectorManipulation.h"


namespace bibseq {

std::vector<readObject> ClustalReader::readClustal(std::string filename, bool processed) {
	std::vector<readObject> ret;
	//table inTab(filename, "whitespace", false);
	if(!bfs::exists(filename)){
		std::stringstream ss;
		ss <<__PRETTY_FUNCTION__ << ", error file " << filename << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::ifstream infile(filename);
	if(!infile){
		std::stringstream ss;
		ss <<__PRETTY_FUNCTION__ << ", error in opening " << filename << "" << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::vector<std::pair<std::string, readObject>> readMap;
	std::string line = "";

	while (bib::files::crossPlatGetline(infile, line)) {
		auto row = tokenizeString(line, "whitespace");
		if(row.size() == 2 && row[0][0] != '*' && !allWhiteSpaceStr(row.front())){
			bool foundMatch = false;
			for (auto & readMapIter : readMap) {
				if (readMapIter.first == row[0]) {
					foundMatch = true;
					readMapIter.second.seqBase_.append(row[1]);
					break;
				}
			}
			if (!foundMatch) {
				readMap.emplace_back(row[0], readObject(seqInfo(row[0], row[1])));
			}
		}
	}

	for (const auto & readMapIter : readMap) {
		ret.emplace_back(readMapIter.second);
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
