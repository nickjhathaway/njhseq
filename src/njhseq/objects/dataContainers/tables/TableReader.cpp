/*
 * TableReader.cpp
 *
 *  Created on: Jan 21, 2018
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
#include "TableReader.hpp"

namespace njhseq {


TableReader::TableReader(const TableIOOpts & tabOpts): tabOpts_(tabOpts){
	//inital header reader
	if("STDIN" == tabOpts_.in_.inFilename_){
		in_ = std::make_unique<InputStream>(tabOpts_.in_);
		if(tabOpts_.hasHeader_){
			std::string currentLine = "";
			njh::files::crossPlatGetline(*in_, currentLine);
			auto toks = tokenizeString(currentLine, tabOpts_.inDelim_, true);
			header_ = table(toks);
		}
	} else {
		njh::files::checkExistenceThrow(tabOpts_.in_.inFilename_);
		std::string currentLine = njh::files::getFirstLine(tabOpts_.in_.inFilename_);
		//std::cout << currentLine << std::endl;
		auto toks = tokenizeString(currentLine, tabOpts_.inDelim_, true);
		VecStr columnNames;
		if (!tabOpts_.hasHeader_) {
			for (const auto i : iter::range(toks.size())) {
				columnNames.emplace_back("col." + leftPadNumStr(i, toks.size()));
			}
		} else {
			columnNames = toks;
		}
		header_ = table(columnNames);
		in_ = std::make_unique<InputStream>(tabOpts_.in_);
		if(tabOpts_.hasHeader_){
			njh::files::crossPlatGetline(*in_, currentLine);
		}
	}
}

void TableReader::setHeaderlessHeader(uint32_t numOfCols) {
	VecStr columnNames;
	for (const auto i : iter::range(numOfCols)) {
		columnNames.emplace_back("col." + leftPadNumStr(i, numOfCols));
	}
	header_ = table(columnNames);
}


bool TableReader::getNextRow(VecStr & row){
	std::string currentLine = "";
	row.clear();
	if(njh::files::crossPlatGetline(*in_, currentLine)){
		row = tokenizeString(currentLine, tabOpts_.inDelim_, true);
		if(row.size() != header_.nCol()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error the row has a different number of columns than the first line" << "\n";
			ss << "rowSize: " << row.size() << ", firstLineSize: " << header_.nCol() << "\n";
			ss << "row: " << currentLine << "\n";
			throw std::runtime_error{ss.str()};
		}
		return true;
	}
	return false;
}

VecStr TableReader::extractCols(const VecStr & row, const VecStr & cols) const{
	VecStr ret;
	ret.reserve(cols.size());
	for(const auto & col : cols){
		ret.emplace_back(row[header_.getColPos(col)]);
	}
	return ret;
}

}  // namespace njhseq
