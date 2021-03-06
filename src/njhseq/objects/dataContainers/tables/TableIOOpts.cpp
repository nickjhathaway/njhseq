/*
 * TableIOOpts.cpp
 *
 *  Created on: May 30, 2016
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

#include "TableIOOpts.hpp"

namespace njhseq {

TableIOOpts TableIOOpts::genTabFileOut(const bfs::path & outFilename, bool header){
	bool zip = njh::endsWith(outFilename.string(), ".gz");
	return TableIOOpts(InOptions(), "\t", OutOptions(outFilename.string(), zip ? ".tab.txt.gz" : ".tab.txt", "tab"), "\t", header);
}

TableIOOpts TableIOOpts::genTabFileIn(const bfs::path & inFilename, bool header){
	bool zip = njh::endsWith(inFilename.string(), ".gz");
	return TableIOOpts(InOptions(inFilename), "\t", OutOptions("", zip ? ".tab.txt.gz" : ".tab.txt", "tab"), "\t", header);
}

TableIOOpts TableIOOpts::genCommaFileOut(const bfs::path & outFilename, bool header){
	bool zip = njh::endsWith(outFilename.string(), ".gz");
	return TableIOOpts(InOptions(), ",", OutOptions(outFilename.string(), zip ? ".csv.gz" : ".csv", "comma"), ",", header);
}

TableIOOpts TableIOOpts::genCommaFileIn(const bfs::path & inFilename, bool header){
	bool zip = njh::endsWith(inFilename.string(), ".gz");
	return TableIOOpts(InOptions(inFilename), ",", OutOptions("", zip ? ".csv.gz" : ".csv", "comma"), ",", header);
}

TableIOOpts::TableIOOpts() :
		IoOptions(OutOptions("", ".txt", "txt")) {

}
TableIOOpts::TableIOOpts(const InOptions & inOpts) :
		IoOptions(inOpts, OutOptions("", ".txt", "txt")) {

}

TableIOOpts::TableIOOpts(const OutOptions & outOpts) {

}

TableIOOpts::TableIOOpts(const InOptions & inOpts, const std::string & inDelim,
		bool header) :
		IoOptions(inOpts, OutOptions("", ".txt", "txt")), inDelim_(inDelim), hasHeader_(
				header) {
	if("tab" == inDelim_){
		inDelim_ = "\t";
	}
}

TableIOOpts::TableIOOpts(const OutOptions & outOpts,
		const std::string & outDelim, bool header) :
		IoOptions(outOpts), outDelim_(outDelim), hasHeader_(header) {

}

TableIOOpts::TableIOOpts(const InOptions & inOpts, const OutOptions & outOpts) :
		IoOptions(inOpts, outOpts) {

}

TableIOOpts::TableIOOpts(const InOptions & inOpts, const std::string & inDelim,
		const OutOptions & outOpts, const std::string & outDelim, bool header) :
		IoOptions(inOpts, outOpts), inDelim_(inDelim), outDelim_(outDelim), hasHeader_(
				header) {
	if("tab" == inDelim_){
		inDelim_ = "\t";
	}
}


}  // namespace njhseq


