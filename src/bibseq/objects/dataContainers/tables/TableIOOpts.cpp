/*
 * TableIOOpts.cpp
 *
 *  Created on: May 30, 2016
 *      Author: nick
 */
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

#include "TableIOOpts.hpp"

namespace bibseq {

TableIOOpts TableIOOpts::genTabFileOut(const bfs::path & outFilename, bool header){
	return TableIOOpts(InOptions(), "\t", OutOptions(outFilename.string(), ".tab.txt", "tab"), "\t", header);
}

TableIOOpts TableIOOpts::genTabFileIn(const bfs::path & inFilename, bool header){
	return TableIOOpts(InOptions(inFilename), "\t", OutOptions("", ".tab.txt", "tab"), "\t", header);
}

TableIOOpts TableIOOpts::genCommaFileOut(const bfs::path & outFilename, bool header){
	return TableIOOpts(InOptions(), ",", OutOptions(outFilename.string(), ".csv", "comma"), ",", header);
}

TableIOOpts TableIOOpts::genCommaFileIn(const bfs::path & inFilename, bool header){
	return TableIOOpts(InOptions(inFilename), ",", OutOptions("", ".csv", "comma"), ",", header);
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

}


}  // namespace bibseq


