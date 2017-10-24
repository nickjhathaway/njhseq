/*
 * TableCache.cpp
 *
 *  Created on: Feb 13, 2016
 *      Author: nick
 */

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
#include "TableCache.hpp"

namespace bibseq {


void TableCache::load() {
	tab_ = table(opts_);
	time_ = bib::files::last_write_time(opts_.in_.inFilename_);
}

TableCache::TableCache(const TableIOOpts & opts) :
	opts_(opts) {
	//load();
}

TableCache::TableCache(const TableIOOpts & opts, const table & inputTable) :
		tab_(inputTable), opts_(opts) {

}


TableCache::TableCache(const TableCache& other) :
	opts_(other.opts_) {
	//load();
}
const table& TableCache::get() {
	update();
	return tab_;
}


bool TableCache::needsUpdate() const {
	return tab_.empty() || time_ != bib::files::last_write_time(opts_.in_.inFilename_);
}

bool TableCache::update(){
	if (needsUpdate()) {
		load();
		return true;
	}
	return false;
}


void TableCache::clearTable(){
	if(!opts_.in_.inExists()){
		//if the TableCache was created with an input table rather than by reading in a table
		//check to see if the infile exist before clearing the table, if it doesn't, write the
		//table so when TableChache::load() reads in it will read the correct table
		auto optsCopy = opts_;
		optsCopy.outDelim_ = opts_.inDelim_;
		optsCopy.hasHeader_ = opts_.hasHeader_;
		optsCopy.out_.outFilename_ = opts_.in_.inFilename_;
		optsCopy.out_.outExtention_ = opts_.in_.inFilename_.string();
		tab_.outPutContents(optsCopy);
	}
	tab_ = table();
}



}  // namespace bibseq

