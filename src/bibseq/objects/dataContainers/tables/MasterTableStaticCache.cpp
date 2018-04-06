/*
 * MasterTableStaticCache.cpp
 *
 *  Created on: Feb 13, 2016
 *      Author: nick
 */

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
#include "MasterTableStaticCache.hpp"

namespace bibseq {

MasterTableStaticCache::MasterTableStaticCache(const TableIOOpts & opts,
		const std::vector<bfs::path> & files,
		bool fill) :
		opts_(opts),fill_(fill), files_(files) {
	bool outOutDated = outUpToDate();
	if (outOutDated) {
		loadTabs();
	} else {
		masterTab_ = table(
				bib::appendAsNeededRet(opts_.out_.outFilename_.string(),
						opts_.out_.outExtention_), opts.outDelim_, opts_.hasHeader_);
	}
}

bool MasterTableStaticCache::outUpToDate() const {
	auto outFilename = bib::appendAsNeededRet(opts_.out_.outFilename_.string(),
			opts_.out_.outExtention_);
	bool outOutDated = false;
	if (bfs::exists(outFilename)) {
		auto outTime = bib::files::last_write_time(outFilename);
		for (const auto & file : files_.files_) {
			if (outTime < bib::files::last_write_time(file.fnp_)) {
				outOutDated = true;
				break;
			}
		}
	} else {
		outOutDated = true;
	}
	return outOutDated;
}

bool MasterTableStaticCache::needsUpdate() {
	return masterTab_.content_.empty() || files_.needsUpdate();
}

bool MasterTableStaticCache::loadTabs(){
	if(needsUpdate()){
		files_.updateTimes();
		masterTab_ = table();
		for (const auto & file : files_.files_) {
			table fileTab(file.fnp_.string(), opts_.inDelim_,
					opts_.hasHeader_);
			if (masterTab_.content_.empty()) {
				masterTab_ = fileTab;
			} else {
				masterTab_.rbind(fileTab, fill_);
			}
		}
		updateTime_ = std::chrono::system_clock::now();
		return true;
	}
	return false;
}

const table& MasterTableStaticCache::get(){
	loadTabs();
	return masterTab_;
}

void MasterTableStaticCache::writeTab(){
	bool outOutDated = outUpToDate();
	if(outOutDated){
		loadTabs();
		masterTab_.outPutContents(opts_);
	}
}

void MasterTableStaticCache::writeTabGz(){
	writeTab();

	bfs::path gzPath = opts_.out_.outName().string() + ".gz";
	if(!bfs::exists(gzPath) || bib::files::firstFileIsOlder(gzPath, opts_.out_.outName())){
		auto gzOpts = IoOptions(InOptions(opts_.out_.outName()), OutOptions(gzPath));
		gzOpts.out_.overWriteFile_ = opts_.out_.overWriteFile_;
		gzZipFile(gzOpts);
	}

}

}  // namespace bibseq


