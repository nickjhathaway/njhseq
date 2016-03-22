#pragma once
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
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"

namespace bibseq {

class DataFileReader {
public:

	DataFileReader(const std::string & filename);

	std::ifstream file_;
	std::string filename_;

	explicit operator bool() const;

	template<typename DATA>
	bool readNextRecord(DATA & record) {
		std::string line = "";
		if (bib::files::crossPlatGetline(file_, line)) {
			while ((line.front() == '#' || line.empty()) && !file_.eof()) {
				bib::files::crossPlatGetline(file_, line);
			}
			if (line.front() == '#' || line.empty()) {
				return false;
			} else {
				record = DATA(line);
				return true;
			}
		} else {
			return false;
		}
	}

	template<typename DATA>
	std::shared_ptr<DATA> readNextRecord() {
		std::string line = "";
		if (bib::files::crossPlatGetline(file_, line)) {
			while ((line.front() == '#' || line.empty()) && !file_.eof()) {
				bib::files::crossPlatGetline(file_, line);
			}
			if (line.front() == '#' || line.empty()) {
				return nullptr;
			} else {
				return std::make_shared<DATA>(line);
			}
		} else {
			return nullptr;
		}
	}

	void checkFileThrow() const;
	virtual ~DataFileReader();

	// no-copy
	DataFileReader(const DataFileReader& other) = delete;
	DataFileReader() = delete;
	DataFileReader& operator=(const DataFileReader& other) = delete;
};

} /* namespace bibseq */

