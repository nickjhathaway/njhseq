#pragma once
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

	explicit operator bool() const ;

	template<typename DATA>
	bool readNextRecord(DATA & record){
		std::string line = "";
		if(std::getline(file_, line)){
			record = DATA(line);
			return true;
		}else{
			return false;
		}
	}

	virtual ~DataFileReader();


	// no-copy
	DataFileReader(const DataFileReader& other) = delete;
	DataFileReader() = delete;
	DataFileReader& operator=(const DataFileReader& other) = delete;
};

} /* namespace bibseq */



