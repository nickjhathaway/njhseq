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
 * DataFileReader.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "../BioDataObject/DataFileReader.hpp"

namespace bibseq {

DataFileReader::DataFileReader(const std::string & filename):filename_(filename){
	file_.open(filename, std::ios::in);
	checkFileThrow();
}

DataFileReader::operator bool() const {
	return file_.is_open();
}
DataFileReader::~DataFileReader(){
	file_.close();
}


void DataFileReader::checkFileThrow() const {
	if (!file_.is_open()) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << ", error in opening: "
				<< bib::bashCT::boldBlack(filename_) << std::endl;
		if (!bib::files::bfs::exists(filename_)) {
			ss << "File: " << bib::bashCT::boldBlack(filename_) << " doesn't exist"
					<< std::endl;
		}
		throw std::runtime_error { bib::bashCT::boldRed(ss.str()) };
	}
}

} /* namespace bibseq */