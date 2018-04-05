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
/*
 * FileWithTime.cpp
 *
 *  Created on: Feb 15, 2016
 *      Author: nick
 */

#include "FileWithTime.hpp"

namespace bibseq {


FileWithTime::FileWithTime(boost::filesystem::path fnp):fnp_(fnp){
	if(!boost::filesystem::exists(fnp_)){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << " file "
			 <<	bib::bashCT::red << fnp_ << " doesn't exist"
			 << bib::bashCT::reset << std::endl;
		throw std::runtime_error{ss.str()};
	}
	setTime();
}

void FileWithTime::setTime(){
	time_ = bib::files::last_write_time(fnp_);
}

bool FileWithTime::checkTime() const{
	if(bib::files::last_write_time(fnp_) > time_){
		return true;
	}
	return false;
}

FilesWithTime::FilesWithTime(const std::vector<boost::filesystem::path> & files){
	for(const auto & fnp : files){
		files_.emplace_back(fnp);
	}
}

bool FilesWithTime::needsUpdate() const {
	bool updateNeeded = false;
	for (const auto & file : files_) {
		if (file.checkTime()) {
			updateNeeded = true;
		}
	}
	return updateNeeded;
}

void FilesWithTime::updateTimes() {
	for (auto & file : files_) {
		file.setTime();
	}
}

}  // namespace bibseq


