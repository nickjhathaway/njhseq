//
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
/*
 * SeqIOOptsWithTime.cpp
 *
 *  Created on: Feb 20, 2016
 *      Author: nick
 */


#include "SeqIOOptsWithTime.hpp"
namespace njhseq {

std::chrono::time_point<std::chrono::system_clock> SeqIOOptsWithTime::getTime() const{
	return time_;
}

void SeqIOOptsWithTime::setTime(const std::chrono::time_point<std::chrono::system_clock> & time){
	time_ = time;
}



SeqIOOptsWithTime::SeqIOOptsWithTime(const SeqIOOptions & opts) :
		opts_(opts) {

}

SeqIOOptsWithTime::SeqIOOptsWithTime(const SeqIOOptsWithTime& other) :
		opts_(other.opts_) {

}

SeqIOOptsWithTime::SeqIOOptsWithTime(const SeqIOOptsWithTime&& other) :
		opts_(other.opts_) {

}

bool SeqIOOptsWithTime::outDated() const {
	bool needsUpdate = true;
	if (opts_.inExists()) {
		auto inTime = njh::files::last_write_time(opts_.firstName_);
		if (time_ == inTime) {
			needsUpdate = false;
		}
	}
	return needsUpdate;
}

}  // namespace njhseq

