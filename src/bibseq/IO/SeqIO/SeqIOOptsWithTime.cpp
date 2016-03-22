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
 * SeqIOOptsWithTime.cpp
 *
 *  Created on: Feb 20, 2016
 *      Author: nick
 */


#include "SeqIOOptsWithTime.hpp"
namespace bibseq {




SeqIOOptsWithTime::SeqIOOptsWithTime(const SeqIOOptions & opts) :
		opts_(opts) {

}

SeqIOOptsWithTime::SeqIOOptsWithTime(const SeqIOOptsWithTime& other) :
		opts_(other.opts_) {

}

bool SeqIOOptsWithTime::outDated() const {
	bool needsUpdate = true;
	if (opts_.inExists()) {
		auto inTime = bib::files::last_write_time(opts_.firstName_);
		if (time_ == inTime) {
			needsUpdate = false;
		}
	}
	return needsUpdate;
}

}  // namespace bibseq

