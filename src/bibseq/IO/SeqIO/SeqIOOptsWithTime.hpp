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
 * SeqIOOptsWithTime.hpp
 *
 *  Created on: Feb 20, 2016
 *      Author: nick
 */

#include "bibseq/IO/SeqIO/SeqInput.hpp"

namespace bibseq {


/**@brief SeqIOOpts with a time point member to keep track if the input file has changed
 *
 */
class SeqIOOptsWithTime {
private:

	std::chrono::time_point<std::chrono::system_clock> time_; /**< time last read */
public:

	const SeqIOOptions opts_; /**< the table input options*/

	/**@brief constructor with the content of the of the file given by opts
	 *
	 * @param opts
	 */
	SeqIOOptsWithTime(const SeqIOOptions & opts);

	/**@brief Copy constructor, just copy the other options, don't copy time
	 *
	 * @param other table cache
	 */
	SeqIOOptsWithTime(const SeqIOOptsWithTime& other);
	template<typename T>
	std::vector<T> get() {
		time_ = bib::files::last_write_time(opts_.firstName_);
		SeqInput reader(opts_);
		return reader.readAllReads<T>();
	}

	/**@brief Get a vector of pointers of type T from the file
	 *
	 * @return
	 */
	template<typename T>
	std::vector<std::shared_ptr<T>> getPtrs() {
		time_ = bib::files::last_write_time(opts_.firstName_);
		SeqInput reader(opts_);
		return reader.readAllReadsPtrs<T>();
	}

	/**@brief Check to see if the file has changed since last get, will not modify the SeqIOOptsWithTime::time_
	 *
	 * @return
	 */
	bool outDated() const;
};



}  // namespace bibseq

