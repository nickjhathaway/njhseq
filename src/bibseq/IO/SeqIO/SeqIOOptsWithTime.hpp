#pragma once
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

	std::chrono::time_point<std::chrono::system_clock> getTime() const;

	void setTime(const std::chrono::time_point<std::chrono::system_clock> & time);

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

	SeqIOOptsWithTime(const SeqIOOptsWithTime&& other);

	/**@brief Read in Whole File
	 *
	 * @return A vector of read objects from the seq io options
	 */
	template<typename T>
	std::vector<T> get() {
		time_ = bib::files::last_write_time(opts_.firstName_);
		SeqInput reader(opts_);
		return reader.readAllReads<T>();
	}

	/**@brief Get a number of reads from a read position (not file position) in the file
	 *
	 * @param pos The nth read to start getting reads from
	 * @param number How many reads after pos to get
	 * @return A vector of reads
	 */
	template<typename T>
	std::vector<T> get(size_t pos, size_t number) {
		time_ = bib::files::last_write_time(opts_.firstName_);
		SeqInput reader(opts_);
		return reader.getReads<T>(pos, number);
	}

	/**@brief Get a vector of pointers of type T from the file
	 *
	 * @return A vector of shared pointers of read objects from the seq io options
	 */
	template<typename T>
	std::vector<std::shared_ptr<T>> getPtrs() {
		time_ = bib::files::last_write_time(opts_.firstName_);
		SeqInput reader(opts_);
		return reader.readAllReadsPtrs<T>();
	}

	/**@brief Get a number of reads from a read position (not file position) in the file as a vector of shared pointers
	 *
	 * @param pos The nth read to start getting reads from
	 * @param number How many reads after pos to get
	 * @return A vector of shared pointers of reads
	 */
	template<typename T>
	std::vector<std::shared_ptr<T>> getPtrs(size_t pos, size_t number) {
		time_ = bib::files::last_write_time(opts_.firstName_);
		SeqInput reader(opts_);
		return reader.getReadsPtrs<T>(pos, number);
	}

	/**@brief Check to see if the file has changed since last get, will not modify the SeqIOOptsWithTime::time_
	 *
	 * @return Whether the file has changed since the last read
	 */
	bool outDated() const;
};



}  // namespace bibseq

