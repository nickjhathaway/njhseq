#pragma once
/*
 * TableCache.hpp
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
#include "bibseq/objects/dataContainers/tables/table.hpp"

namespace bibseq {

class TableCache {
private:


	table tab_; /**< the table content */
	std::chrono::time_point<std::chrono::system_clock> time_; /**< time last read */

	/**@brief Load the content of the file and log the time the file was last edited
	 *
	 */
	void load() ;

public:
	const TableIOOpts opts_; /**< the table input options*/

	/**@brief constructor with the content of the of the file given by opts
	 *
	 * @param opts
	 */
	TableCache(const TableIOOpts & opts);

	/**@brief construct the table with the the input table and options for writing
	 *
	 * @param opts The table input/output options
	 * @param inputTable The input table
	 */
	TableCache(const TableIOOpts & opts, const table & inputTable);

	/**@brief Copy constructor
	 *
	 * @param other table cache
	 */
	TableCache(const TableCache& other);

	/**@brief Get the content of the file and update as needed
	 *
	 * @return the current content of the file
	 */
	const table& get() ;

	/**@brief Check to see if the file has changed since
	 *
	 * @return
	 */
	bool needsUpdate() const ;

	/**@brief Update file and return whether the file had to be updated
	 *
	 * @return Whether the file needed to be updated
	 */
	bool update();

	/**@brief Write the table, and make sure the input options are what were written
	 *
	 */
	void clearTable();

};

}  // namespace bibseq



