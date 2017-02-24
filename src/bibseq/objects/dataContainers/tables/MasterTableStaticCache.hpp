#pragma once
/*
 * MasterTableStaticCache.hpp
 *
 *  Created on: Feb 13, 2016
 *      Author: nick
 */

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

#include "bibseq/objects/dataContainers/tables/table.hpp"
#include "bibseq/IO/FileWithTime.hpp"

namespace bibseq {

class MasterTableStaticCache {
public:
	MasterTableStaticCache(const TableIOOpts & opts,
			const std::vector<boost::filesystem::path> & files);

	TableIOOpts opts_;
private:
	FilesWithTime files_;
	std::chrono::time_point<std::chrono::system_clock> updateTime_;
	table masterTab_;
private:
	/**@brief checks to see if MasterTableStaticCache::masterTab_ needs to be updated
	 * @return true if the master table needs to be updated
	 */
	bool needsUpdate();
public:
	/**@brief loads in tables if it needs to
	 * @return whether the tables needed to be loaded in
	 */
	bool loadTabs();

	/**@brief get an up to date table
	 * @return a constant reference to the master table
	 */
	const table& get();

	/**@brief write the table,
	 * loads in the input files if it needs to and will write only if the out is out of date
	 */
	void writeTab();

	/**@brief write the table,gzipped
	 * loads in the input files if it needs to and will write only if the out is out of date
	 */
	void writeTabGz();

	/**@brief Checks to see if the out file is newer than input files and if it exists
	 * not safe if overwriting the out file and the out file wasn't created from these input files
	 * as it could be newer than the input files but not contain their contents
	 * @return whether the out file is up to date
	 */
	bool outUpToDate() const;

};

}  // namespace bibseq





