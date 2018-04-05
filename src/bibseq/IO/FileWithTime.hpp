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
 * FileWithTime.hpp
 *
 *  Created on: Feb 15, 2016
 *      Author: nick
 */



#include "bibseq/utils.h"

namespace bibseq {

/**@brief Simple class with boost filesystem path and time modified when last checked
 *
 */
class FileWithTime{
public:
	/**@brief Construct with filepath, fnp must exist will throw otherwise
	 *
	 * @param fnp Filepath name
	 */
	FileWithTime(boost::filesystem::path fnp);
	const boost::filesystem::path fnp_; /**< path to the file*/
private:
	std::chrono::time_point<std::chrono::system_clock> time_ /**< time file was modified last time checked*/;
public:
	/**@brief Check the modification time and return true if it has been modified since last time, update FileWithTime::time_
	 *
	 * @return Whether the file has been modified since the last time it was checked
	 */
	bool checkTime() const;

	/**@brief Set the time to the current last write time of fnp_;
	 *
	 */
	void setTime();
};

class FilesWithTime {
public:

	FilesWithTime(const std::vector<boost::filesystem::path> & files);
	/**@brief a vector of the files with their last write times when the last time they were set
	 *
	 */
	std::vector<FileWithTime> files_;
	/**@brief Check to see if the files have been modified since their last usage
	 *
	 * @return true if any of the files have changed
	 */
	bool needsUpdate() const;
	/**@brief update the last write times of the files
	 *
	 */
	void updateTimes();
};

}  // namespace bibseq
