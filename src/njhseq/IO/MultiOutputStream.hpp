#pragma once
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
 * MultiOutputStream.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: nick
 */

#include "njhseq/IO/OutputStreamWrap.hpp"


namespace njhseq {



/**@brief A class to hold multiple OutputStreamWraps and to help with keeping a limit on the number of files opened at once
 *
 */
class MultiOutputStream {
private:
	std::unordered_map<std::string, std::unique_ptr<OutputStreamWrap>> outStreams_; /**< Map to hold all the OutputStreamWrap */
public:
	/**@brief Check to see if stream with uid exists
	 *
	 * @param uid The uid to check for
	 * @return bool whether the stream exists
	 */
	bool containsStream(const std::string & uid) const;
	/**@brief create a stream with the uid and options
	 *
	 * @param uid the uid for the stream
	 * @param opts The options for the stream
	 */
	void addStream(const std::string & uid, const OutOptions & opts);


	/**@brief Close all the output read files for all the OutputStreamWraps held in MultiOutputStream::outStreams_
	 *
	 */
	void closeOutAll();

	/**@brief Close all the output read files for all the OutputStreamWraps for reopening (their option will be changed to append to append to the file now on reopning) held in MultiOutputStream::outStreams_
	 *
	 */
	void closeOutForReopeningAll();

	/**@brief Write the given line for stream with uid, will fail if stream outputs has not been opened
	 *
	 * @param uid The uid of the stream to write to
	 * @param line The line to write
	 * @param flush whether or not to flush the stream after the write
	 */

	void write(const std::string & uid, const std::string & line, bool flush = false);

	/**@brief Write the given seqInfo for stream with uid, will open the the stream's outputs if they are not currently opened
	 *
	 * Will keep track of how many files are open and keep only MultiOutputStream::outOpenLimit_ number of files open
	 * @param uid The uid of the stream to write to
	 * @param line The line to write
	 */

	void openWrite(const std::string & uid, const std::string & line, bool flush = false);

	void openWrite(const std::string & uid, const std::vector<std::string> & lines, bool flush = false);


	/**@brief Open uid writer if it isn't open and update priority queue
	 *
	 * @param uid The uid of the writer
	 */
	void openOut(const std::string & uid);




	/**@brief get the limit for open file numbers
	 *
	 * @return the limit for open files
	 */
	uint32_t getOpenLimit() const;

	/**@brief Set the number of files allowed open
	 *
	 * @param limit The new limit for files to be kept open
	 */
	void setOpenLimit(uint32_t limit);

	/**@brief Check for the OutputStreamWrap and throw an exception if it isn't found
	 *
	 * @param uid the uid of the OutputStreamWrap
	 */
	void containsStreamThrow(const std::string & uid)const;

	friend class MiltiOutpuStreamCache;
private:


#if defined( __APPLE__ ) || defined( __APPLE_CC__ ) || defined( macintosh ) || defined( __MACH__ )
	uint32_t outOpenLimit_ = 200; /**< The maximum number of files to be kept open at one time*/
#else
	uint32_t outOpenLimit_ = 1000; /**< The maximum number of files to be kept open at one time */
#endif

	uint32_t outCurrentlyOpen_ = 0; /**< The number of files currently open*/
	std::deque<std::string> outsOpen_; /**< A deque to keep track of which files to close next*/
	std::unordered_map<std::string, uint32_t> openPriorityCounts_; /**< A map to keep track of how times a uid is in deque to give priority to keep files that were more recently written to*/

	/**@brief Determine the next file that shoudl bed closed and close it
	 *
	 */
	void closeNext();

	std::mutex mut_;/**< mutex to lock the class*/

};



}  // namespace njhseq

