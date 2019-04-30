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
 * MultiSeqIO.hpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */

#include "njhseq/IO/SeqIO/SeqIO.hpp"


namespace njhseq {



/**@brief A class to hold multiple njhseq::readObjectIOOpt and to help with keeping a limit on the number of files opned at once
 *
 */
class MultiSeqIO {
private:
	std::unordered_map<std::string, std::unique_ptr<SeqIO>> readIos_; /**< Map to hold all the readObjectIOs */
public:
	/**@brief Check to see if reader with uid exists
	 *
	 * @param uid The uid to check for
	 * @return bool whether the reader exists
	 */
	bool containsReader(const std::string & uid) const;
	/**@brief create a reader with the uid and options
	 *
	 * @param uid the uid for the reader
	 * @param opts The options for the reader
	 */
	void addReader(const std::string & uid, const SeqIOOptions & opts);

	/**@brief Read the next read from the given reader
	 *
	 * @param uid the uid of the reader
	 * @param read the readObject to store the info for the next reader
	 * @return bool whether there was a next read to read
	 */
	template<typename T>
	bool readNextRead(const std::string & uid, T & read) {
		containsReaderThrow(uid);
		std::lock_guard<std::mutex> lock(readIos_.at(uid)->mut_);
		return readIos_.at(uid)->in_.readNextRead(read);
	}

	/**@brief Read the next read from the given reader and return a shared pointer to it's location
	 *
	 * @param uid the uid of the reader
	 * @return a shared pointer to the info for the next read that was read, will be nullptr if no read to read
	 */
	template<typename T>
	std::shared_ptr<T> readNextRead(const std::string & uid) {
		containsReaderThrow(uid);
		std::lock_guard<std::mutex> lock(readIos_.at(uid)->mut_);
		return readIos_.at(uid)->in_.readNextRead<T>();
	}

	/**@brief Close all the input read files for all the readObjectIOs held in MultiReadObjectIOOpt::readIos_
	 *
	 */
	void closeInAll();

	/**@brief Close all the output read files for all the readObjectIOs held in MultiReadObjectIOOpt::readIos_
	 *
	 */
	void closeOutAll();

	/**@brief Close all the output read files for all the readObjectIOs for reopening (their option will be changed to append to append to the file now on reopning) held in MultiReadObjectIOOpt::readIos_
	 *
	 */
	void closeOutForReopeningAll();

	/**@brief Write the given seq object for reader with uid, will fail if reader outputs has not been opened
	 *
	 * @param uid The uid of the reader to write to
	 * @param read The seqInfo to write
	 */
	template<typename T>
	void write(const std::string & uid, const T & seq) {
		containsReaderThrow(uid);
		std::lock_guard<std::mutex> lock(readIos_.at(uid)->mut_);
		readIos_.at(uid)->out_.write(seq);
	}

	/**@brief Write the given seqInfo for reader with uid, will open the the reader's outputs if they are not currently opened
	 *
	 * Will keep track of how many files are open and keep only MultiReadObjectIOOpt::outOpenLimit_ number of files open
	 * @param uid The uid of the reader to write to
	 * @param read The seq object to write
	 */
	template<typename T>
	void openWrite(const std::string & uid, const T & seq) {
		openOut(uid);
		std::lock_guard<std::mutex> lock(readIos_.at(uid)->mut_);

		auto readIO = readIos_.find(uid);

		readIO->second->out_.writeNoCheck(seq);
	}

	/**@brief Open uid writer if it isn't open and update priority queue
	 *
	 * @param uid The uid of the writer
	 */
	void openOut(const std::string & uid);

	/**@brief Write the given line for reader with uid, will open the the reader's outputs if they are not currently opened
	 *
	 * Will keep track of how many files are open and keep only MultiReadObjectIOOpt::outOpenLimit_ number of files open
	 * @param uid The uid of the reader to write to
	 * @param line The line to write
	 */
	void openWrite(const std::string & uid, const std::string & line);

	/**@brief
	 *
	 * @param uid Uid of the reader
	 * @param read the sff read to write the flow for
	 */
	void openWriteFlow(const std::string & uid, const sffObject & read);

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

	template<typename T>
	friend class MultiSeqOutCache;
private:
	/**@brief Check for the SeqIO and throw an exception if it isn't found
	 *
	 * @param uid the uid of the SeqIO
	 */
	void containsReaderThrow(const std::string & uid)const;

#if defined( __APPLE__ ) || defined( __APPLE_CC__ ) || defined( macintosh ) || defined( __MACH__ )
	uint32_t outOpenLimit_ = 200; /**< The maximum number of files to be kept open at one time*/
#else
	uint32_t outOpenLimit_ = 968; /**< The maximum number of files to be kept open at one time */
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

