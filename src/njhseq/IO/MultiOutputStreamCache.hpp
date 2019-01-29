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
 * MultiOutputStreamCache.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: nick
 */

#include "njhseq/IO/MultiOutputStream.hpp"

namespace njhseq {

/**@brief A class to cache reads rather than writing at once to limit io usage
 *
 */

class MultiOutputStreamCache {
public:

	/**@brief create a outputStream with the uid and options
	 *
	 * @param uid the uid for the reader
	 * @param opts The options for the reader
	 */
	void addOutputStream(const std::string & uid, const OutOptions & opts) {
		writers_.addStream(uid, opts);
	}

	/**@brief Add a read to cache for the uid SeqIO
	 *
	 * This will add the read to the cache and if it hits the cache limit
	 * it will write the cache and then clear it
	 *
	 * @param uid the uid of the SeqIO
	 * @param read the read to add
	 */
	void add(const std::string & uid, const std::string & line) {
		writers_.containsStreamThrow(uid);
		cache_[uid].push_back(line);
		++cacheSize_;
		if (cacheSize_ == cacheLimit_) {
			writeCache();
		}
	}
	/**@brief Add a reads to cache for the uid SeqIO
	 *
	 * This will add the reads to the cache and if it hits the cache limit
	 * it will write the cache and then clear it
	 *
	 * @param uid the uid of the SeqIO
	 * @param reads the reads to add
	 */
	void add(const std::string & uid, const std::vector<std::string> & lines) {
		writers_.containsStreamThrow(uid);
		for (const auto & line : lines) {
			cache_[uid].push_back(line);
			++cacheSize_;
			if (cacheSize_ == cacheLimit_) {
				writeCache();
			}
		}
	}
	/**@brief Write cache and then clear it
	 *
	 */
	void writeCache() {
		for (const auto & lines : cache_) {
			if(lines.second.empty()){
				continue;
			}
			//std::cout << __PRETTY_FUNCTION__ << std::endl;
			writers_.openWrite(lines.first, lines.second);
			//std::cout << __PRETTY_FUNCTION__ << std::endl;
		}
		for (auto & reads : cache_) {
			reads.second.clear();
		}
		cacheSize_ = 0;
	}

	/**@brief Close all the output read files for all the readObjectIOs held in MultiOutputStreamCache::writer_
	 *
	 */
	void closeOutAll() {
		writeCache();
		writers_.closeOutAll();
	}

	/**@brief Close all the output read files for all the readObjectIOs held in MultiOutputStreamCache::writer_ for reopening (their option will be changed to append to append to the file now on re-opening)
	 *
	 */
	void closeOutForReopeningAll() {
		writeCache();
		writers_.closeOutForReopeningAll();
	}

	/**@brief Get the current cache limit
	 *
	 * @return
	 */
	uint32_t getCacheLimit() const {
		return cacheLimit_;
	}

	/**@brief Set the cache limit
	 *
	 * @param cacheLimit Set the new cache limit
	 */
	void setCacheLimit(uint32_t cacheLimit) {
		cacheLimit_ = cacheLimit;
	}

	/**@brief Set the file open limit for MultiOutputStreamCache::writer_
	 *
	 * @param fileOpenLimit the new file open limit
	 */
	void setOpenLimit(uint32_t fileOpenLimit){
		writers_.setOpenLimit(fileOpenLimit);
	}

	/**@brief On destruction, close out writers which will also write the cache
	 *
	 */
	~MultiOutputStreamCache(){
		closeOutForReopeningAll();
	}

	void containsReaderThrow(const std::string & uid) const {
		writers_.containsStreamThrow(uid);
	}

	bool containsReader(const std::string & uid) const {
		return writers_.containsStream(uid);
	}

private:

	std::unordered_map<std::string, std::vector<std::string>> cache_; /**< The cache of lines*/
	uint32_t cacheLimit_ = 10000;/**< The cache limit*/
	uint32_t cacheSize_ = 0;/**< The current cache size*/
	MultiOutputStream writers_;/**< The MultiOutputStream responsible for writing the cache*/
};

}  // namespace njhseq



