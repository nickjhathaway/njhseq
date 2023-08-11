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
 * MultiSeqOutCache.hpp
 *
 *  Created on: Feb 17, 2016
 *      Author: nick
 */

#include "njhseq/IO/SeqIO/MultiSeqIO.hpp"
//#include <shared_mutex>

namespace njhseq {

/**@brief A class to cache reads rather than writing at once to limit io usage
 *
 */
template<typename T>
class MultiSeqOutCache {
public:

	/**@brief create a reader with the uid and options
	 *
	 * @param uid the uid for the reader
	 * @param opts The options for the reader
	 */
	void addReader(const std::string & uid, const SeqIOOptions & opts) {
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		writers_.addReader(uid, opts);
	}

	/**@brief Add a read to cache for the uid SeqIO
	 *
	 * This will add the read to the cache and if it hits the cache limit
	 * it will write the cache and then clear it
	 *
	 * @param uid the uid of the SeqIO
	 * @param read the read to add
	 */
	void add(const std::string & uid, const T & read) {
		writers_.containsReaderThrow(uid);
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		cache_[uid].push_back(read);
		++cacheSize_;
		if (cacheSize_ >= cacheLimit_) {
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
	void add(const std::string & uid, const std::vector<T> & reads) {
		writers_.containsReaderThrow(uid);
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		for (const auto & read : reads) {
			cache_[uid].push_back(read);
			++cacheSize_;
			if (cacheSize_ >= cacheLimit_) {
				writeCache();
			}
		}
	}
	/**@brief Write cache and then clear it
	 *
	 */
	void writeCache() {
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		writeCacheLockFree();
	}

	void writeCacheLockFree() {
		for (const auto & reads : cache_) {
			if(reads.second.empty()){
				continue;
			}
			//std::cout << __PRETTY_FUNCTION__ << std::endl;
			writers_.openWrite(reads.first, reads.second);
			//std::cout << __PRETTY_FUNCTION__ << std::endl;
		}
		for (auto & reads : cache_) {
			reads.second.clear();
		}
		cacheSize_ = 0;
	}

	/**@brief Close all the output read files for all the readObjectIOs held in MultiSeqOutCache::writer_
	 *
	 */
	void closeOutAll() {
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		writeCache();
		writers_.closeOutAll();
	}

	/**@brief Close all the output read files for all the readObjectIOs held in MultiSeqOutCache::writer_ for reopening (their option will be changed to append to append to the file now on re-opening)
	 *
	 */
	void closeOutForReopeningAll() {
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		writeCache();
		writers_.closeOutForReopeningAll();
	}

	/**@brief Get the current cache limit
	 *
	 * @return the cache limit
	 */
	uint32_t getCacheLimit() const {
		//std::shared_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		return cacheLimit_;
	}

	/**@brief Set the cache limit
	 *
	 * @param cacheLimit Set the new cache limit
	 */
	void setCacheLimit(uint32_t cacheLimit) {
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		cacheLimit_ = cacheLimit;
	}

	/**@brief Set the file open limit for MultiSeqOutCache::writer_
	 *
	 * @param fileOpenLimit the new file open limit
	 */
	void setOpenLimit(uint32_t fileOpenLimit){
		//std::unique_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		writers_.setOpenLimit(fileOpenLimit);
	}

	/**@brief On destruction, close out writers which will also write the cache
	 *
	 */
	~MultiSeqOutCache(){
		closeOutForReopeningAll();
	}

	void containsReaderThrow(const std::string & uid) const {
		//std::shared_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		writers_.containsReaderThrow(uid);
	}

	bool containsReader(const std::string & uid) const {
		//std::shared_lock<std::shared_timed_mutex> lock(mut_,std::defer_lock);
		return writers_.containsReader(uid);
	}

private:

	std::unordered_map<std::string, std::vector<T>> cache_; /**< The cache of reads*/
	uint32_t cacheLimit_ = 50000;/**< The cache limit*/
	uint32_t cacheSize_ = 0;/**< The current cache size*/
	MultiSeqIO writers_;/**< The MultiSeqIO responsible for writing the cache*/

	//std::shared_timed_mutex mut_; /**< to make thread safe*/
};

}  // namespace njhseq



