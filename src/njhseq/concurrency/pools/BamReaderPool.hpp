#pragma once
/*
 * BamReaderPool.hpp
 *
 *  Created on: Jan 7, 2016
 */

/*
 Copyright 2014 Arjan van der Velde, vandervelde.ag [at] gmail
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 http://www.apache.org/licenses/LICENSE-2.0
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

#pragma once

#include <vector>
#include <mutex>
#include <boost/filesystem.hpp>

#include <api/BamReader.h> // bamtools
#include <njhcpp.h>
#include "njhseq/utils.h"
#include "njhseq/concurrency/ConcurrentQueue.hpp"




namespace njhseq {
namespace concurrent {


typedef std::shared_ptr<BamTools::BamReader> PooledReader;

// Pool of BamTools::BamReaders
class BamReaderPool {
private:
	bfs::path file_;	// BAM file to operate on
	const uint32_t size_; // pool size
	bool open_; // pool status (open/closed)
	std::vector<BamTools::BamReader> readers_; // vector of BamReaders
	concurrent::ConcurrentQueue<PooledReader> queue_; // "ready-queue"
	std::mutex poolmtx_;
	bool closing_; // state indicator for closing down

	// close all readers (no locking)
	void closeBamFile_nolock();

	// push reader onto the ready-queue
	void pushReader(BamTools::BamReader& reader);

public:
	BamReaderPool(const bfs::path& file, const size_t size);
	BamReaderPool(const size_t size);
	BamReaderPool(const BamReaderPool& that);
	virtual ~BamReaderPool();

	// open BAM file (file_)
	void openBamFile();

	// open BAM file
	void openBamFile(const bfs::path& file);

	// close BAM file
	void closeBamFile();

	// request a reader from the pool (wait on queue)
	PooledReader popReader();

	// accessor for file_
	const bfs::path& getFile() const;
};

class BamReaderPoolException: public njh::err::Exception {
public:
	BamReaderPoolException(const std::string s) :
			Exception(s) {

	}
	BamReaderPoolException() :
			Exception("General BamReaderPool exception.") {
	}
};

} // namepsace concurrent
} // namespace njhseq
