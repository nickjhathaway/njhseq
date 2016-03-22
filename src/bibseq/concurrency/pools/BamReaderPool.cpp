/*
 * BamReaderPool.cpp
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

#include "BamReaderPool.hpp"

namespace bibseq {
namespace concurrent {

BamReaderPool::BamReaderPool(const bfs::path& file, const size_t size) :
		file_(file), size_(size), open_(false), closing_(false) {
	readers_.resize(size_);
}

BamReaderPool::BamReaderPool(const size_t size) :
		file_(bfs::path()), size_(size), open_(false), closing_(false) {
	readers_.resize(size_);
}

BamReaderPool::BamReaderPool(const BamReaderPool& that) :
		file_(that.file_), size_(that.size_), open_(false), closing_(false) {
	readers_.resize(size_);
}

BamReaderPool::~BamReaderPool() {
	std::lock_guard<std::mutex> lock(poolmtx_); // GUARD
	closeBamFile_nolock();
}

void BamReaderPool::openBamFile() {
	openBamFile(file_);
}

void BamReaderPool::openBamFile(const bfs::path& file) {
	if (file.empty()) {
		throw BamReaderPoolException("No BAM file specified.");
	}
	std::lock_guard<std::mutex> lock(poolmtx_); // GUARD
	closeBamFile_nolock();
	if (!open_) {
		// try to open all readers
		for (BamTools::BamReader& reader : readers_) {
			if (!reader.Open(file.c_str())) {
				throw BamReaderPoolException(
						"Unable to open BAM file " + file.string() + ":\n\t"
								+ reader.GetErrorString());
			}
		}
		// try to open index
		for (BamTools::BamReader& reader : readers_) {
			if (!reader.LocateIndex()) {
				throw BamReaderPoolException(
						"Unable to open BAM index file " + file.string() + ":\n\t"
								+ reader.GetErrorString());
			}
		}
		// all are open, now push them onto the queue
		for (BamTools::BamReader& reader : readers_) {
			pushReader(reader);
		}
		file_ = file;
		open_ = true;
	} else {
		throw BamReaderPoolException("Failed to close BAM reader pool.");
	}
}

void BamReaderPool::closeBamFile() {
	std::lock_guard<std::mutex> lock(poolmtx_); // GUARD
	closeBamFile_nolock();
}

void BamReaderPool::closeBamFile_nolock() {
	if (open_) {
		closing_ = true;
		for (size_t i = 0; i < size_; ++i) {
			PooledReader reader_ptr;
			queue_.waitPop(reader_ptr);
			reader_ptr->Close();
		}
	}
	open_ = false;
	closing_ = false;
}

void BamReaderPool::pushReader(BamTools::BamReader& reader) {
	if (!closing_) {
		BamReaderPool* p = this;
		queue_.push(
				std::shared_ptr<BamTools::BamReader>(&reader,
						[p](BamTools::BamReader* reader)
						{
							p->pushReader(*reader);
						}));
	}
}

PooledReader BamReaderPool::popReader() {
	PooledReader reader_ptr;
	queue_.waitPop(reader_ptr);
	return reader_ptr;
}

const bfs::path& BamReaderPool::getFile() const {
	return file_;
}

} // namespace concurrent
} // namespace bamcuts

