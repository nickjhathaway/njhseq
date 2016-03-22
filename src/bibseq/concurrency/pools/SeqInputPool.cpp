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
/*
 * SeqIO.cpp
 *
 *  Created on: Jan 17, 2016
 *  Author: nick hathaway
 */

#include "SeqInputPool.hpp"

namespace bibseq {
namespace concurrent {

SeqInputPool::SeqInputPool(const SeqIOOptions & opts, const size_t size) :
		opts_(opts), size_(size), open_(false), closing_(false) {
	for(uint32_t i = 0; i < size_; ++i){
		readers_.emplace_back(opts_);
	}
}

SeqInputPool::SeqInputPool(const SeqInputPool& that) :
		opts_(that.opts_), size_(that.size_), open_(false), closing_(false) {
	for(uint32_t i = 0; i < size_; ++i){
		readers_.emplace_back(opts_);
	}
}

SeqInputPool::~SeqInputPool() {
	std::lock_guard<std::mutex> lock(poolmtx_); // GUARD
	closeSeqFile_nolock();
}

void SeqInputPool::openSeqFile() {

	std::lock_guard<std::mutex> lock(poolmtx_); // GUARD
	closeSeqFile_nolock();
	if (!open_) {
		// try to open all readers
		for (SeqInput& reader : readers_) {
			reader.openIn();
		}
		// all are open, now push them onto the queue
		for (SeqInput& reader : readers_) {
			pushReader(reader);
		}
		open_ = true;
	} else {
		throw SeqInputPoolException("Failed to close BAM reader pool.");
	}
}



void SeqInputPool::closeSeqFile() {
	std::lock_guard<std::mutex> lock(poolmtx_); // GUARD
	closeSeqFile_nolock();
}


void SeqInputPool::closeSeqFile_nolock() {
	if (open_) {
		closing_ = true;
		for (size_t i = 0; i < size_; ++i) {
			PooledSeqReader reader_ptr;
			queue_.waitPop(reader_ptr);
			reader_ptr->closeIn();
		}
	}
	open_ = false;
	closing_ = false;
}

void SeqInputPool::pushReader(SeqInput& reader) {
	if (!closing_) {
		SeqInputPool* p = this;
		queue_.push(
				std::shared_ptr<SeqInput>(&reader,
						[p](SeqInput* reader)
						{
							p->pushReader(*reader);
						}));
	}
}

PooledSeqReader SeqInputPool::popReader() {
	PooledSeqReader reader_ptr;
	queue_.waitPop(reader_ptr);
	return reader_ptr;
}

const SeqIOOptions& SeqInputPool::getOpts() const {
	return opts_;
}

} // namespace concurrent
} // namespace bamcuts

