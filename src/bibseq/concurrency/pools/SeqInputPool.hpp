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
 * SeqIO.hpp
 *
 *  Created on: Jan 17, 2016
 *  Author: nick hathaway
 */



#pragma once


#include <bibcpp.h>
#include "bibseq/concurrency/ConcurrentQueue.hpp"
#include "bibseq/IO/SeqIO/SeqIO.hpp"



namespace bibseq {
namespace concurrent {

typedef std::shared_ptr<SeqInput> PooledSeqReader;

// Pool of bibseq::SeqInput
class SeqInputPool {
private:
	SeqIOOptions opts_;	// Sequence file to operate on
	const uint32_t size_; // pool size
	bool open_; // pool status (open/closed)
	std::vector<SeqInput> readers_; // vector of BamReaders
	concurrent::ConcurrentQueue<PooledSeqReader> queue_; // "ready-queue"
	std::mutex poolmtx_;
	bool closing_; // state indicator for closing down

	// close all readers (no locking)
	void closeSeqFile_nolock();

	// push reader onto the ready-queue
	void pushReader(SeqInput& reader);

public:
	SeqInputPool(const SeqIOOptions & opts, const size_t size);
	SeqInputPool(const SeqInputPool& that);
	virtual ~SeqInputPool();

	// open seq file
	void openSeqFile();

	// close seq file
	void closeSeqFile();

	// request a reader from the pool (wait on queue)
	PooledSeqReader popReader();

	// accessor for file_
	const SeqIOOptions& getOpts() const;
};

class SeqInputPoolException: public bib::err::Exception {
public:
	SeqInputPoolException(const std::string s) :
			Exception(s) {

	}
	SeqInputPoolException() :
			Exception("General SeqInputPool exception.") {
	}
};

} // namepsace concurrent
} // namespace bibseq
