#pragma once
/*
 * BamReaderPool.hpp
 *
 *  Created on: Jan 7, 2016
 *      Author: nick
 */

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

#pragma once

#include <vector>
#include <mutex>

#include "njhseq/alignment/aligner.h"

#include "njhseq/concurrency/ConcurrentQueue.hpp"

namespace njhseq {
namespace concurrent {

typedef std::shared_ptr<aligner> PooledAligner;

// Pool of BamTools::BamReaders
class AlignerPool {
private:
	//aligner specifics
	uint64_t startingMaxLen_;
	gapScoringParameters gapInfo_;
	substituteMatrix scoring_;

	const uint32_t size_; // pool size

	std::vector<aligner> aligners_; // vector of aligners
	concurrent::ConcurrentQueue<PooledAligner> queue_; // "ready-queue"
	std::mutex poolmtx_;


	void destoryAlignersNoLock();

	// push aligner onto the ready-queue
	void pushAligner(aligner& alignerObj);
	bool closing_ = false;
public:
	AlignerPool(uint64_t startingMaxLen, const gapScoringParameters & gapInfo,
			const substituteMatrix & scoring, const size_t size);
	AlignerPool(const aligner & alignerObj, const size_t size);
	AlignerPool(const size_t size);
	AlignerPool(const AlignerPool& that);
	virtual ~AlignerPool();

	std::string outAlnDir_ = "";
	std::string inAlnDir_ = "";

	// request a aligner from the pool (wait on queue)
	PooledAligner popAligner();

	void initAligners();

	void destoryAligners();

};

class AlignerPoolException: public njh::err::Exception {
public:
	AlignerPoolException(const std::string s) :
			Exception(s) {

	}
	AlignerPoolException() :
			Exception("General AlignerPool exception.") {
	}
};

} // namespace concurrent
} // namespace bamcuts
