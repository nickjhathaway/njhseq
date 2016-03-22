#pragma once
/*
 * BamReaderPool.hpp
 *
 *  Created on: Jan 7, 2016
 *      Author: nick
 */

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

#pragma once

#include <vector>
#include <mutex>

#include "bibseq/alignment/aligner.h"

#include "bibseq/concurrency/ConcurrentQueue.hpp"

namespace bibseq {
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

class AlignerPoolException: public bib::err::Exception {
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
