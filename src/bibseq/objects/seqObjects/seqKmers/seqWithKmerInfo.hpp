#pragma once
/*
 * seqWithKmerInfo.hpp
 *
 *  Created on: Jul 29, 2014
 *      Author: nickhathaway
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
#include "njhseq/objects/seqObjects/BaseObjects/baseReadObject.hpp"
#include "njhseq/helpers/seqUtil.hpp"
#include "njhseq/objects/kmer/kmerInfo.hpp"

namespace njhseq {




/**@brief A simple seq object with a njhseq::kmerInfo memeber to keep track of kmers
 *
 */
class seqWithKmerInfo : public baseReadObject{
public:
	/**@b Basic construct with seq info object
	 *
	 * @param info seq info object (contains seq, qual, count, fraction, name)
	 */
	seqWithKmerInfo(const seqInfo & info );

	/**@b Basic construct with seq info object and set kmers
	 *
	 * @param info seq info object (contains seq, qual, count, fraction, name)
	 */

	seqWithKmerInfo(const seqInfo & info, uint32_t kLength, bool setReverse);

	kmerInfo kInfo_; /**< Kmer information holder */

	/**@b set kmer information
	 *
	 * @param kLength The length of the kmer substring
	 * @param setReverse Whether to also set information for the reverse complement of the seq as well
	 */
	void setKmers(uint32_t kLength, bool setReverse);
	/**@b Compare kmers between two reads
	 *
	 * @param read The other read to compare to
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmers(const seqWithKmerInfo & read) const;
	/**@b Compare kmers of this seq against the reverse complement kmers of the other read
	 *
	 * @param read The other read to compare to
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmersRevComp(const seqWithKmerInfo & read) const;
	/**@b Compare a sub set of kmers to the full kmers of the read
	 *
	 * @param read The read to compare against
	 * @param startPos Starting position to compare against
	 * @param windowSize The size of the window to use
	 * @return A pair, first is the number of kmers shared and second is num of kmers shared divided by max possible
	 */
	std::pair<uint32_t, double> compareSubKmersToFull(const seqWithKmerInfo & read,
			uint32_t startPos, uint32_t windowSize) const;
	/**@b compare kmers only if they were found between startPos and startPos + windowSize - kLen_
	 *
	 * @param read The other read to compare to
	 * @param startPos The starting position to count from
	 * @param windowSize The size of the window
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmers(const seqWithKmerInfo & read,
			uint32_t startPos, uint32_t windowSize) const;
	/**@b Get a sliding window kmer comparison
	 *
	 * @param read The other read to compare to
	 * @param windowSize The size of the window
	 * @param windowStepSize The size of the step to start the comparison at the next window
	 * @return A vector of std::pairs where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::unordered_map<size_t, std::pair<uint32_t, double>> slideCompareKmers(const seqWithKmerInfo & read,
			uint32_t windowSize, uint32_t windowStepSize) const ;
	/**@b Get a sliding window of kmer comparison to the full kmers
	 *
	 * @param read The read compare it's full kmers to
	 * @param windowSize The size of the window
	 * @param windowStepSize The step size of the window
	 * @return A vector of the various kmers shared in the window
	 */
	std::unordered_map<size_t, std::pair<uint32_t, double>> slideCompareSubKmersToFull(
			const seqWithKmerInfo & read, uint32_t windowSize,
			uint32_t windowStepSize) const;

};






} /* namespace njh */


