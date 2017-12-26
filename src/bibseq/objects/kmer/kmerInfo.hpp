#pragma once
//
//  kmer.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 11/30/12.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//
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
#include "bibseq/objects/kmer/kmer.hpp"

namespace bibseq {

/**@brief class to hold kmer info for a single read or sequence
 *
 */
class kmerInfo {

public:
	/**@brief empty constructor
	 *
	 */
	kmerInfo();
	/**@brief Construct with kmer info from seq of kLength
	 *
	 * @param seq The seq to index the kmers for
	 * @param kLength the kmer length
	 * @param setReverse whether to do the reverse complement as well
	 */
	kmerInfo(const std::string & seq, uint32_t kLength, bool setReverse);

	/**@brief Construct with kmer info from seq of kLength
	 *
	 * @param seq The seq to index the kmers for
	 * @param kLength the kmer length
	 * @param setReverse whether to do the reverse complement as well
	 * @param pos the position to set setting kmers from
	 * @param len the length from pos to set kmers from
	 *
	 */
	kmerInfo(const std::string & seq, uint32_t kLength,
			size_t pos, uint32_t len, bool setReverse);

	/**@brief holder kmer information, key is kmer str and value is kmer class to hold counts and position information
	 *
	 */
	std::unordered_map<std::string, kmer> kmers_;
	/**@brief same as above but to hold reverse complement kmer information
	 *
	 */
	std::unordered_map<std::string, kmer> kmersRevComp_;
	/**@brief the length of kmer being held in kmers_ and kmersRevComp_
	 *
	 */
	uint32_t kLen_;

	uint64_t seqLen_;

	size_t seqPos_ = 0;

	bool infoSet_ = false;

	void setKmers(const std::string & seq, uint32_t kLength, bool setReverse);

	void setKmersFromPortion(const std::string & seq, uint32_t kLength,
			size_t pos, uint32_t len, bool setReverse);



	/**@brief Compare kmers between two reads
	 *
	 * @param read The other read to compare to
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmers(const kmerInfo & info) const;
	/**@brief Compare kmers of this seq against the reverse complement kmers of the other read
	 *
	 * @param read The other read to compare to
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmersRevComp(const kmerInfo & info) const;
	/**@brief Compare a sub set of kmers to the full kmers of the read
	 *
	 * @param read The read to compare against
	 * @param startPos Starting position to compare against
	 * @param windowSize The size of the window to use
	 * @return A pair, first is the number of kmers shared and second is num of kmers shared divided by max possible
	 */
	std::pair<uint32_t, double> compareSubKmersToFull(const kmerInfo & info,
			uint32_t startPos, uint32_t windowSize) const;
	/**@brief compare kmers only if they were found between startPos and startPos + windowSize - kLen_
	 *
	 * @param read The other read to compare to
	 * @param startPos The starting position to count from
	 * @param windowSize The size of the window
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmers(const kmerInfo & info,
			uint32_t startPos, uint32_t windowSize) const;

	/**@brief compare kmers only if they were found between startPos and startPos + windowSize - kLen_
	 *
	 * @param read The other read to compare to
	 * @param startPos The starting position to count from this sequence
	 * @param otherStartPos The starting position to count from in the other sequence
	 * @param windowSize The size of the window
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmers(const kmerInfo & info,
			uint32_t startPos,uint32_t otherStartPos, uint32_t windowSize) const;

	/**@brief Get a sliding window kmer comparison
	 *
	 * @param read The other read to compare to
	 * @param windowSize The size of the window
	 * @param windowStepSize The size of the step to start the comparison at the next window
	 * @return A map of std::pairs where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers, key is position
	 */
	std::unordered_map<size_t, std::pair<uint32_t, double>> slideCompareKmers(
			const kmerInfo & info, uint32_t windowSize,
			uint32_t windowStepSize) const;
	/**@brief Get a sliding window of kmer comparison to the full kmers of another read
	 *
	 * @param read The read compare it's full kmers to
	 * @param windowSize The size of the window
	 * @param windowStepSize The step size of the window
	 * @return A map of the various kmers shared in the window
	 */
	std::unordered_map<size_t,std::pair<uint32_t, double>> slideCompareSubKmersToFull(
			const kmerInfo & info, uint32_t windowSize,
			uint32_t windowStepSize) const;

	/**@brief Get a sliding window of kmer comparisons to sub windows of other kmers
	 *
	 * @param read The read compare it's full kmers to
	 * @param windowSize The size of the window
	 * @param windowStepSize The step size of the window
	 * @return A map of maps of the various kmers shared in the window, key 1 is position in this info, key 2 is position of other info
	 */
	std::unordered_map<size_t,std::unordered_map<size_t,std::pair<uint32_t, double>>> slideCompareSubKmersToSubKmers(
			const kmerInfo & info, uint32_t windowSize,
			uint32_t windowStepSize) const;


	static uint32_t getMinimumNonRedundant(const std::string & seq);



};
}  // namespace bibseq


