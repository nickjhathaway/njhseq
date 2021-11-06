#pragma once
//
//  kmer.hpp
//
//  Created by Nick Hathaway on 11/30/12.
//
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
#include "njhseq/objects/kmer/kmer.hpp"

namespace njhseq {

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
	uint32_t kLen_ = std::numeric_limits<uint32_t>::max();

	uint64_t seqLen_{0};

	size_t seqPos_ = 0;

	bool infoSet_ = false;

	void setKmers(const std::string & seq, uint32_t kLength, bool setReverse);
	void updateKmers(const std::string & seq, bool setReverse);

	void setKmersFromPortion(const std::string & seq, uint32_t kLength,
			size_t pos, uint32_t len, bool setReverse);


	struct DetailedKmerDist{
		uint32_t totalShared_{0};//!< total kmers shared between the two datasets
		uint32_t totalKmersIn1_{0};//!< total number of kmers in 1
		uint32_t totalKmersIn2_{0};//!< total number of kmers in 2

		uint32_t totalUniqShared_{0};//!< total number of unique kmers shared
		uint32_t totalUniqKmersIn1_{0};//!< total number of unique kmers found in 1
		uint32_t totalUniqKmersIn2_{0};//!< total number of unique kmers found in 2

		uint32_t totalUniqBetween_{0};//!< total number of unique kmers between both sets

		/**
		 * @fn double getDistTotalShared()const
		 * @brief take the total number shared divided by the sum of total kmers in both datasets
		 *
		 * @return a number between 0 and 1, 0 being no kmers shared and 1 being all kmers between the datasets are shared
		 */
		double getDistTotalShared() const;

		/**
		 * @fn double getDistTotalSharedAdjusted()const
		 * @brief take the total number shared divided by minimum number of kmers in either dataset
		 *
		 * @return a number between 0 and 1, 0 being no kmers shared and 1 being all kmers from the smaller sequence can be found within the larger one
		 */
		double getDistTotalSharedLenAdjusted() const;

		/**
		 * @fn double getDistUniqueShared()const
		 * @brief take the total number of unqiue kmers (regardless of their within seq counts) shared between the datasets and divide by the total unique between the datasets
		 *
		 * @return a number between 0 and 1, 0 being no kmers shared, 1 being all unique kmers are shared between the two datasets
		 */
		double getDistUniqueShared() const;

		/**
		 * @fn double getDistUniqueSharedLenAdjusted()const
		 * @brief take the total number of unqiue kmers (regardless of their within seq counts) shared between the datasets and divide by the total unique between the datasets
		 *
		 * @return a number between 0 and 1, 0 being no kmers shared, 1 being all unique kmers within the smaller unique dataset can be found within the large unqiue dataset
		 */
		double getDistUniqueSharedLenAdjusted() const;




	};


	/**@brief Compare kmers between two reads
	 *
	 * @param read The other read to compare to
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmers(const kmerInfo & info) const;

	/**
		 * @fn DetailedKmerDist compareKmersDetailed(const kmerInfo&)const
	 * @brief Compare kmers of this seq against the kmers of the other read
	 *
	 * @param info The other info to compare to
	 * @return a detailed structure of distances
	 */
	DetailedKmerDist compareKmersDetailed(const kmerInfo & info) const;

	/**@brief Compare kmers of this seq against the reverse complement kmers of the other read
	 *
	 * @param read The other read to compare to
	 * @return a std::pair where first is the number of kmers shared
	 * and second is the fraction of kmers shared of the maximum possible shared kmers
	 */
	std::pair<uint32_t, double> compareKmersRevComp(const kmerInfo & info) const;

	/**
		 * @fn DetailedKmerDist compareKmersRevCompDetailed(const kmerInfo&)const
	 * @brief Get the detailed breakdown of kmers being shared between the reverse complement kmers of other seq compared to this seq
	 *
	 * @param info The other info to compare it's reverse complement seqs to
	 * @return a detailed structure of distances
	 */
	DetailedKmerDist compareKmersRevCompDetailed(const kmerInfo & info) const;

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


	/**@biref Get the smallest k-mer length so there are no repeated kmers
	 *
	 * @param seq the sequence to search
	 * @return the length
	 */
	static uint32_t getMinimumNonRedundant(const std::string & seq);

	/**@brief compute kmer entropy
	 *
	 * @return value ranging from 0 (all one kmer) to 2 (even amount of all kmers possible for given length)
	 */
	double computeKmerEntropy() const;



};
}  // namespace njhseq


