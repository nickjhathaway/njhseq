//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
/*
 * seqWithKmerInfo.hpp
 *
 *  Created on: Jul 29, 2014
 *      Author: nickhathaway
 */

#include "bibseq/objects/seqObjects/baseReadObject.hpp"
#include "bibseq/helpers/seqUtil.hpp"
#include "bibseq/objects/helperObjects/kmer.hpp"

namespace bibseq {

class seqWithKmerInfo : public baseReadObject{
public:
	/**@b Basic construct with seq info object
	 *
	 * @param info seq info object (contains seq, qual, count, fraction, name)
	 */
	seqWithKmerInfo(const seqInfo & info ): baseReadObject(info){}

	/**@b Basic construct with seq info object and set kmers
	 *
	 * @param info seq info object (contains seq, qual, count, fraction, name)
	 */

	seqWithKmerInfo(const seqInfo & info, uint32_t kLength, bool setReverse):
		baseReadObject(info), kLen_(kLength){
		setKmers(kLength, setReverse);
	}

	/**@b holder kmer information, key is kmer str and value is kmer class to hold counts and position information
	 *
	 */
	std::unordered_map<std::string, kmer> kmers_;
	/**@b same as above but to hold reverse complement kmer information
	 *
	 */
	std::unordered_map<std::string, kmer> kmersRevComp_;
	/**@b the length of kmer being held in kmers_ and kmersRevComp_
	 *
	 */
	uint32_t kLen_ = 1;
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
	std::vector<std::pair<uint32_t, double>> slideCompareKmers(const seqWithKmerInfo & read,
			uint32_t windowSize, uint32_t windowStepSize) const ;

	virtual ~seqWithKmerInfo(){}
};

void allSetKmers(std::vector<std::unique_ptr<seqWithKmerInfo>> & reads,
		uint32_t kLength, bool setReverse);

void allSetKmers(std::vector<seqWithKmerInfo> & reads,
		uint32_t kLength, bool setReverse);

template<typename T>
std::vector<std::unique_ptr<seqWithKmerInfo>> createKmerReadVec(const std::vector<T> & reads){
	std::vector<std::unique_ptr<seqWithKmerInfo>> ret;
	for(const auto & read : reads){
		ret.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
	}
	return ret;
}

template<typename T>
std::vector<std::unique_ptr<seqWithKmerInfo>> createKmerReadVec(const std::vector<T> & reads,
		uint32_t kLength, bool setReverse){
	std::vector<std::unique_ptr<seqWithKmerInfo>> ret;
	for(const auto & read : reads){
		ret.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_, kLength, setReverse));
	}
	return ret;
}



class kmerCluster {
public:
	//constructor
	kmerCluster(std::unique_ptr<seqWithKmerInfo> & firstRead):
		mainRead_(std::forward<std::unique_ptr<seqWithKmerInfo>>(firstRead)){

	}
	//members
	std::unique_ptr<seqWithKmerInfo> mainRead_;
	std::vector<std::unique_ptr<seqWithKmerInfo>> reads_;

	//functions
	bool compareRead(std::unique_ptr<seqWithKmerInfo> & read,
			double cutOff, bool checkComplement);

	void writeInfo(std::ofstream & out) const;
};


class kmerClusterPos {
public:
	//constructor
	kmerClusterPos(std::unique_ptr<seqWithKmerInfo> & firstRead, uint64_t firstReadPos):
		mainRead_(std::forward<std::unique_ptr<seqWithKmerInfo>>(firstRead)),
				readPositions_({firstReadPos}){

	}
	//members
	std::unique_ptr<seqWithKmerInfo> mainRead_;
	std::vector<uint64_t> readPositions_;

	//functions
	bool compareRead(std::unique_ptr<seqWithKmerInfo> & read, uint64_t readPos,
			double cutOff, bool checkComplement);
};


} /* namespace bib */


