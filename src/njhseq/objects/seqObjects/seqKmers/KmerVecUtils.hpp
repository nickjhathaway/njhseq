#pragma once
/*
 * KmerVecUtils.hpp
 *
 *  Created on: May 24, 2016
 *      Author: nick
 */
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


#include "njhseq/objects/seqObjects/seqKmers/seqWithKmerInfo.hpp"
#include "njhseq/objects/seqObjects/seqKmers/KmerCluster.hpp"
#include "njhseq/IO/SeqIO/SeqIOOptions.hpp"


namespace njhseq {

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



template<typename READ>
std::vector<kmerCluster> greedyKmerSimCluster(const std::vector<READ> & inReads,
		uint32_t kLength, double kmerSimCutOff, bool checkComplement,
		bool verbose) {
	auto reads = createKmerReadVec(inReads, kLength, checkComplement);
	std::vector<kmerCluster> kClusters;
	for (const auto readPos : iter::range(reads.size())) {
		if (verbose) {
			std::cout << "currently on " << readPos << " of " << reads.size() << "\r";
		}
		bool foundMatch = false;
		for (auto & kClus : kClusters) {
			foundMatch = kClus.compareRead(reads[readPos], kmerSimCutOff,
					checkComplement);
			if (foundMatch) {
				break;
			}
		}
		if (!foundMatch) {
			kClusters.emplace_back(kmerCluster(reads[readPos]));
		}
	}
	if (verbose) {
		std::cout << std::endl;
	}
	return kClusters;
}

std::vector<kmerCluster> greedyKmerSimCluster(const SeqIOOptions & inReadsOpts,
		uint32_t kLength, double kmerSimCutOff, bool checkComplement,
		bool verbose);

}  // namespace njhseq
