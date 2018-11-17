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


#include "njhseq/objects/seqObjects/seqKmers/seqWithKmerInfo.hpp"
#include "njhseq/IO/SeqIO/SeqIOOptions.hpp"

namespace njhseq {

class kmerCluster {
public:
	//constructor
	kmerCluster(std::unique_ptr<seqWithKmerInfo> & firstRead);
	//members
	std::unique_ptr<seqWithKmerInfo> mainRead_;
	std::vector<std::unique_ptr<seqWithKmerInfo>> reads_;

	//functions
	bool compareRead(std::unique_ptr<seqWithKmerInfo> & read, double cutOff,
			bool checkComplement);

	void writeInfo(const SeqIOOptions & outOpts) const;
};

class kmerClusterPos {
public:
	//constructor
	kmerClusterPos(std::unique_ptr<seqWithKmerInfo> & firstRead,
			uint64_t firstReadPos);
	//members
	std::unique_ptr<seqWithKmerInfo> mainRead_;
	std::vector<uint64_t> readPositions_;

	//functions
	bool compareRead(std::unique_ptr<seqWithKmerInfo> & read, uint64_t readPos,
			double cutOff, bool checkComplement);
};

} /* namespace njh */

