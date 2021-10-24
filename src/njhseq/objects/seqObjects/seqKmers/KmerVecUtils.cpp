/*
 * KmerVecUtils.cpp
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

#include "KmerVecUtils.hpp"
#include "njhseq/IO/SeqIO/SeqInput.hpp"

namespace njhseq {





std::vector<kmerCluster> greedyKmerSimCluster(const SeqIOOptions & inReadsOpts,
		uint32_t kLength, double kmerSimCutOff, bool checkComplement,
		bool verbose) {

	SeqInput reader(inReadsOpts);
	reader.openIn();
	return greedyKmerSimCluster(reader.readAllReads<readObject>(), kLength,
			kmerSimCutOff, checkComplement, verbose);
}

}  // namespace njhseq
