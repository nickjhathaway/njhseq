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

#include "consensusHelper.hpp"


namespace bibseq {


void consensusHelper::genConsensusFromCounters(seqInfo & info,
		const std::map<uint32_t, charCounter> & counters,
		const std::map<uint32_t, std::map<uint32_t, charCounter>> & insertions,
		const std::map<int32_t, charCounter> & beginningGap) {

	info.seq_.clear();
	info.qual_.clear();
	// first deal with any gaps in the beginning
	double fortyPercent = 0.40 * info.cnt_;
	for (const auto & bCount : beginningGap) {
		uint32_t bestQuality = 0;
		char bestBase = ' ';
		bCount.second.getBest(bestBase, bestQuality);
		if (bestBase == '-' || bCount.second.getTotalCount() < fortyPercent) {
			continue;
		}
		info.seq_.push_back(bestBase);
		info.qual_.emplace_back(bestQuality);
	}

	//read.seqBase_.outPutFastq(std::cout);
	// the iterators to over the letter counter maps
	for (const auto & count : counters) {

		uint32_t bestQuality = 0;
		char bestBase = ' ';
		// if there is an insertion look at those if there is a majority of reads
		// with that insertion

		auto search = insertions.find(count.first);
		if (search != insertions.end()) {

			for (auto & counterInsert : search->second) {
				bestQuality = 0;
				bestBase = ' ';
				counterInsert.second.getBest(bestBase, bestQuality,
						std::round(info.cnt_));
				if (bestBase == ' ') {
					continue;
				} else {
					info.seq_.push_back(bestBase);
					info.qual_.emplace_back(bestQuality);
				}
			}
		}

		count.second.getBest(bestBase, bestQuality);
		if (bestBase == '-' || count.second.getTotalCount() < fortyPercent) {
			continue;
		}

		info.seq_.push_back(bestBase);
		info.qual_.emplace_back(bestQuality);
	}

}


}  // namespace bibseq
