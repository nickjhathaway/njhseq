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

#include "consensusHelper.hpp"


namespace njhseq {


void consensusHelper::genConsensusFromCounters(seqInfo & info,
		const std::map<uint32_t, charCounter> & counters,
		const std::map<uint32_t, std::map<uint32_t, charCounter>> & insertions,
		const std::map<int32_t, charCounter> & beginningGap) {

	info.seq_.clear();
	info.qual_.clear();
	//for debugging;
	/*
	bool print = false;
	if(njh::containsSubString(info.name_, "lib1_Minor.00_seq.0001_5")){
		print = true;
	}*/
	// first deal with any gaps in the beginning
	double fortyPercent = 0.30 * info.cnt_;
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
	uint32_t countAdjustedToBeatForInserts = fortyPercent * 2;
	if(info.cnt_ < 10){
		countAdjustedToBeatForInserts = info.cnt_;
	}
	for (const auto & count : counters) {
		uint32_t bestQuality = 0;
		char bestBase = ' ';
		// if there is an insertion look at those if there is a majority of reads
		// with that insertion
//		std::cout << count.first << std::endl;
		auto search = insertions.find(count.first);
		if (search != insertions.end()) {
			for (auto & counterInsert : search->second) {
				bestQuality = 0;
				bestBase = ' ';
				char bestBaseTest = ' ';
				uint32_t bestQualTest = 0;
				counterInsert.second.getBest(bestBaseTest, bestQualTest);
//				std::cout << "\t"<< "bestBaseTest: " << bestBaseTest << ": " << counterInsert.second.chars_[bestBaseTest]<< " total:" << std::round(info.cnt_) << std::endl;
				counterInsert.second.getBest(bestBase, bestQuality, countAdjustedToBeatForInserts);

				if (bestBase == ' ') {
					continue;
				} else {
					info.seq_.push_back(bestBase);
					info.qual_.emplace_back(bestQuality);
				}
			}
		}
		count.second.getBest(bestBase, bestQuality);
		//for debugging;
		/*
		if(print && count.first == 1472){
			std::cout << __FILE__ << "  " << __LINE__ << "  " << __PRETTY_FUNCTION__ << std::endl;
			std::cout << info.name_ << std::endl;
			std::cout << "bestBase " << bestBase << std::endl;
			std::cout << "bestQuality " << bestQuality << std::endl;
			std::cout << "count.second.getTotalCount() " << count.second.getTotalCount() << std::endl;
			std::cout << "fortyPercent " << fortyPercent << std::endl;
		}*/
		if (bestBase == '-' || count.second.getTotalCount() < fortyPercent) {
			continue;
		}
//		std::cout << count.first << " bestBase: " << bestBase << " " << count.second.chars_[bestBase] << "/"<< count.second.getTotalCount()<< std::endl;
		info.seq_.push_back(bestBase);
		info.qual_.emplace_back(bestQuality);
	}

}


}  // namespace njhseq
