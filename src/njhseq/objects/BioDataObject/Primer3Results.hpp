#pragma once

/*
 * Primer3Results.hpp
 *
 *  Created on: Aug 16, 2017
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
#include "njhseq/common.h"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"

namespace njhseq {
class Primer3Results {
public:
	struct region{
		uint32_t start_ = std::numeric_limits<uint32_t>::max();
		uint32_t size_ = std::numeric_limits<uint32_t>::max();

		Json::Value toJson() const;
	};


	class Primer {
	public:
		std::string name_ = "";
		double penalty_ = -1;
		std::string seq_ = ""; //seq in the file
		double tm_ = -1;
		double gc_percent_ = -1; // between 0 - 100
		double self_any_th_ = -1;
		double self_end_th_ = -1;
		double hairpin_th_ = -1;
		double end_stability_ = -1;
		bool right_ = false;

		VecStr problems_;

		region originalPos_; //zero based
		region forwardOrientationPos_; //zero based, for right primers the position given is the end position, this is the start in the forward orientation

		Json::Value toJson() const;
		GenomicRegion genRegion(const std::string & templateId) const;
	};


	class PrimerPair{
	public:
		std::string name_;
		Primer left_;
		Primer right_;
		double compl_any_th_{std::numeric_limits<double>::max()};
		double compl_end_th_{std::numeric_limits<double>::max()};
		double penalty_{std::numeric_limits<double>::max()};
		uint32_t product_size_{std::numeric_limits<uint32_t>::max()};

		Json::Value toJson() const;

		uint32_t getStart() const;
		uint32_t getEnd() const;

		bool overlaps(const PrimerPair & other, uint32_t minOverlap =1)const;
		uint32_t getOverlapLen(const PrimerPair & other)const;

	};



	std::string sequence_id_;
	std::string sequence_template_;
	std::string primer_task_{"generic"};
	std::vector<region> sequence_target_;
	std::vector<region> sequence_excluded_region_;
	uint32_t primer_left_num_returned_ = 0;
	uint32_t primer_right_num_returned_ = 0;
	uint32_t primer_internal_num_returned_ = 0;
	uint32_t primer_pair_num_returned_ = 0;
	VecStr warnings_;

	std::string primer_left_explain_;
	std::string primer_right_explain_;
	std::string primer_pair_explain_;

	std::unordered_multimap<std::string, std::string> misc_;//!< everything else that's not specifically handled

	std::vector<std::shared_ptr<PrimerPair>> primerPairs_;

	std::unordered_map<std::string, GenomicRegion> genPrimersRegions() const;

	std::shared_ptr<PrimerPair> getLowestPenaltyPair() const;// use only for "generic" mood
	std::vector<std::shared_ptr<PrimerPair>> getLowestPenaltyNonOverlappingPairs() const;// use only for "generic" mood

	Json::Value toJson() const;
	static std::vector<std::shared_ptr<Primer3Results>> parsePrimer3OutputResults(const bfs::path & input, bool ignoreUnexpected = false);
};


}  // namespace njhseq


