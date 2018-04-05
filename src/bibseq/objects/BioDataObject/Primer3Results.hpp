#pragma once

/*
 * Primer3Results.hpp
 *
 *  Created on: Aug 16, 2017
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bibseq/common.h"
#include "bibseq/objects/BioDataObject/GenomicRegion.hpp"

namespace bibseq {
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
		double penality_ = -1;
		std::string originalSeq_ = ""; //seq in the file
		std::string fowardOrientationSeq_ = ""; // for right primer this is the reverse complient of the seq in the file so it's in the forward direction
		double tm_ = -1;
		double gc_percent_ = -1; // between 0 - 100
		double self_any_th_ = -1;
		double self_end_th_ = -1;
		double hairpin_th_ = -1;
		double end_stability_ = -1;
		bool right_ = false;
		region originalPos_; //zero based
		region forwardOrientationPos_; //zero based, for right primers the position given is the end position, this is the start in the forward orientation
		Json::Value toJson() const;

		GenomicRegion genRegion(const std::string & templateId) const;
	};

	std::string sequence_id_;
	std::string sequence_template_;
	std::vector<region> sequence_target_;
	std::vector<region> sequence_excluded_region_;
	uint32_t primer_left_num_returned_ = 0;
	uint32_t primer_right_num_returned_ = 0;
	uint32_t primer_internal_num_returned_ = 0;
	uint32_t primer_pair_num_returned_ = 0;


	std::unordered_map<std::string, std::shared_ptr<Primer>> primers_;

	std::unordered_map<std::string, GenomicRegion> genPrimersRegions() const;

	Json::Value toJson() const;
	static std::vector<std::shared_ptr<Primer3Results>> parsePrimer3OutputResults(
			const bfs::path & input);
};


}  // namespace bibseq


