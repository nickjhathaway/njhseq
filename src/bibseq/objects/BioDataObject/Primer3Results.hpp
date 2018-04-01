#pragma once

/*
 * Primer3Results.hpp
 *
 *  Created on: Aug 16, 2017
 *      Author: nick
 */

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


