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
#include "njhseq/programUtils/seqSetUp.hpp"

namespace njhseq {


class Primer3Runner {
public:

	struct Primer3Options{
		uint32_t PRIMER_MAX_SIZE = 35;
		uint32_t PRIMER_MIN_SIZE = 18;
		uint32_t PRIMER_OPT_SIZE = 33;
		uint32_t PRIMER_OPT_GC_PERCENT = 50;
		bfs::path primer3ConfigPath = "/usr/local/Cellar/primer3/2.4.0/share/primer3/primer3_config/";

		uint32_t minSize = 175;
		uint32_t maxSize = 290;
		uint32_t PRIMER_NUM_RETURN = 50;

		std::string task = "generic";


		//--PRIMER_MAX_SIZE,--PRIMER_MIN_SIZE,--PRIMER_OPT_SIZE
		void setPrimerSizeOpts(seqSetUp & setUp);

		//--minAmpSize,--maxAmpSize
		void setInsertSizeOptions(seqSetUp & setUp);

		//--primer3ConfigPath
		void setPrimaryOptions(seqSetUp & setUp);

		//--PRIMER_NUM_RETURN
		void setReturnOptions(seqSetUp & setUp);

	};



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

	class PrimerPair {
	public:

		PrimerPair(const Primer & left, const Primer & right): left_(left), right_(right){
			if(right.forwardOrientationPos_.start_ < left.forwardOrientationPos_.start_ + left.forwardOrientationPos_.size_){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "right position, "
						<<right_.forwardOrientationPos_.start_ << ":" << right_.forwardOrientationPos_.start_ + right_.forwardOrientationPos_.size_
						<< "if before the left primer, "
						<< left_.forwardOrientationPos_.start_ << ":" << left_.forwardOrientationPos_.start_ + left_.forwardOrientationPos_.size_<< "\n";
				throw std::runtime_error{ss.str()};
			}
			name_ = njh::pasteAsStr(left_.name_, "--", right_.name_);
		}
		PrimerPair(){

		}
		std::string name_;
		Primer left_;
		Primer right_;

		virtual Json::Value toJson() const;

		uint32_t getStart() const;
		uint32_t getEnd() const;

		bool overlaps(const PrimerPair &other, uint32_t minOverlap = 1) const;
		uint32_t getOverlapLen(const PrimerPair &other) const;




		virtual ~PrimerPair() {
		}
	};

	class PrimerPairGeneric: public PrimerPair {
	public:
		double compl_any_th_ { std::numeric_limits<double>::max() };
		double compl_end_th_ { std::numeric_limits<double>::max() };
		double penalty_ { std::numeric_limits<double>::max() };
		uint32_t product_size_ { std::numeric_limits<uint32_t>::max() };

		virtual Json::Value toJson() const;
		virtual ~PrimerPairGeneric() {
		}

	};

	class Primer3Results {
	public:



		std::string sequence_id_;
		std::string sequence_template_;
		std::string primer_task_;
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

		virtual ~Primer3Results() {
		}

		virtual Json::Value toJson() const;

	};

	class Primer3ResultsGeneric : public Primer3Results {
	public:

		virtual ~Primer3ResultsGeneric() {
		}

		virtual Json::Value toJson() const;
		std::vector<std::shared_ptr<Primer3Runner::PrimerPairGeneric>> primerPairs_;

		std::unordered_map<std::string, GenomicRegion> genPrimersRegions() const;

		std::shared_ptr<Primer3Runner::PrimerPairGeneric> getLowestPenaltyPair() const;// use only for "generic" mood
		std::vector<std::shared_ptr<Primer3Runner::PrimerPairGeneric>> getLowestPenaltyNonOverlappingPairs() const;// use only for "generic" mood


		static std::vector<std::shared_ptr<Primer3ResultsGeneric>> parsePrimer3OutputResults(const bfs::path & input, bool ignoreUnexpected = false);
	};


	class Primer3ResultsList : public Primer3Results {
	public:

		virtual ~Primer3ResultsList() {
		}

		virtual Json::Value toJson() const;
		std::vector<std::shared_ptr<Primer>> leftPrimers_;
		std::vector<std::shared_ptr<Primer>> rightPrimers_;

		std::unordered_map<std::string, GenomicRegion> genPrimersRegions() const;

		static std::vector<std::shared_ptr<Primer3ResultsList>> parsePrimer3OutputResults(const bfs::path & input, bool ignoreUnexpected = false);

	};





};





}  // namespace njhseq


