#pragma once
/*
 * StrandBiasKmerCorrector.hpp
 *
 *  Created on: Oct 16, 2026
 *      Author: nick
 */
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2026 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"

namespace njhseq {

/**@brief Parameters for detecting and correcting strand-restricted sequencing artifacts
 * (Illumina sequence-specific errors) at the k-mer level during BAM extraction.
 *
 * A k-mer is corrected to its Hamming-1 partner when the k-mer is seen on essentially one
 * sequencing strand only while the partner is dominant and strand balanced, with depth on the
 * strand the candidate is missing from. Corrections are applied base-wise to reads before they
 * are written out, so no extra files are produced.
 */
struct StrandBiasKmerCorrectionPars {
	bool on_{true};                               //!< run the correction (turn off with --noStrandBiasKmerCorrection in callers)
	uint32_t kLen_{25};                           //!< k-mer length, must be <= 32
	uint32_t minTotal_{5};                        //!< a candidate k-mer needs at least this many reads
	double oneSidedFrac_{0.05};                   //!< candidate minor-strand count must be <= max(1, oneSidedFrac * total)
	uint32_t minPartnerMissingStrandDepth_{20};   //!< partner must have at least this many reads on the strand the candidate is missing from
	double ratioHard_{0.15};                      //!< (candidate VAF on missing strand)/(candidate VAF on present strand) must be below this
	double fisherPhred_{30};                      //!< two-sided Fisher exact on [partner fwd, cand fwd; partner rev, cand rev], Phred scaled, must exceed this
	double partnerMargin_{2};                     //!< best partner total must be at least this many fold the second best partner
	uint32_t minVotes_{2};                        //!< a base is edited when at least this many k-mers covering it agree on the edit and none disagree
	uint32_t maxMinorStrandReads_{1};             //!< block an edit if any k-mer covering the base has more than this many reads on its minor strand
	uint32_t maxPasses_{5};                       //!< tally/correct passes, each pass re-tallies on the corrected reads; stops early when a pass proposes no new candidate k-mers
	bool padRegionByK_{true};                     //!< tally reads from the region padded by k on both sides so edge positions keep partner depth

	Json::Value toJson() const;
};

/**@brief Holds the per-pass correction maps and applies them to reads/alignments.
 *
 * Sequences are always handled in genome plus orientation; strand is the original sequencing
 * strand (BAM flag 0x10 for mapped reads, the opposite of the mate's strand for unmapped mates).
 */
class StrandBiasKmerCorrector {
public:
	struct StrandCounts {
		uint32_t fwd_{0};
		uint32_t rev_{0};
		uint32_t total() const { return fwd_ + rev_; }
		uint32_t minor() const { return std::min(fwd_, rev_); }
	};
	using KmerTable = std::unordered_map<uint64_t, StrandCounts>;

	struct Edit {
		uint32_t offset_{0}; //!< offset within the k-mer of the base to change
		char to_{'N'};       //!< base to change it to (the partner's base)
	};

	struct Pass {
		std::unordered_map<uint64_t, Edit> corrections_; //!< candidate k-mer hash -> edit
		std::unordered_set<uint64_t> balanced_;          //!< k-mers with more than maxMinorStrandReads_ on their minor strand (evidence against an edit)
		uint32_t oneSidedKmers_{0};                      //!< number of one-sided k-mers considered
	};

	explicit StrandBiasKmerCorrector(const StrandBiasKmerCorrectionPars & pars);

	const StrandBiasKmerCorrectionPars pars_;
	std::vector<Pass> passes_;
	std::unordered_set<uint64_t> candidatesSeen_; //!< every candidate k-mer proposed by any pass so far; a later pass only proposes new ones, which is what makes the pass loop converge

	//orientation and eligibility of alignments
	/**@brief true for primary alignments that are mapped or are an unmapped read with a mapped mate*/
	static bool alnEligible(const BamTools::BamAlignment & aln);
	/**@brief whether QueryBases are in genome plus orientation; mapped reads always are, an
	 * unmapped mate only when its 0x10 flag differs from its mate's (the SAM convention stores an
	 * unmapped mate with the mate's strand flag while a proper pair's mates are on opposite strands) */
	static bool alnIsPlusOriented(const BamTools::BamAlignment & aln);
	/**@brief original sequencing strand, true for reverse; for an unmapped mate the opposite of its mate*/
	static bool alnSeqStrandIsReverse(const BamTools::BamAlignment & aln);

	//k-mer hashing (2 bits per base, k <= 32)
	static bool hashKmer(const std::string & seq, uint32_t start, uint32_t k, uint64_t & out);
	static std::string unhash(uint64_t h, uint32_t k);
	static double fisherTwoSidedPhred(uint32_t a, uint32_t b, uint32_t c, uint32_t d);

	//pass building
	void tallySeq(KmerTable & table, const std::string & plusSeq, bool reverseStrand) const;
	Pass findCorrections(const KmerTable & table) const;
	void addPass(Pass pass);
	bool hasCorrections() const;
	uint32_t numPasses() const;

	//application
	/**@brief apply all passes in order to a plus oriented sequence, returns number of bases edited*/
	uint32_t applyToSeq(std::string & plusSeq) const; //!< sweeps the pass maps until the read is stable, at most maxPasses_ sweeps
	/**@brief apply to an alignment's QueryBases taking care of orientation, returns number of bases edited*/
	uint32_t correct(BamTools::BamAlignment & aln) const;

private:
	uint32_t applyPass(const Pass & pass, std::string & plusSeq) const;
};

}  // namespace njhseq
