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

#include <array>
#include <limits>
#include <map>
#include <set>

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
	//!@brief Scale that veto with the k-mer's own depth, blocking only when the minor strand holds
	//! more than max(maxMinorStrandReads_, maxMinorStrandFrac_ * total) reads; 0 keeps the plain
	//! absolute threshold.
	//!
	//! The veto exists to stop an edit where a covering k-mer is genuinely present on both strands, so
	//! what matters is the share of that k-mer's reads on its minor strand, not the raw count. Held
	//! absolute, it tightens as coverage rises: a k-mer that is 97% one sided carries about one minor
	//! strand read at 90x and about four at 650x, so beyond a few hundred fold the artifact k-mer
	//! vetoes its own correction and the correction quietly stops working.
	//!
	//! The fraction has to stay below the share a genuinely two-sided allele carries, which is usually
	//! 0.3 to 0.5, so that a real allele still vetoes an edit; above about 0.25 real alleles start
	//! being corrected away. A very small fraction achieves little, since the artifact k-mer is itself
	//! the low-depth one and a few percent of its own reads is close to the absolute threshold it
	//! replaces.
	double maxMinorStrandFrac_{0.20};
	uint32_t maxPasses_{5};                       //!< tally/correct passes, each pass re-tallies on the corrected reads; stops early when a pass proposes no new candidate k-mers
	bool padRegionByK_{true};                     //!< tally reads from the region padded by k on both sides so edge positions keep partner depth

	//!@brief Quality-informed relaxation. A sequence-specific error degrades the whole strand in its
	//! window, so the reference base at the artifact position is itself much lower quality on the
	//! artifact strand than on the other strand. That drop is independent evidence of an artifact, so
	//! where it is present the two guards that lose power at low depth -- partner depth and the Fisher
	//! test -- are relaxed. One-sidedness, the VAF ratio, partner dominance and the minor-strand
	//! evidence guard are never relaxed.
	bool useQualityRelax_{true};
	//!@brief Require the quality collapse rather than only relaxing with it. Useful where strand
	//! restriction alone is not sufficient evidence, such as multigene families, in which a genuine
	//! rare paralog can have one sided coverage. In single-copy regions the strand and depth evidence
	//! stands on its own, so this is off by default.
	bool requireQualityDrop_{false};

	//!@brief Also correct an adjacent two-base (doublet) artifact, where the same miscall lands on two
	//! neighbouring positions on the same strand.
	//!
	//! The single-base criterion cannot reach these. A k-mer covering such a pair differs from the
	//! true k-mer at two positions, so it has no Hamming-1 partner and is never proposed; only a
	//! window ending between the two positions sees one of them alone, which is too few windows to
	//! reach minVotes_. The single-base criterion is not inert there either, since it still removes
	//! the reads carrying only one of the two miscalls, which can sharpen the remainder into a
	//! consistent false haplotype.
	//!
	//! A doublet partner is considered only when no Hamming-1 partner can explain the k-mer, and the
	//! partner margin is taken against the whole Hamming-1 neighbourhood so a single-base explanation
	//! always wins. A doublet is then held to stricter evidence than a single: the depth and Fisher
	//! guards must pass on their own with no quality relaxation, and the quality collapse is required
	//! at both offsets.
	bool doublet_{true};
	//!@brief Require the same substitution at both offsets of a doublet.
	//!
	//! Off by default, because requiring it misses a real artifact class: adjacent pairs where one base
	//! changes one way and its neighbour another, such as T to C alongside G to A. What keeps a doublet
	//! correction safe is one-sidedness, partner dominance and the quality collapse required at both
	//! offsets, not the substitution matching, so requiring it costs sensitivity without adding
	//! protection. Turn it on to correct only pairs that change the same way.
	bool doubletRequireSameSubstitution_{false};

	//!@brief Frequency-only mode: an alternative criterion that ignores strand entirely and corrects a
	//! rare k-mer to a dominant one-off neighbour purely on a depth ratio. Intended for data where the
	//! strand axis is unusable -- a stranded RNA-seq library, where strand collapses onto mate identity
	//! -- and where the errors are in the template (reverse transcription, RNA) rather than the
	//! basecall, so they carry normal quality and the strand/quality criteria cannot see them.
	//!
	//! It is NOT selective: it cannot distinguish an error from a real allele at the same frequency, so
	//! it erases real variants below roughly 1/(minDepthRatio_ + 1) of their locus, which loses minor
	//! strains and rare paralogs below that fraction. Off by default, and deliberately exclusive with
	//! the strand criterion rather than combined with it.
	bool frequencyOnly_{false};
	double minDepthRatio_{10};                  //!< partner must be at least this many fold the candidate
	uint32_t freqMinCandidateTotal_{1};         //!< a candidate may be as rare as a single read
	uint32_t freqMinPartnerDepth_{10};          //!< but the partner needs real support before we trust the ratio
	//!@brief Absolute ceiling on a frequency-mode candidate, regardless of the ratio; 0 disables it.
	//! The depth ratio cannot separate a sequencing error from a real low-copy variant, and in a
	//! tandem-repeat gene the two are indistinguishable by it, since a genuine repeat-copy variant
	//! sits at a fraction of the gene's coverage and differs by one base in similar context. Error
	//! mass is overwhelmingly singletons and doubletons, so a low ceiling keeps the correction while
	//! leaving repeat-copy variation alone. Depth is in reads; see freqCollapseDuplicates_.
	uint32_t freqMaxCandidateTotal_{3};
	//!@brief Frequency mode: count a duplicate fragment once in k-mer depth.
	//! A fragment is a duplicate when the full sequences of all its mates equal an earlier fragment's.
	//! Mates of a fragment still count separately.
	bool freqCollapseDuplicates_{true};
	//!@brief Frequency mode: also correct a candidate over freqMaxCandidateTotal_, or with a depth ratio
	//! between freqQualityMinDepthRatio_ and minDepthRatio_, when its differing base looks miscalled:
	//! its mean quality is at least qualDropMin_ below the partner's (partner has qualMinReads_ calls),
	//! its share of calls >= freqQualityHighQual_ is at most freqQualityMaxHighQualShare_ times the
	//! partner's share, and it has at most freqQualityMaxConfidentCalls_ such calls. With doublet_, an
	//! adjacent two-base candidate is admitted the same way, with the test required at both offsets.
	bool freqQualityOverride_{true};
	double freqQualityMinDepthRatio_{3};        //!< ratio required instead of minDepthRatio_ when the quality drop is seen
	uint32_t freqQualityHighQual_{30};          //!< Phred at or above which a call is confident
	double freqQualityMaxHighQualShare_{0.5};   //!< candidate's confident share must be at most this fraction of the partner's
	uint32_t freqQualityMaxConfidentCalls_{3};  //!< maximum confident calls the candidate may have
	double qualDropMin_{8};                       //!< required Phred drop at the candidate's offset: other strand minus artifact strand, or in frequency mode partner base minus candidate base
	uint32_t qualMinReads_{5};                    //!< partner reads needed on each strand (in frequency mode, in total) before a drop is trusted
	uint32_t relaxedPartnerMissingStrandDepth_{8};//!< partner depth required instead of minPartnerMissingStrandDepth_ when a drop is seen
	double relaxedFisherPhred_{10};               //!< Fisher Phred required instead of fisherPhred_ when a drop is seen

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
		uint32_t offset_{0};  //!< offset within the k-mer of the base to change
		char to_{'N'};        //!< base to change it to (the partner's base)
		uint32_t offset2_{0}; //!< second offset, for an adjacent doublet
		char to2_{'N'};       //!< 'N' when this is a single-base edit
		bool isDoublet() const { return 'N' != to2_; }
	};

	/**@brief A proposed correction before the guards that quality can relax have been resolved.*/
	struct Candidate {
		uint64_t kmer_{0};
		uint64_t partner_{0};
		uint32_t offset_{0};
		char to_{'N'};
		bool missingIsFwd_{false};            //!< strand the candidate k-mer is (nearly) absent from
		uint32_t partnerDepthMissing_{0};
		double fisherPhred_{0};
		bool failedDepth_{false};             //!< relaxable
		bool failedFisher_{false};            //!< relaxable
		uint32_t offset2_{0};                 //!< second offset, for an adjacent doublet
		char to2_{'N'};                       //!< 'N' when this is a single-base candidate
		bool needsQuality_{false};            //!< frequency mode: accepted only if the quality test passes
		bool isDoublet() const { return 'N' != to2_; }
	};

	/**@brief One region read held in genome plus orientation, so passes can re-tally without re-reading the BAM.*/
	struct BufferedRead {
		std::string seq_;
		std::vector<uint8_t> quals_;
		bool reverseStrand_{false};
		uint64_t fragment_{std::numeric_limits<uint64_t>::max()}; //!< shared by the two mates of a pair; max for a read with no mate buffered
	};

	//!@brief (k-mer, offset) -> {n fwd, n rev, quality sum fwd, quality sum rev} of the base at that offset;
	//! keyed by the partner k-mer, and in frequency mode by the candidate k-mer as well; slots 4 and 5
	//! count the calls at or above freqQualityHighQual_ {fwd, rev}
	using QualTable = std::map<std::pair<uint64_t, uint32_t>, std::array<double, 6>>;

	struct Pass {
		std::unordered_map<uint64_t, Edit> corrections_; //!< candidate k-mer hash -> edit
		std::unordered_set<uint64_t> balanced_;          //!< k-mers whose minor strand exceeds the veto threshold (evidence against an edit)
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
	/**@brief tally buffered reads; in frequency mode with freqCollapseDuplicates_, duplicate fragments are skipped*/
	void tallyReads(KmerTable & table, const std::vector<BufferedRead> & reads) const;
	/**@brief every one-sided k-mer with a dominant partner, with the relaxable guards recorded rather than applied*/
	std::vector<Candidate> findCandidates(const KmerTable & table) const;
	/**@brief mean quality of the partner's base at each candidate's offset, per strand, from the buffered reads*/
	QualTable collectOffsetQuality(const std::vector<BufferedRead> & reads, const std::vector<Candidate> & candidates) const;
	/**@brief turn candidates into a pass, accepting those whose only failures are relaxable and that show a quality drop*/
	Pass buildPass(const std::vector<Candidate> & candidates, const KmerTable & table, const QualTable & quals, uint32_t & qualityRescued) const;
	Pass findCorrections(const KmerTable & table) const; //!< convenience: findCandidates + buildPass with no quality information
	/**@brief Frequency mode: remove corrections that would edit a base other windows show to be real.
	 *
	 * A correction accepted on depth alone is removed when, in any buffered read, another window covering
	 * the base it edits is above freqMaxCandidateTotal_ and is not itself a correction in this pass.
	 * Quality-backed corrections that edit the same base of the same read are grouped into one variant;
	 * the variant's calls are pooled once per read and, if the pooled calls fail the confident-call tests
	 * against the partners, all of the variant's corrections are removed.
	 * Returns the number of corrections removed.
	 */
	uint32_t vetoSupportedAlleles(Pass & pass, const std::vector<Candidate> & candidates, const KmerTable & table,
			const QualTable & quals, const std::vector<BufferedRead> & reads) const;
	void addPass(Pass pass);
	bool hasCorrections() const;
	uint32_t numPasses() const;
	uint32_t numQualityRescued() const { return qualityRescued_; }
	void addQualityRescued(uint32_t n) { qualityRescued_ += n; }

	//application
	/**@brief apply all passes in order to a plus oriented sequence, returns number of bases edited*/
	uint32_t applyToSeq(std::string & plusSeq) const; //!< sweeps the pass maps until the read is stable, at most maxPasses_ sweeps
	/**@brief apply to an alignment's QueryBases taking care of orientation, returns number of bases edited*/
	uint32_t correct(BamTools::BamAlignment & aln) const;

private:
	uint32_t qualityRescued_{0}; //!< corrections admitted only because of a quality drop
	uint32_t applyPass(const Pass & pass, std::string & plusSeq) const;
};

}  // namespace njhseq
