#pragma once
/*
 * SoftClipFilter.hpp
 *
 *  Created on: Nov 24, 2026
 *      Author: nicholashathaway
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
#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

namespace njhseq {

/**@brief Trim an alignment's sequence to a region, returned in the plus strand orientation of the reference
 *
 * Bases are cut where the alignment crosses the region's start or end. A soft clip can be kept as long as it would
 * still sit inside the region if it continued off the end of the alignment, so a leading clip keeps at most
 * (alignment start - region start) bases, the ones next to the alignment, and a trailing clip at most
 * (region end - alignment end). An aligner often shows an insertion near the end of a read, e.g. extra copies of a
 * tandem repeat, as a soft clip, and those bases are part of the region.
 *
 * @param keepClip decides per clip (leading is true for the clip at the start of the alignment) whether it's kept,
 * an empty function drops all soft clipped bases
 */
seqInfo trimAlnSeqToRegion(const BamTools::BamAlignment & bAln,
		const GenomicRegion & region,
		const BamTools::RefVector & refData,
		const std::function<bool(const BamTools::BamAlignment &, bool)> & keepClip);

/**@brief keep all soft clips that fit in the region when keepInternalSoftClips is true, otherwise drop them all
 */
seqInfo trimAlnSeqToRegion(const BamTools::BamAlignment & bAln,
		const GenomicRegion & region,
		const BamTools::RefVector & refData,
		bool keepInternalSoftClips);

/**@brief Is a soft clip a tandem repeat extension, i.e. more copies of a repeat that the alignment ends in
 *
 * Works on a clip oriented to sit to the left of the junction, a trailing clip has to be reversed (not reverse
 * complemented) along with its anchor and the reference window.
 *
 * For each period p up to maxPeriod_ the anchor side of the junction has to carry the period for at least a unit,
 * then the clip is walked away from the junction comparing each base to the base p further on. The walk stops when
 * the mismatch rate in the trailing window of 2p bases is over maxMismatchFrac_. The period with the longest tract
 * is kept. Whatever of the clip is left past the tract, e.g. the flank on the far side of the repeat when the aligner
 * anchored on the near side, has to be no more than freeRemainder_ bases or match the reference window ungapped at
 * flankIdentity_ or better.
 *
 * Low complexity sequence matches itself at almost any shift, so the tract also has to have a base entropy of at
 * least minEntropy_ (bits), otherwise a poly-T or (AT)n chimera reads as a repeat extension.
 */
class SoftClipTandemExtension {
public:
	struct Pars {
		uint32_t maxPeriod_{30};
		double maxMismatchFrac_{0.25};
		uint32_t minTract_{8}; /**< the fewest clip bases the periodic tract has to cover */
		uint32_t freeRemainder_{5}; /**< clip bases past the tract allowed without matching the reference */
		double flankIdentity_{0.85};
		uint32_t refWindow_{250}; /**< bases of reference past the clip length on each side of the junction to search the remainder in */
		double minEntropy_{1.4};
		uint32_t minClipLen_{10}; /**< clips shorter than this aren't checked, they'd pass a hard cut off anyways */
		Json::Value toJson() const;
	};

	struct ClipResult {
		bool pass_{false};
		uint32_t period_{0};
		uint32_t tractLen_{0};
		double entropy_{0};
		bool usedFlank_{false};
	};

	/**@brief clip sits to the left of anchor, junction is between them */
	static ClipResult classifyClip(const std::string & clip,
			const std::string & anchor,
			const std::string & refWindow,
			const Pars & pars);

	static uint32_t periodicTractLen(const std::string & clip, const std::string & anchor, uint32_t period, const Pars & pars);
	static double bestUngappedIdentity(const std::string & query, const std::string & target);
	static double baseEntropy(const std::string & seq);
};

/**@brief How reads with more soft clipping than the hard cut off are handled */
struct SoftClipFilterPars {
	enum class Mode {
		hardCutOff, /**< filter them all */
		saElsewhere, /**< keep them unless a supplementary alignment (SA tag) lands outside the region, padded by saWindow_ */
		saElsewhereOrInverted, /**< saElsewhere and also filter a supplementary alignment in the region on the opposite strand (a foldback) */
		tandem /**< keep them only when every clip is a tandem repeat extension, needs genome2bit_ */
	};
	Mode mode_{Mode::hardCutOff};
	/**@brief which soft clipped bases are kept when a read is trimmed to the region */
	enum class TrimClips {
		drop, /**< none */
		explained, /**< clips that are tandem repeat extensions (SoftClipTandemExtension), needs genome2bit_ */
		all /**< all that would sit inside the region */
	};
	TrimClips trimClips_{TrimClips::explained};
	uint32_t saWindow_{1000};
	SoftClipTandemExtension::Pars tandem_;
	bfs::path genome2bit_; /**< needed by Mode::tandem */

	static Mode modeFromStr(const std::string & mode);
	static TrimClips trimClipsFromStr(const std::string & trimClips);
	static std::string trimClipsToStr(TrimClips trimClips);
	static VecStr trimClipsOptions();
	static std::string modeToStr(Mode mode);
	static VecStr modeOptions();
	Json::Value toJson() const;
};

/**@brief Decide whether an alignment passes soft clip filtering, rescuing some of what fails the hard cut off by SoftClipFilterPars::mode_ */
class SoftClipFilter {
public:
	SoftClipFilter(const SoftClipFilterPars & pars,
			uint32_t hardCutOff,
			double percentageCutOff,
			const BamTools::RefVector & refData);

	bool passes(const BamTools::BamAlignment & bAln, const GenomicRegion & region);

	/**@brief whether to keep a clip when trimming to the region, by SoftClipFilterPars::trimClips_
	 * @param leading the clip at the start of the alignment, otherwise the one at the end
	 */
	bool keepClipWhenTrimming(const BamTools::BamAlignment & bAln, bool leading);

	/**@brief is the leading or trailing clip a tandem repeat extension, false without a genome or when the clip is shorter than tandem_.minClipLen_
	 */
	bool clipIsTandemExtension(const BamTools::BamAlignment & bAln, bool leading);

	uint64_t failedHardCutOff_{0};
	uint64_t rescued_{0};

	/**@brief leading and trailing soft clip lengths, skipping over any hard clips */
	static std::pair<uint32_t, uint32_t> softClipLens(const BamTools::BamAlignment & bAln);

private:
	SoftClipFilterPars pars_;
	uint32_t hardCutOff_;
	double percentageCutOff_;
	const BamTools::RefVector & refData_;
	std::unique_ptr<TwoBit::TwoBitFile> tReader_;

	bool passesHardCutOff(const BamTools::BamAlignment & bAln) const;
	bool passesSA(const BamTools::BamAlignment & bAln, const GenomicRegion & region, bool inverseIsChimeric) const;
	bool passesTandem(const BamTools::BamAlignment & bAln);
	std::string getRef(const std::string & chrom, uint32_t start, uint32_t end);
};

}  // namespace njhseq
