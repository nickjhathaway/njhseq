/*
 * StrandBiasKmerCorrector.cpp
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
#include "StrandBiasKmerCorrector.hpp"
#include "njhseq/helpers/seqUtil.hpp"

#include <boost/math/distributions/hypergeometric.hpp>

namespace njhseq {

Json::Value StrandBiasKmerCorrectionPars::toJson() const {
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["on_"] = njh::json::toJson(on_);
	ret["kLen_"] = njh::json::toJson(kLen_);
	ret["minTotal_"] = njh::json::toJson(minTotal_);
	ret["oneSidedFrac_"] = njh::json::toJson(oneSidedFrac_);
	ret["minPartnerMissingStrandDepth_"] = njh::json::toJson(minPartnerMissingStrandDepth_);
	ret["ratioHard_"] = njh::json::toJson(ratioHard_);
	ret["fisherPhred_"] = njh::json::toJson(fisherPhred_);
	ret["partnerMargin_"] = njh::json::toJson(partnerMargin_);
	ret["minVotes_"] = njh::json::toJson(minVotes_);
	ret["maxMinorStrandReads_"] = njh::json::toJson(maxMinorStrandReads_);
	ret["maxPasses_"] = njh::json::toJson(maxPasses_);
	ret["padRegionByK_"] = njh::json::toJson(padRegionByK_);
	return ret;
}

StrandBiasKmerCorrector::StrandBiasKmerCorrector(const StrandBiasKmerCorrectionPars & pars) :
		pars_(pars) {
	if (pars_.kLen_ < 5 || pars_.kLen_ > 32) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error kLen_ must be between 5 and 32, was " << pars_.kLen_ << "\n";
		throw std::runtime_error { ss.str() };
	}
}

bool StrandBiasKmerCorrector::alnEligible(const BamTools::BamAlignment & aln) {
	if (!aln.IsPrimaryAlignment()) {
		return false;
	}
	if (aln.IsMapped()) {
		return true;
	}
	return aln.IsPaired() && aln.IsMateMapped();
}

bool StrandBiasKmerCorrector::alnIsPlusOriented(const BamTools::BamAlignment & aln) {
	if (aln.IsMapped()) {
		return true;
	}
	//unmapped mate: stored with the mate's strand flag by convention, the pair's reads are on
	//opposite strands, so bases are plus oriented only when the flags differ
	return aln.IsReverseStrand() != aln.IsMateReverseStrand();
}

bool StrandBiasKmerCorrector::alnSeqStrandIsReverse(const BamTools::BamAlignment & aln) {
	if (aln.IsMapped()) {
		return aln.IsReverseStrand();
	}
	return !aln.IsMateReverseStrand();
}

namespace {
inline int baseCode(char c) {
	switch (c) {
	case 'A': case 'a': return 0;
	case 'C': case 'c': return 1;
	case 'G': case 'g': return 2;
	case 'T': case 't': return 3;
	default: return -1;
	}
}
const char CODE_BASES[4] = { 'A', 'C', 'G', 'T' };
}  // namespace

bool StrandBiasKmerCorrector::hashKmer(const std::string & seq, uint32_t start, uint32_t k, uint64_t & out) {
	uint64_t h = 0;
	for (uint32_t pos = start; pos < start + k; ++pos) {
		int code = baseCode(seq[pos]);
		if (code < 0) {
			return false;
		}
		h = (h << 2) | static_cast<uint64_t>(code);
	}
	out = h;
	return true;
}

std::string StrandBiasKmerCorrector::unhash(uint64_t h, uint32_t k) {
	std::string ret(k, 'N');
	for (uint32_t pos = 0; pos < k; ++pos) {
		ret[k - 1 - pos] = CODE_BASES[h & 3];
		h >>= 2;
	}
	return ret;
}

double StrandBiasKmerCorrector::fisherTwoSidedPhred(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
	// table [[a, b], [c, d]]; hypergeometric on a with row total a+b, column total a+c, grand total N
	const uint64_t N = static_cast<uint64_t>(a) + b + c + d;
	if (0 == N) {
		return 0;
	}
	const uint64_t r = static_cast<uint64_t>(a) + b;
	const uint64_t n = static_cast<uint64_t>(a) + c;
	boost::math::hypergeometric_distribution<double> dist(r, n, N);
	const uint64_t lower = (r + n > N) ? (r + n - N) : 0;
	const uint64_t upper = std::min(r, n);
	const double pObs = boost::math::pdf(dist, a);
	double p = 0;
	for (uint64_t x = lower; x <= upper; ++x) {
		double px = boost::math::pdf(dist, x);
		if (px <= pObs * (1 + 1e-7)) {
			p += px;
		}
	}
	p = std::min(1.0, p);
	if (p <= 1e-300) {
		return 3000;
	}
	return -10 * std::log10(p);
}

void StrandBiasKmerCorrector::tallySeq(KmerTable & table, const std::string & plusSeq, bool reverseStrand) const {
	const uint32_t k = pars_.kLen_;
	if (plusSeq.size() < k) {
		return;
	}
	uint64_t h = 0;
	for (uint32_t pos = 0; pos + k <= plusSeq.size(); ++pos) {
		if (!hashKmer(plusSeq, pos, k, h)) {
			continue;
		}
		auto & counts = table[h];
		if (reverseStrand) {
			++counts.rev_;
		} else {
			++counts.fwd_;
		}
	}
}

StrandBiasKmerCorrector::Pass StrandBiasKmerCorrector::findCorrections(const KmerTable & table) const {
	Pass ret;
	const uint32_t k = pars_.kLen_;
	for (const auto & entry : table) {
		const auto & cand = entry.second;
		if (cand.minor() > pars_.maxMinorStrandReads_) {
			ret.balanced_.insert(entry.first);
		}
		const uint32_t tot = cand.total();
		if (tot < pars_.minTotal_) {
			continue;
		}
		if (cand.minor() > std::max<double>(1, pars_.oneSidedFrac_ * tot)) {
			continue;
		}
		if (njh::in(entry.first, candidatesSeen_)) {
			continue; //already proposed by an earlier pass; it was either applied (and is gone from the reads) or blocked at application, re-proposing it would never converge
		}
		++ret.oneSidedKmers_;
		const bool missingIsFwd = cand.fwd_ <= cand.rev_;
		//find the best Hamming-1 partner
		uint64_t bestHash = 0;
		StrandCounts best;
		uint32_t bestOffset = 0;
		char bestBase = 'N';
		uint32_t secondTotal = 0;
		for (uint32_t i = 0; i < k; ++i) {
			const uint32_t shift = 2 * (k - 1 - i);
			const uint64_t origCode = (entry.first >> shift) & 3;
			for (uint64_t code = 0; code < 4; ++code) {
				if (code == origCode) {
					continue;
				}
				const uint64_t neighbor = (entry.first & ~(static_cast<uint64_t>(3) << shift)) | (code << shift);
				auto search = table.find(neighbor);
				if (search == table.end()) {
					continue;
				}
				const uint32_t ntot = search->second.total();
				if (ntot > best.total()) {
					secondTotal = best.total();
					best = search->second;
					bestHash = neighbor;
					bestOffset = i;
					bestBase = CODE_BASES[code];
				} else if (ntot > secondTotal) {
					secondTotal = ntot;
				}
			}
		}
		if ('N' == bestBase) {
			continue; //no partner
		}
		if (best.total() <= tot) {
			continue; //partner not dominant
		}
		if (secondTotal > 0 && best.total() < pars_.partnerMargin_ * secondTotal) {
			continue; //ambiguous partner
		}
		const uint32_t partnerMissing = missingIsFwd ? best.fwd_ : best.rev_;
		const uint32_t partnerPresent = missingIsFwd ? best.rev_ : best.fwd_;
		const uint32_t candMissing = missingIsFwd ? cand.fwd_ : cand.rev_;
		const uint32_t candPresent = missingIsFwd ? cand.rev_ : cand.fwd_;
		if (partnerMissing < pars_.minPartnerMissingStrandDepth_) {
			continue; //cannot tell artifact from one sided coverage
		}
		const double vafMissing = (candMissing + partnerMissing) > 0 ? candMissing / static_cast<double>(candMissing + partnerMissing) : 0;
		const double vafPresent = (candPresent + partnerPresent) > 0 ? candPresent / static_cast<double>(candPresent + partnerPresent) : 0;
		const double ratio = vafPresent > 0 ? vafMissing / vafPresent : 1.0;
		if (ratio >= pars_.ratioHard_) {
			continue;
		}
		if (fisherTwoSidedPhred(best.fwd_, cand.fwd_, best.rev_, cand.rev_) <= pars_.fisherPhred_) {
			continue;
		}
		(void) bestHash;
		ret.corrections_[entry.first] = Edit { bestOffset, bestBase };
	}
	return ret;
}

void StrandBiasKmerCorrector::addPass(Pass pass) {
	for (const auto & correction : pass.corrections_) {
		candidatesSeen_.insert(correction.first);
	}
	passes_.emplace_back(std::move(pass));
}

bool StrandBiasKmerCorrector::hasCorrections() const {
	for (const auto & pass : passes_) {
		if (!pass.corrections_.empty()) {
			return true;
		}
	}
	return false;
}

uint32_t StrandBiasKmerCorrector::numPasses() const {
	return passes_.size();
}

uint32_t StrandBiasKmerCorrector::applyPass(const Pass & pass, std::string & plusSeq) const {
	const uint32_t k = pars_.kLen_;
	if (pass.corrections_.empty() || plusSeq.size() < k) {
		return 0;
	}
	//collect votes per base position; votes[pos] = map of base -> count
	std::unordered_map<uint32_t, std::unordered_map<char, uint32_t>> votes;
	uint64_t h = 0;
	const uint32_t nWindows = plusSeq.size() - k + 1;
	for (uint32_t pos = 0; pos < nWindows; ++pos) {
		if (!hashKmer(plusSeq, pos, k, h)) {
			continue;
		}
		auto search = pass.corrections_.find(h);
		if (search != pass.corrections_.end()) {
			++votes[pos + search->second.offset_][search->second.to_];
		}
	}
	if (votes.empty()) {
		return 0;
	}
	uint32_t edited = 0;
	std::string out = plusSeq;
	for (const auto & vote : votes) {
		if (vote.second.size() != 1) {
			continue; //conflicting edits
		}
		const auto & edit = *vote.second.begin();
		if (edit.second < pars_.minVotes_) {
			continue;
		}
		const uint32_t basePos = vote.first;
		//evidence against: any covering window that is strand balanced
		bool blocked = false;
		const uint32_t firstWindow = basePos + 1 >= k ? basePos + 1 - k : 0;
		const uint32_t lastWindow = std::min(basePos, nWindows - 1);
		for (uint32_t w = firstWindow; w <= lastWindow && !blocked; ++w) {
			if (hashKmer(plusSeq, w, k, h) && njh::in(h, pass.balanced_)) {
				blocked = true;
			}
		}
		if (blocked) {
			continue;
		}
		out[basePos] = edit.first;
		++edited;
	}
	if (edited > 0) {
		plusSeq = out;
	}
	return edited;
}

uint32_t StrandBiasKmerCorrector::applyToSeq(std::string & plusSeq) const {
	//apply the pass maps repeatedly until the read stops changing: fixing one artifact base can turn a
	//neighbouring window from two mismatches into a one-off k-mer that an existing map already covers
	//(clustered artifacts), so a single sweep would leave those behind. Bounded by maxPasses_ sweeps.
	uint32_t edited = 0;
	for (uint32_t sweep = 0; sweep < std::max<uint32_t>(1, pars_.maxPasses_); ++sweep) {
		uint32_t editedThisSweep = 0;
		for (const auto & pass : passes_) {
			editedThisSweep += applyPass(pass, plusSeq);
		}
		edited += editedThisSweep;
		if (0 == editedThisSweep) {
			break;
		}
	}
	return edited;
}

uint32_t StrandBiasKmerCorrector::correct(BamTools::BamAlignment & aln) const {
	if (passes_.empty() || !alnEligible(aln) || aln.QueryBases.size() < pars_.kLen_) {
		return 0;
	}
	if (alnIsPlusOriented(aln)) {
		return applyToSeq(aln.QueryBases);
	}
	std::string plusSeq = seqUtil::reverseComplement(aln.QueryBases, "DNA");
	const uint32_t edited = applyToSeq(plusSeq);
	if (edited > 0) {
		aln.QueryBases = seqUtil::reverseComplement(plusSeq, "DNA");
	}
	return edited;
}

}  // namespace njhseq
