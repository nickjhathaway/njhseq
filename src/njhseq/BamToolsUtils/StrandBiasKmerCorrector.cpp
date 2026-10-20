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
#include <functional>

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
	ret["maxMinorStrandFrac_"] = njh::json::toJson(maxMinorStrandFrac_);
	ret["maxPasses_"] = njh::json::toJson(maxPasses_);
	ret["padRegionByK_"] = njh::json::toJson(padRegionByK_);
	ret["useQualityRelax_"] = njh::json::toJson(useQualityRelax_);
	ret["requireQualityDrop_"] = njh::json::toJson(requireQualityDrop_);
	ret["doublet_"] = njh::json::toJson(doublet_);
	ret["doubletRequireSameSubstitution_"] = njh::json::toJson(doubletRequireSameSubstitution_);
	ret["frequencyOnly_"] = njh::json::toJson(frequencyOnly_);
	ret["minDepthRatio_"] = njh::json::toJson(minDepthRatio_);
	ret["freqMinCandidateTotal_"] = njh::json::toJson(freqMinCandidateTotal_);
	ret["freqMinPartnerDepth_"] = njh::json::toJson(freqMinPartnerDepth_);
	ret["freqMaxCandidateTotal_"] = njh::json::toJson(freqMaxCandidateTotal_);
	ret["freqCollapseDuplicates_"] = njh::json::toJson(freqCollapseDuplicates_);
	ret["freqQualityOverride_"] = njh::json::toJson(freqQualityOverride_);
	ret["freqQualityMinDepthRatio_"] = njh::json::toJson(freqQualityMinDepthRatio_);
	ret["freqQualityHighQual_"] = njh::json::toJson(freqQualityHighQual_);
	ret["freqQualityMaxHighQualShare_"] = njh::json::toJson(freqQualityMaxHighQualShare_);
	ret["freqQualityMaxConfidentCalls_"] = njh::json::toJson(freqQualityMaxConfidentCalls_);
	ret["qualDropMin_"] = njh::json::toJson(qualDropMin_);
	ret["qualMinReads_"] = njh::json::toJson(qualMinReads_);
	ret["relaxedPartnerMissingStrandDepth_"] = njh::json::toJson(relaxedPartnerMissingStrandDepth_);
	ret["relaxedFisherPhred_"] = njh::json::toJson(relaxedFisherPhred_);
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

void StrandBiasKmerCorrector::tallyReads(KmerTable & table, const std::vector<BufferedRead> & reads) const {
	if (!(pars_.frequencyOnly_ && pars_.freqCollapseDuplicates_)) {
		for (const auto & rd : reads) {
			tallySeq(table, rd.seq_, rd.reverseStrand_);
		}
		return;
	}
	//group the mates of each pair, keeping buffer order so which duplicate is counted is deterministic
	std::vector<std::vector<uint32_t>> fragments;
	std::unordered_map<uint64_t, uint32_t> fragmentIndex;
	for (uint32_t i = 0; i < reads.size(); ++i) {
		if (std::numeric_limits<uint64_t>::max() == reads[i].fragment_) {
			fragments.emplace_back(std::vector<uint32_t>{i});
			continue;
		}
		auto search = fragmentIndex.find(reads[i].fragment_);
		if (search == fragmentIndex.end()) {
			fragmentIndex.emplace(reads[i].fragment_, fragments.size());
			fragments.emplace_back(std::vector<uint32_t>{i});
		} else {
			fragments[search->second].emplace_back(i);
		}
	}
	std::unordered_set<std::string> seenFragments;
	for (const auto & fragment : fragments) {
		std::string key;
		for (const auto idx : fragment) {
			key += reads[idx].seq_;
			key += reads[idx].reverseStrand_ ? "-" : "+";
		}
		if (!seenFragments.emplace(std::move(key)).second) {
			continue; //a duplicate of an earlier fragment
		}
		for (const auto idx : fragment) {
			tallySeq(table, reads[idx].seq_, reads[idx].reverseStrand_);
		}
	}
}

std::vector<StrandBiasKmerCorrector::Candidate> StrandBiasKmerCorrector::findCandidates(const KmerTable & table) const {
	std::vector<Candidate> ret;
	const uint32_t k = pars_.kLen_;
	for (const auto & entry : table) {
		const auto & cand = entry.second;
		const uint32_t tot = cand.total();
		bool overCap = false;
		if (pars_.frequencyOnly_) {
			if (tot < pars_.freqMinCandidateTotal_) {
				continue;
			}
			overCap = pars_.freqMaxCandidateTotal_ > 0 && tot > pars_.freqMaxCandidateTotal_;
			if (overCap && !pars_.freqQualityOverride_) {
				continue; //too well supported to call an error; a repeat-copy variant looks like this
			}
		} else {
			if (tot < pars_.minTotal_) {
				continue;
			}
			if (cand.minor() > std::max<double>(1, pars_.oneSidedFrac_ * tot)) {
				continue;
			}
		}
		if (njh::in(entry.first, candidatesSeen_)) {
			continue; //already proposed by an earlier pass; re-proposing it would never converge
		}
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
					secondTotal = std::max(secondTotal, best.total());
					best = search->second;
					bestHash = neighbor;
					bestOffset = i;
					bestBase = CODE_BASES[code];
				} else if (ntot > secondTotal) {
					secondTotal = ntot;
				}
			}
		}
		//A one-off explanation always wins, so a doublet partner is only looked for when no single
		//neighbour is dominant. The margin is then taken against the whole single-off neighbourhood,
		//not just the other doublets, so a plausible single explanation still blocks the doublet.
		const bool singleUsable = ('N' != bestBase) && best.total() > tot
				&& !(secondTotal > 0 && best.total() < pars_.partnerMargin_ * secondTotal);
		//in frequency mode a doublet is only reachable through the quality override (see freqQualityOverride_)
		if (!singleUsable && pars_.doublet_ && (!pars_.frequencyOnly_ || pars_.freqQualityOverride_)) {
			const uint32_t bestSingleTotal = ('N' == bestBase) ? 0 : best.total();
			uint64_t dHash = 0;
			StrandCounts dBest;
			uint32_t dOffset = 0;
			char dBase = 'N';
			char dBase2 = 'N';
			uint32_t dSecond = bestSingleTotal; //any single neighbour competes with the doublet
			for (uint32_t i = 0; i + 1 < k; ++i) {
				const uint32_t shiftA = 2 * (k - 1 - i);
				const uint32_t shiftB = 2 * (k - 2 - i);
				const uint64_t origA = (entry.first >> shiftA) & 3;
				const uint64_t origB = (entry.first >> shiftB) & 3;
				for (uint64_t codeA = 0; codeA < 4; ++codeA) {
					if (codeA == origA) {
						continue;
					}
					for (uint64_t codeB = 0; codeB < 4; ++codeB) {
						if (codeB == origB) {
							continue;
						}
						if (pars_.doubletRequireSameSubstitution_ && !(origA == origB && codeA == codeB)) {
							continue; //the artifact changes both bases the same way
						}
						uint64_t neighbor = entry.first;
						neighbor &= ~(static_cast<uint64_t>(3) << shiftA);
						neighbor &= ~(static_cast<uint64_t>(3) << shiftB);
						neighbor |= (codeA << shiftA) | (codeB << shiftB);
						auto search = table.find(neighbor);
						if (search == table.end()) {
							continue;
						}
						const uint32_t ntot = search->second.total();
						if (ntot > dBest.total()) {
							dSecond = std::max(dSecond, dBest.total());
							dBest = search->second;
							dHash = neighbor;
							dOffset = i;
							dBase = CODE_BASES[codeA];
							dBase2 = CODE_BASES[codeB];
						} else if (ntot > dSecond) {
							dSecond = ntot;
						}
					}
				}
			}
			const bool doubletUsable = 'N' != dBase && dBest.total() > tot
					&& !(dSecond > 0 && dBest.total() < pars_.partnerMargin_ * dSecond);
			if (doubletUsable && pars_.frequencyOnly_) {
				//no depth-only path: a frequency-mode doublet always needs the quality drop at both offsets
				if (dBest.total() >= pars_.freqMinPartnerDepth_
						&& dBest.total() >= pars_.freqQualityMinDepthRatio_ * tot) {
					Candidate dc;
					dc.kmer_ = entry.first;
					dc.partner_ = dHash;
					dc.offset_ = dOffset;
					dc.to_ = dBase;
					dc.offset2_ = dOffset + 1;
					dc.to2_ = dBase2;
					dc.missingIsFwd_ = missingIsFwd;
					dc.needsQuality_ = true;
					ret.emplace_back(dc);
				}
			} else if (doubletUsable) {
				const uint32_t partnerMissing = missingIsFwd ? dBest.fwd_ : dBest.rev_;
				const uint32_t partnerPresent = missingIsFwd ? dBest.rev_ : dBest.fwd_;
				const uint32_t candMissing = missingIsFwd ? cand.fwd_ : cand.rev_;
				const uint32_t candPresent = missingIsFwd ? cand.rev_ : cand.fwd_;
				const double vafMissing = (candMissing + partnerMissing) > 0 ? candMissing / static_cast<double>(candMissing + partnerMissing) : 0;
				const double vafPresent = (candPresent + partnerPresent) > 0 ? candPresent / static_cast<double>(candPresent + partnerPresent) : 0;
				const double dRatio = vafPresent > 0 ? vafMissing / vafPresent : 1.0;
				if (dRatio < pars_.ratioHard_) {
					Candidate dc;
					dc.kmer_ = entry.first;
					dc.partner_ = dHash;
					dc.offset_ = dOffset;
					dc.to_ = dBase;
					dc.offset2_ = dOffset + 1;
					dc.to2_ = dBase2;
					dc.missingIsFwd_ = missingIsFwd;
					dc.partnerDepthMissing_ = partnerMissing;
					dc.fisherPhred_ = fisherTwoSidedPhred(dBest.fwd_, cand.fwd_, dBest.rev_, cand.rev_);
					//a doublet is never admitted by relaxation, so these are hard failures for it
					dc.failedDepth_ = partnerMissing < pars_.minPartnerMissingStrandDepth_;
					dc.failedFisher_ = dc.fisherPhred_ <= pars_.fisherPhred_;
					ret.emplace_back(dc);
				}
			}
		}
		if ('N' == bestBase) {
			continue; //no partner
		}
		//guards that are never relaxed
		if (best.total() <= tot) {
			continue; //partner not dominant
		}
		if (secondTotal > 0 && best.total() < pars_.partnerMargin_ * secondTotal) {
			continue; //ambiguous partner: we cannot tell which neighbour to correct toward
		}
		if (pars_.frequencyOnly_) {
			//depth ratio alone; quality is consulted only when the cap or the ratio refuses the candidate
			if (best.total() < pars_.freqMinPartnerDepth_) {
				continue;
			}
			const bool needsQuality = overCap || best.total() < pars_.minDepthRatio_ * tot;
			if (needsQuality && !(pars_.freqQualityOverride_
					&& best.total() >= pars_.freqQualityMinDepthRatio_ * tot)) {
				continue;
			}
			Candidate fc;
			fc.needsQuality_ = needsQuality;
			fc.kmer_ = entry.first;
			fc.partner_ = bestHash;
			fc.offset_ = bestOffset;
			fc.to_ = bestBase;
			fc.missingIsFwd_ = missingIsFwd;
			fc.partnerDepthMissing_ = missingIsFwd ? best.fwd_ : best.rev_;
			fc.fisherPhred_ = 0;
			fc.failedDepth_ = false;
			fc.failedFisher_ = false;
			ret.emplace_back(fc);
			continue;
		}
		const uint32_t partnerMissing = missingIsFwd ? best.fwd_ : best.rev_;
		const uint32_t partnerPresent = missingIsFwd ? best.rev_ : best.fwd_;
		const uint32_t candMissing = missingIsFwd ? cand.fwd_ : cand.rev_;
		const uint32_t candPresent = missingIsFwd ? cand.rev_ : cand.fwd_;
		const double vafMissing = (candMissing + partnerMissing) > 0 ? candMissing / static_cast<double>(candMissing + partnerMissing) : 0;
		const double vafPresent = (candPresent + partnerPresent) > 0 ? candPresent / static_cast<double>(candPresent + partnerPresent) : 0;
		const double ratio = vafPresent > 0 ? vafMissing / vafPresent : 1.0;
		if (ratio >= pars_.ratioHard_) {
			continue;
		}
		Candidate c;
		c.kmer_ = entry.first;
		c.partner_ = bestHash;
		c.offset_ = bestOffset;
		c.to_ = bestBase;
		c.missingIsFwd_ = missingIsFwd;
		c.partnerDepthMissing_ = partnerMissing;
		c.fisherPhred_ = fisherTwoSidedPhred(best.fwd_, cand.fwd_, best.rev_, cand.rev_);
		c.failedDepth_ = partnerMissing < pars_.minPartnerMissingStrandDepth_;
		c.failedFisher_ = c.fisherPhred_ <= pars_.fisherPhred_;
		ret.emplace_back(c);
	}
	return ret;
}

StrandBiasKmerCorrector::QualTable StrandBiasKmerCorrector::collectOffsetQuality(
		const std::vector<BufferedRead> & reads, const std::vector<Candidate> & candidates) const {
	QualTable ret;
	if (candidates.empty()) {
		return ret;
	}
	//only the candidates whose acceptance depends on quality need measuring
	std::unordered_map<uint64_t, std::set<uint32_t>> wanted;
	for (const auto & c : candidates) {
		if (pars_.frequencyOnly_) {
			//frequency mode compares the candidate's own base with the partner's, so both k-mers are measured
			if (c.needsQuality_) {
				for (const auto kmer : {c.kmer_, c.partner_}) {
					wanted[kmer].emplace(c.offset_);
					if (c.isDoublet()) {
						wanted[kmer].emplace(c.offset2_);
					}
				}
			}
			continue;
		}
		if (c.isDoublet()) {
			//a doublet is only ever accepted with a collapse at both of its offsets, so both are needed
			wanted[c.partner_].emplace(c.offset_);
			wanted[c.partner_].emplace(c.offset2_);
			continue;
		}
		if (pars_.requireQualityDrop_ || c.failedDepth_ || c.failedFisher_) {
			wanted[c.partner_].emplace(c.offset_);
		}
	}
	if (wanted.empty()) {
		return ret;
	}
	const uint32_t k = pars_.kLen_;
	uint64_t h = 0;
	for (const auto & rd : reads) {
		if (rd.seq_.size() < k || rd.quals_.size() != rd.seq_.size()) {
			continue;
		}
		const uint32_t nWindows = rd.seq_.size() - k + 1;
		for (uint32_t pos = 0; pos < nWindows; ++pos) {
			if (!hashKmer(rd.seq_, pos, k, h)) {
				continue;
			}
			auto search = wanted.find(h);
			if (search == wanted.end()) {
				continue;
			}
			for (const auto offset : search->second) {
				auto & slot = ret[{h, offset}];
				const uint8_t q = rd.quals_[pos + offset];
				slot[rd.reverseStrand_ ? 1 : 0] += 1;
				slot[rd.reverseStrand_ ? 3 : 2] += q;
				if (q >= pars_.freqQualityHighQual_) {
					slot[rd.reverseStrand_ ? 5 : 4] += 1;
				}
			}
		}
	}
	return ret;
}

StrandBiasKmerCorrector::Pass StrandBiasKmerCorrector::buildPass(const std::vector<Candidate> & candidates,
		const KmerTable & table, const QualTable & quals, uint32_t & qualityRescued) const {
	Pass ret;
	if (!pars_.frequencyOnly_) {
		//evidence that a base is genuinely present on both strands; meaningless without a strand axis
		for (const auto & entry : table) {
			//the threshold scales with the k-mer's depth, so the veto means the same thing at 90x as
			//at 650x rather than tightening as coverage rises
			const double allowedMinor = std::max<double>(pars_.maxMinorStrandReads_,
					pars_.maxMinorStrandFrac_ * entry.second.total());
			if (entry.second.minor() > allowedMinor) {
				ret.balanced_.insert(entry.first);
			}
		}
	}
	ret.oneSidedKmers_ = candidates.size();
	//quality collapse at one offset of a partner k-mer, measured from the buffered reads
	auto hasCollapseAt = [this, &quals](const Candidate & c, uint32_t offset) {
		auto search = quals.find({c.partner_, offset});
		if (search == quals.end()) {
			return false;
		}
		const auto & q = search->second;
		const double nMissing = c.missingIsFwd_ ? q[0] : q[1];
		const double nPresent = c.missingIsFwd_ ? q[1] : q[0];
		if (nMissing < pars_.qualMinReads_ || nPresent < pars_.qualMinReads_) {
			return false;
		}
		//the artifact strand is the one the candidate is present on, so a positive drop means the
		//reference base is degraded exactly where the candidate appears
		const double drop = (c.missingIsFwd_ ? q[2] : q[3]) / nMissing - (c.missingIsFwd_ ? q[3] : q[2]) / nPresent;
		return drop >= pars_.qualDropMin_;
	};
	//frequency-mode quality test at one offset: the candidate's base against the partner's
	auto hasFreqDropAt = [this, &quals](const Candidate & c, uint32_t offset) {
		auto candSearch = quals.find({c.kmer_, offset});
		auto partnerSearch = quals.find({c.partner_, offset});
		if (candSearch == quals.end() || partnerSearch == quals.end()) {
			return false;
		}
		const auto & cq = candSearch->second;
		const auto & pq = partnerSearch->second;
		const double nCand = cq[0] + cq[1];
		const double nPartner = pq[0] + pq[1];
		if (nCand < 1 || nPartner < pars_.qualMinReads_) {
			return false;
		}
		if ((pq[2] + pq[3]) / nPartner - (cq[2] + cq[3]) / nCand < pars_.qualDropMin_) {
			return false;
		}
		if (cq[4] + cq[5] > pars_.freqQualityMaxConfidentCalls_) {
			return false;
		}
		const double candHighShare = (cq[4] + cq[5]) / nCand;
		const double partnerHighShare = (pq[4] + pq[5]) / nPartner;
		return candHighShare <= pars_.freqQualityMaxHighQualShare_ * partnerHighShare;
	};
	for (const auto & c : candidates) {
		if (pars_.frequencyOnly_) {
			//the depth ratio was already required in findCandidates; quality only where it was relaxed
			if (c.needsQuality_) {
				if (!hasFreqDropAt(c, c.offset_) || (c.isDoublet() && !hasFreqDropAt(c, c.offset2_))) {
					continue;
				}
				++qualityRescued;
			}
			ret.corrections_[c.kmer_] = Edit { c.offset_, c.to_, c.offset2_, c.to2_ };
			continue;
		}
		if (c.isDoublet()) {
			//stricter than a single on purpose: the depth and Fisher guards must pass on their own, with
			//no relaxation path, and the quality collapse is required at both offsets
			if (c.failedDepth_ || c.failedFisher_) {
				continue;
			}
			if (!hasCollapseAt(c, c.offset_) || !hasCollapseAt(c, c.offset2_)) {
				continue;
			}
			ret.corrections_[c.kmer_] = Edit { c.offset_, c.to_, c.offset2_, c.to2_ };
			continue;
		}
		bool accept = !c.failedDepth_ && !c.failedFisher_;
		if (pars_.requireQualityDrop_) {
			//quality is a precondition, not a relaxation: a candidate is corrected only when the
			//reference base at its offset is measurably worse on the strand carrying the candidate
			accept = false;
			auto search = quals.find({c.partner_, c.offset_});
			if (search != quals.end()) {
				const auto & q = search->second;
				const double nMissing = c.missingIsFwd_ ? q[0] : q[1];
				const double nPresent = c.missingIsFwd_ ? q[1] : q[0];
				if (nMissing >= pars_.qualMinReads_ && nPresent >= pars_.qualMinReads_) {
					const double drop = (c.missingIsFwd_ ? q[2] : q[3]) / nMissing - (c.missingIsFwd_ ? q[3] : q[2]) / nPresent;
					if (drop >= pars_.qualDropMin_ && !c.failedDepth_ && !c.failedFisher_) {
						accept = true;
					} else if (drop >= pars_.qualDropMin_ && pars_.useQualityRelax_
							&& c.partnerDepthMissing_ >= pars_.relaxedPartnerMissingStrandDepth_
							&& c.fisherPhred_ > pars_.relaxedFisherPhred_) {
						accept = true;
						++qualityRescued;
					}
				}
			}
			if (accept) {
				ret.corrections_[c.kmer_] = Edit { c.offset_, c.to_, 0, 'N' };
			}
			continue;
		}
		if (!accept && pars_.useQualityRelax_) {
			auto search = quals.find({c.partner_, c.offset_});
			if (search != quals.end()) {
				const auto & q = search->second;
				const double nMissing = c.missingIsFwd_ ? q[0] : q[1];
				const double nPresent = c.missingIsFwd_ ? q[1] : q[0];
				if (nMissing >= pars_.qualMinReads_ && nPresent >= pars_.qualMinReads_) {
					const double meanMissing = (c.missingIsFwd_ ? q[2] : q[3]) / nMissing;
					const double meanPresent = (c.missingIsFwd_ ? q[3] : q[2]) / nPresent;
					//the artifact strand is the one the candidate is present on, so a positive drop
					//means the reference base is degraded exactly where the candidate appears
					const double drop = meanMissing - meanPresent;
					if (drop >= pars_.qualDropMin_
							&& c.partnerDepthMissing_ >= pars_.relaxedPartnerMissingStrandDepth_
							&& c.fisherPhred_ > pars_.relaxedFisherPhred_) {
						accept = true;
						++qualityRescued;
					}
				}
			}
		}
		if (accept) {
			ret.corrections_[c.kmer_] = Edit { c.offset_, c.to_, 0, 'N' };
		}
	}
	return ret;
}

StrandBiasKmerCorrector::Pass StrandBiasKmerCorrector::findCorrections(const KmerTable & table) const {
	uint32_t rescued = 0;
	return buildPass(findCandidates(table), table, QualTable{}, rescued);
}

uint32_t StrandBiasKmerCorrector::vetoSupportedAlleles(Pass & pass, const std::vector<Candidate> & candidates,
		const KmerTable & table, const QualTable & quals, const std::vector<BufferedRead> & reads) const {
	if (!pars_.frequencyOnly_ || pass.corrections_.empty()) {
		return 0;
	}
	std::unordered_map<uint64_t, const Candidate *> qualityBacked;
	for (const auto & c : candidates) {
		if (c.needsQuality_ && njh::in(c.kmer_, pass.corrections_)) {
			qualityBacked.emplace(c.kmer_, &c);
		}
	}
	const uint32_t k = pars_.kLen_;
	auto supported = [this, &table, &pass](uint64_t h) {
		if (njh::in(h, pass.corrections_)) {
			return false;
		}
		auto search = table.find(h);
		return search != table.end() && search->second.total() > pars_.freqMaxCandidateTotal_;
	};
	//quality-backed corrections that edit the same base of the same read are one variant
	std::unordered_map<uint64_t, uint64_t> parent;
	std::function<uint64_t(uint64_t)> findRoot = [&parent, &findRoot](uint64_t h) -> uint64_t {
		auto search = parent.find(h);
		if (search == parent.end() || search->second == h) {
			return h;
		}
		const uint64_t root = findRoot(search->second);
		parent[h] = root;
		return root;
	};
	for (const auto & qb : qualityBacked) {
		parent.emplace(qb.first, qb.first);
	}
	//per read: the bases quality-backed corrections edit, and which corrections edit them
	struct ReadEdits {
		uint32_t readIdx;
		std::vector<std::pair<uint32_t, uint64_t>> baseKmers;
	};
	std::vector<ReadEdits> readEdits;
	std::unordered_set<uint64_t> vetoed;
	std::vector<uint64_t> hashes;
	std::vector<bool> valid;
	for (uint32_t readIdx = 0; readIdx < reads.size(); ++readIdx) {
		const auto & rd = reads[readIdx];
		if (rd.seq_.size() < k) {
			continue;
		}
		const uint32_t nWindows = rd.seq_.size() - k + 1;
		hashes.assign(nWindows, 0);
		valid.assign(nWindows, false);
		bool hasCorrection = false;
		for (uint32_t pos = 0; pos < nWindows; ++pos) {
			uint64_t h = 0;
			valid[pos] = hashKmer(rd.seq_, pos, k, h);
			hashes[pos] = h;
			hasCorrection = hasCorrection || (valid[pos] && njh::in(h, pass.corrections_));
		}
		if (!hasCorrection) {
			continue;
		}
		ReadEdits edits{readIdx, {}};
		for (uint32_t pos = 0; pos < nWindows; ++pos) {
			if (!valid[pos]) {
				continue;
			}
			auto search = pass.corrections_.find(hashes[pos]);
			if (search == pass.corrections_.end()) {
				continue;
			}
			std::vector<uint32_t> editedBases{pos + search->second.offset_};
			if (search->second.isDoublet()) {
				editedBases.emplace_back(pos + search->second.offset2_);
			}
			if (njh::in(hashes[pos], qualityBacked)) {
				edits.baseKmers.emplace_back(editedBases.front(), hashes[pos]);
				continue;
			}
			//depth-only correction: removed if a covering window is above the cap and not itself a correction
			if (0 == pars_.freqMaxCandidateTotal_ || njh::in(hashes[pos], vetoed)) {
				continue;
			}
			bool veto = false;
			for (const auto basePos : editedBases) {
				const uint32_t firstWindow = basePos + 1 >= k ? basePos + 1 - k : 0;
				const uint32_t lastWindow = std::min(basePos, nWindows - 1);
				for (uint32_t w = firstWindow; w <= lastWindow && !veto; ++w) {
					veto = valid[w] && supported(hashes[w]);
				}
			}
			if (veto) {
				vetoed.emplace(hashes[pos]);
			}
		}
		if (!edits.baseKmers.empty()) {
			std::sort(edits.baseKmers.begin(), edits.baseKmers.end());
			for (uint32_t i = 1; i < edits.baseKmers.size(); ++i) {
				if (edits.baseKmers[i].first == edits.baseKmers[i - 1].first) {
					const uint64_t a = findRoot(edits.baseKmers[i].second);
					const uint64_t b = findRoot(edits.baseKmers[i - 1].second);
					if (a != b) {
						parent[a] = b;
					}
				}
			}
			readEdits.emplace_back(std::move(edits));
		}
	}
	//pool each variant's calls once per read and re-test them against the partners
	std::unordered_map<uint64_t, std::array<double, 2>> variantCalls; //root -> {calls, confident calls}
	std::unordered_set<uint64_t> countedThisRead;
	for (const auto & edits : readEdits) {
		const auto & rd = reads[edits.readIdx];
		countedThisRead.clear();
		for (const auto & baseKmer : edits.baseKmers) {
			const uint64_t root = findRoot(baseKmer.second);
			if (!countedThisRead.emplace(root).second || baseKmer.first >= rd.quals_.size()) {
				continue;
			}
			auto & calls = variantCalls[root];
			calls[0] += 1;
			if (rd.quals_[baseKmer.first] >= pars_.freqQualityHighQual_) {
				calls[1] += 1;
			}
		}
	}
	std::unordered_map<uint64_t, std::array<double, 2>> partnerCalls; //root -> {calls, confident calls}
	for (const auto & qb : qualityBacked) {
		auto search = quals.find({qb.second->partner_, qb.second->offset_});
		if (search == quals.end()) {
			continue;
		}
		auto & calls = partnerCalls[findRoot(qb.first)];
		calls[0] += search->second[0] + search->second[1];
		calls[1] += search->second[4] + search->second[5];
	}
	std::unordered_set<uint64_t> realVariants;
	for (const auto & variant : variantCalls) {
		auto partner = partnerCalls.find(variant.first);
		if (partner == partnerCalls.end() || partner->second[0] <= 0 || variant.second[0] <= 0) {
			continue;
		}
		const double candShare = variant.second[1] / variant.second[0];
		const double partnerShare = partner->second[1] / partner->second[0];
		if (variant.second[1] > pars_.freqQualityMaxConfidentCalls_
				|| candShare > pars_.freqQualityMaxHighQualShare_ * partnerShare) {
			realVariants.emplace(variant.first);
		}
	}
	for (const auto & qb : qualityBacked) {
		if (njh::in(findRoot(qb.first), realVariants)) {
			vetoed.emplace(qb.first);
		}
	}
	for (const auto h : vetoed) {
		pass.corrections_.erase(h);
	}
	return vetoed.size();
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
			if (search->second.isDoublet()) {
				++votes[pos + search->second.offset2_][search->second.to2_];
			}
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
