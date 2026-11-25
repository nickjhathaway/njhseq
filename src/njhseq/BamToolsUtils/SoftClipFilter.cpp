/*
 * SoftClipFilter.cpp
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
#include "SoftClipFilter.hpp"
#include "njhseq/alignment/alignerUtils/alignerUtils.hpp"
#include "njhseq/alignment/aligner/alignCalc.hpp"

namespace njhseq {

namespace {

std::vector<uint8_t> qualsFromBam(const BamTools::BamAlignment & bAln, uint32_t start, uint32_t len) {
	std::vector<uint8_t> ret;
	ret.reserve(len);
	bool hasQuals = bAln.Qualities.size() == bAln.QueryBases.size();
	for (uint32_t pos = start; pos < start + len; ++pos) {
		ret.emplace_back(hasQuals ? static_cast<uint8_t>(bAln.Qualities[pos] - SangerQualOffset) : 40);
	}
	return ret;
}

std::string reversed(const std::string & seq) {
	return std::string(seq.rbegin(), seq.rend());
}

}  // namespace

seqInfo trimAlnSeqToRegion(const BamTools::BamAlignment & bAln,
		const GenomicRegion & region,
		const BamTools::RefVector & refData,
		bool keepInternalSoftClips) {
	std::function<bool(const BamTools::BamAlignment &, bool)> keepAll = [](const BamTools::BamAlignment &, bool) { return true; };
	return trimAlnSeqToRegion(bAln, region, refData, keepInternalSoftClips ? keepAll : std::function<bool(const BamTools::BamAlignment &, bool)>{});
}

seqInfo trimAlnSeqToRegion(const BamTools::BamAlignment & bAln,
		const GenomicRegion & region,
		const BamTools::RefVector & refData,
		const std::function<bool(const BamTools::BamAlignment &, bool)> & keepClip) {
	seqInfo querySeq = bamAlnToSeqInfo(bAln, true);
	GenomicRegion balnRegion(bAln, refData);
	if (region.end_ <= balnRegion.start_ || region.start_ >= balnRegion.end_) {
		//doesn't overlap the region, nothing sensible to trim to
		return querySeq;
	}
	seqInfo holderSeq(balnRegion.uid_, std::string(balnRegion.getLen(), 'N'));
	auto alnInfo = bamAlnToAlnInfoLocal(bAln);
	alignCalc::rearrangeLocal(holderSeq.seq_, querySeq.seq_, '-', alnInfo.begin()->second);
	alignCalc::rearrangeLocal(holderSeq.qual_, querySeq.qual_, 0, alnInfo.begin()->second);

	uint32_t startAln = 0;
	if (region.start_ > balnRegion.start_) {
		startAln = getAlnPosForRealPos(holderSeq.seq_, region.start_ - balnRegion.start_);
	}
	uint32_t endAln = len(holderSeq);
	if (region.end_ < balnRegion.end_) {
		endAln = getAlnPosForRealPos(holderSeq.seq_, region.end_ - balnRegion.start_ - 1) + 1;
	}
	auto ret = querySeq.getSubRead(startAln, endAln - startAln);
	ret.removeGaps();

	if (keepClip) {
		auto clipLens = SoftClipFilter::softClipLens(bAln);
		uint32_t keepLeading = 0;
		if (clipLens.first > 0 && region.start_ <= balnRegion.start_ && keepClip(bAln, true)) {
			keepLeading = std::min<uint32_t>(clipLens.first, balnRegion.start_ - region.start_);
		}
		uint32_t keepTrailing = 0;
		if (clipLens.second > 0 && region.end_ >= balnRegion.end_ && keepClip(bAln, false)) {
			keepTrailing = std::min<uint32_t>(clipLens.second, region.end_ - balnRegion.end_);
		}
		if (keepLeading > 0 || keepTrailing > 0) {
			const uint32_t qLen = bAln.QueryBases.size();
			std::string seq = bAln.QueryBases.substr(clipLens.first - keepLeading, keepLeading) + ret.seq_
					+ bAln.QueryBases.substr(qLen - clipLens.second, keepTrailing);
			auto quals = qualsFromBam(bAln, clipLens.first - keepLeading, keepLeading);
			quals.insert(quals.end(), ret.qual_.begin(), ret.qual_.end());
			auto trailingQuals = qualsFromBam(bAln, qLen - clipLens.second, keepTrailing);
			quals.insert(quals.end(), trailingQuals.begin(), trailingQuals.end());
			ret.seq_ = seq;
			ret.qual_ = quals;
		}
	}
	return ret;
}

Json::Value SoftClipTandemExtension::Pars::toJson() const {
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["maxPeriod_"] = njh::json::toJson(maxPeriod_);
	ret["maxMismatchFrac_"] = njh::json::toJson(maxMismatchFrac_);
	ret["minTract_"] = njh::json::toJson(minTract_);
	ret["freeRemainder_"] = njh::json::toJson(freeRemainder_);
	ret["flankIdentity_"] = njh::json::toJson(flankIdentity_);
	ret["refWindow_"] = njh::json::toJson(refWindow_);
	ret["minEntropy_"] = njh::json::toJson(minEntropy_);
	ret["minClipLen_"] = njh::json::toJson(minClipLen_);
	return ret;
}

uint32_t SoftClipTandemExtension::periodicTractLen(const std::string & clip,
		const std::string & anchor,
		uint32_t period,
		const Pars & pars) {
	const std::string s = clip + anchor;
	const uint32_t c = clip.size();
	//the anchor side of the junction has to carry the period for at least a unit
	const uint32_t nAnchor = std::max(period, pars.minTract_);
	if (c + nAnchor + period > s.size()) {
		return 0;
	}
	uint32_t anchorMismatches = 0;
	for (uint32_t i = 0; i < nAnchor; ++i) {
		if (s[c + i] != s[c + i + period]) {
			++anchorMismatches;
		}
	}
	if (anchorMismatches / static_cast<double>(nAnchor) > pars.maxMismatchFrac_) {
		return 0;
	}
	//walk away from the junction into the clip
	const uint32_t window = std::max<uint32_t>(2 * period, 8);
	std::deque<bool> hist;
	uint32_t windowMismatches = 0;
	uint32_t lastGood = 0;
	for (uint32_t k = 0; k < c; ++k) {
		const uint32_t i = c - 1 - k;
		const bool bad = s[i] != s[i + period];
		hist.push_back(bad);
		windowMismatches += bad;
		if (hist.size() > window) {
			windowMismatches -= hist.front();
			hist.pop_front();
		}
		if (k + 1 >= window && windowMismatches / static_cast<double>(hist.size()) > pars.maxMismatchFrac_) {
			break;
		}
		if (!bad) {
			lastGood = k + 1;
		}
	}
	return lastGood;
}

double SoftClipTandemExtension::bestUngappedIdentity(const std::string & query, const std::string & target) {
	if (query.empty() || target.size() < query.size()) {
		return 0;
	}
	uint32_t best = 0;
	for (uint32_t start = 0; start + query.size() <= target.size(); ++start) {
		uint32_t matches = 0;
		for (uint32_t pos = 0; pos < query.size(); ++pos) {
			if (query[pos] == target[start + pos]) {
				++matches;
			}
		}
		best = std::max(best, matches);
		if (best == query.size()) {
			break;
		}
	}
	return best / static_cast<double>(query.size());
}

double SoftClipTandemExtension::baseEntropy(const std::string & seq) {
	if (seq.empty()) {
		return 0;
	}
	std::unordered_map<char, uint32_t> counts;
	for (const auto base : seq) {
		++counts[base];
	}
	double ret = 0;
	for (const auto & count : counts) {
		const double frac = count.second / static_cast<double>(seq.size());
		ret -= frac * std::log2(frac);
	}
	return ret;
}

SoftClipTandemExtension::ClipResult SoftClipTandemExtension::classifyClip(const std::string & clip,
		const std::string & anchor,
		const std::string & refWindow,
		const Pars & pars) {
	ClipResult ret;
	for (uint32_t period = 1; period <= pars.maxPeriod_; ++period) {
		const auto tract = periodicTractLen(clip, anchor, period, pars);
		if (tract > ret.tractLen_) {
			ret.tractLen_ = tract;
			ret.period_ = period;
		}
	}
	if (ret.tractLen_ < pars.minTract_) {
		return ret;
	}
	ret.entropy_ = baseEntropy(clip.substr(clip.size() - ret.tractLen_));
	if (ret.entropy_ < pars.minEntropy_) {
		return ret;
	}
	const std::string remainder = clip.substr(0, clip.size() - ret.tractLen_);
	if (remainder.size() <= pars.freeRemainder_) {
		ret.pass_ = true;
	} else if (bestUngappedIdentity(remainder, refWindow) >= pars.flankIdentity_) {
		ret.pass_ = true;
		ret.usedFlank_ = true;
	}
	return ret;
}

SoftClipFilterPars::Mode SoftClipFilterPars::modeFromStr(const std::string & mode) {
	if ("hardCutOff" == mode) {
		return Mode::hardCutOff;
	} else if ("saElsewhere" == mode) {
		return Mode::saElsewhere;
	} else if ("saElsewhereOrInverted" == mode) {
		return Mode::saElsewhereOrInverted;
	} else if ("tandem" == mode) {
		return Mode::tandem;
	}
	std::stringstream ss;
	ss << __PRETTY_FUNCTION__ << ", error unknown soft clip filter mode: " << mode
			<< ", options are: " << njh::conToStr(modeOptions(), ", ") << "\n";
	throw std::runtime_error{ss.str()};
}

std::string SoftClipFilterPars::modeToStr(Mode mode) {
	switch (mode) {
		case Mode::hardCutOff:
			return "hardCutOff";
		case Mode::saElsewhere:
			return "saElsewhere";
		case Mode::saElsewhereOrInverted:
			return "saElsewhereOrInverted";
		case Mode::tandem:
			return "tandem";
	}
	return "hardCutOff";
}

SoftClipFilterPars::TrimClips SoftClipFilterPars::trimClipsFromStr(const std::string & trimClips) {
	if ("drop" == trimClips) {
		return TrimClips::drop;
	} else if ("explained" == trimClips) {
		return TrimClips::explained;
	} else if ("all" == trimClips) {
		return TrimClips::all;
	}
	std::stringstream ss;
	ss << __PRETTY_FUNCTION__ << ", error unknown trim soft clips option: " << trimClips
			<< ", options are: " << njh::conToStr(trimClipsOptions(), ", ") << "\n";
	throw std::runtime_error{ss.str()};
}

std::string SoftClipFilterPars::trimClipsToStr(TrimClips trimClips) {
	switch (trimClips) {
		case TrimClips::drop:
			return "drop";
		case TrimClips::explained:
			return "explained";
		case TrimClips::all:
			return "all";
	}
	return "explained";
}

VecStr SoftClipFilterPars::trimClipsOptions() {
	return VecStr{"drop", "explained", "all"};
}

VecStr SoftClipFilterPars::modeOptions() {
	return VecStr{"hardCutOff", "saElsewhere", "saElsewhereOrInverted", "tandem"};
}

Json::Value SoftClipFilterPars::toJson() const {
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["mode_"] = modeToStr(mode_);
	ret["trimClips_"] = trimClipsToStr(trimClips_);
	ret["saWindow_"] = njh::json::toJson(saWindow_);
	ret["tandem_"] = tandem_.toJson();
	ret["genome2bit_"] = njh::json::toJson(genome2bit_);
	return ret;
}

SoftClipFilter::SoftClipFilter(const SoftClipFilterPars & pars,
		uint32_t hardCutOff,
		double percentageCutOff,
		const BamTools::RefVector & refData) :
		pars_(pars), hardCutOff_(hardCutOff), percentageCutOff_(percentageCutOff), refData_(refData) {
	const bool haveGenome = !pars_.genome2bit_.empty() && bfs::exists(pars_.genome2bit_);
	if (SoftClipFilterPars::Mode::tandem == pars_.mode_ && !haveGenome) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error soft clip filter mode tandem needs a genome 2bit file, got: "
				<< (pars_.genome2bit_.empty() ? std::string("none") : pars_.genome2bit_.string()) << "\n";
		throw std::runtime_error{ss.str()};
	}
	//without a genome no clip counts as explained, so trimming with TrimClips::explained drops them all
	if (haveGenome) {
		tReader_ = std::make_unique<TwoBit::TwoBitFile>(pars_.genome2bit_);
	}
}

std::pair<uint32_t, uint32_t> SoftClipFilter::softClipLens(const BamTools::BamAlignment & bAln) {
	std::pair<uint32_t, uint32_t> ret{0, 0};
	for (auto it = bAln.CigarData.begin(); it != bAln.CigarData.end() && ('S' == it->Type || 'H' == it->Type); ++it) {
		if ('S' == it->Type) {
			ret.first += it->Length;
		}
	}
	for (auto it = bAln.CigarData.rbegin(); it != bAln.CigarData.rend() && ('S' == it->Type || 'H' == it->Type); ++it) {
		if ('S' == it->Type) {
			ret.second += it->Length;
		}
	}
	return ret;
}

bool SoftClipFilter::passesHardCutOff(const BamTools::BamAlignment & bAln) const {
	const auto amount = getSoftClipAmount(bAln);
	return amount < hardCutOff_ && amount / static_cast<double>(bAln.QueryBases.size()) < percentageCutOff_;
}

bool SoftClipFilter::passesSA(const BamTools::BamAlignment & bAln, const GenomicRegion & region, bool inverseIsChimeric) const {
	std::string sa;
	if (!bAln.GetTag("SA", sa) || sa.empty()) {
		return true;
	}
	const uint32_t padStart = region.start_ > pars_.saWindow_ ? region.start_ - pars_.saWindow_ : 0;
	const uint32_t padEnd = region.end_ + pars_.saWindow_;
	for (const auto & entry : njh::tokenizeString(sa, ";")) {
		if (entry.empty()) {
			continue;
		}
		auto toks = njh::tokenizeString(entry, ",");
		if (toks.size() < 3) {
			continue;
		}
		if (toks[0] != region.chrom_) {
			return false;
		}
		//SA positions are 1-based
		const auto pos = njh::StrToNumConverter::stoToNum<uint32_t>(toks[1]);
		if (pos == 0 || pos - 1 < padStart || pos - 1 >= padEnd) {
			return false;
		}
		if (inverseIsChimeric && (toks[2] == "-") != bAln.IsReverseStrand()) {
			return false;
		}
	}
	return true;
}

std::string SoftClipFilter::getRef(const std::string & chrom, uint32_t start, uint32_t end) {
	std::string ret;
	(*tReader_)[chrom]->getSequence(ret, start, end);
	njh::strToUpper(ret);
	return ret;
}

bool SoftClipFilter::clipIsTandemExtension(const BamTools::BamAlignment & bAln, bool leading) {
	if (!tReader_ || !bAln.IsMapped()) {
		return false;
	}
	const auto clipLens = softClipLens(bAln);
	const uint32_t c = leading ? clipLens.first : clipLens.second;
	const auto & tPars = pars_.tandem_;
	if (c < tPars.minClipLen_) {
		return false;
	}
	const auto & seq = bAln.QueryBases;
	const auto & chrom = refData_[bAln.RefID].RefName;
	const uint32_t chromLen = refData_[bAln.RefID].RefLength;
	if (leading) {
		const uint32_t alnStart = bAln.Position;
		const uint32_t refStart = alnStart > c + tPars.refWindow_ ? alnStart - c - tPars.refWindow_ : 0;
		const uint32_t refEnd = std::min(chromLen, alnStart + tPars.refWindow_);
		return SoftClipTandemExtension::classifyClip(seq.substr(0, c), seq.substr(c), getRef(chrom, refStart, refEnd), tPars).pass_;
	}
	const uint32_t alnEnd = bAln.GetEndPosition();
	const uint32_t refStart = alnEnd > tPars.refWindow_ ? alnEnd - tPars.refWindow_ : 0;
	const uint32_t refEnd = std::min(chromLen, alnEnd + c + tPars.refWindow_);
	//reverse so the clip sits to the left of the junction, the reference window is reversed too (not complemented)
	const std::string rSeq = reversed(seq);
	return SoftClipTandemExtension::classifyClip(rSeq.substr(0, c), rSeq.substr(c), reversed(getRef(chrom, refStart, refEnd)), tPars).pass_;
}

bool SoftClipFilter::passesTandem(const BamTools::BamAlignment & bAln) {
	const auto clipLens = softClipLens(bAln);
	const auto & tPars = pars_.tandem_;
	bool checkedAny = false;
	if (clipLens.first >= tPars.minClipLen_) {
		checkedAny = true;
		if (!clipIsTandemExtension(bAln, true)) {
			return false;
		}
	}
	if (clipLens.second >= tPars.minClipLen_) {
		checkedAny = true;
		if (!clipIsTandemExtension(bAln, false)) {
			return false;
		}
	}
	return checkedAny;
}

bool SoftClipFilter::keepClipWhenTrimming(const BamTools::BamAlignment & bAln, bool leading) {
	switch (pars_.trimClips_) {
		case SoftClipFilterPars::TrimClips::drop:
			return false;
		case SoftClipFilterPars::TrimClips::all:
			return true;
		case SoftClipFilterPars::TrimClips::explained:
			return clipIsTandemExtension(bAln, leading);
	}
	return false;
}

bool SoftClipFilter::passes(const BamTools::BamAlignment & bAln, const GenomicRegion & region) {
	if (passesHardCutOff(bAln)) {
		return true;
	}
	++failedHardCutOff_;
	bool rescue = false;
	switch (pars_.mode_) {
		case SoftClipFilterPars::Mode::hardCutOff:
			rescue = false;
			break;
		case SoftClipFilterPars::Mode::saElsewhere:
			rescue = passesSA(bAln, region, false);
			break;
		case SoftClipFilterPars::Mode::saElsewhereOrInverted:
			rescue = passesSA(bAln, region, true);
			break;
		case SoftClipFilterPars::Mode::tandem:
			rescue = bAln.IsMapped() && passesTandem(bAln);
			break;
	}
	if (rescue) {
		++rescued_;
	}
	return rescue;
}

}  // namespace njhseq
