/*
 * HmmerUtility.cpp
 *
 *  Created on: Jan 16, 2022
 *      Author: nick
 */


#include "HmmerUtility.hpp"


namespace njhseq {



seqInfo HmmerUtility::trimSeqToHmmDomainHitEnv(const seqInfo &seq,
		const HmmerDomainHitTab &dom,
		bool mark) {
	auto region = dom.genBed6_env();
	return trimSeqToHmmDomainHit(seq, region, dom, mark);
}

seqInfo HmmerUtility::trimSeqToHmmDomainHitAln(const seqInfo &seq,
		const HmmerDomainHitTab &dom,
		bool mark) {
	auto region = dom.genBed6_aln();
	return trimSeqToHmmDomainHit(seq, region, dom, mark);
}

seqInfo HmmerUtility::trimSeqToHmmDomainHit(const seqInfo &seq,
		const Bed6RecordCore &region,
		const HmmerDomainHitTab &dom, bool mark) {

	auto subSeq = seq.getSubRead(region.chromStart_, region.length());
	if (mark) {
		MetaDataInName meta;
		if (MetaDataInName::nameHasMetaData(subSeq.name_)) {
			meta = MetaDataInName(subSeq.name_);
		}
		meta.addMeta("modelName", dom.queryName_, true);
		meta.addMeta("hmmTo", dom.hmmTo_, true);
		meta.addMeta("hmmFrom", dom.zeroBasedHmmFrom(), true);
		meta.addMeta("hmmCovered", dom.modelCoverage(), true);
		meta.addMeta("modelAccuracy", dom.acc_, true);
		meta.addMeta("score", dom.domainScore_, true);
		meta.addMeta("evalue", dom.domain_i_evalue_, true);
		meta.addMeta("trimStart", region.chromStart_, true);
		meta.addMeta("trimEnd", region.chromEnd_, true);
		meta.addMeta("trimLen", region.length(), true);

		meta.resetMetaInName(subSeq.name_);
	}
	return subSeq;
}
seqInfo HmmerUtility::trimSeqToHmmSeqHitEnv(const seqInfo &seq,
		const HmmerSeqHitTab &dom, bool mark){
	auto region = dom.genBed6_env();
	return trimSeqToHmmSeqHit(seq, region, dom, mark);
}

seqInfo HmmerUtility::trimSeqToHmmSeqHitAlign(const seqInfo &seq,
		const HmmerSeqHitTab &dom, bool mark){
	auto region = dom.genBed6_aln();
	return trimSeqToHmmSeqHit(seq, region, dom, mark);
}

seqInfo HmmerUtility::trimSeqToHmmSeqHit(const seqInfo &seq,
		const Bed6RecordCore &region,
		const HmmerSeqHitTab &dom,
		bool mark){
	seqInfo subSeq = seq.getSubRead(region.chromStart_, region.length());
	if (mark) {
		MetaDataInName meta;
		if (MetaDataInName::nameHasMetaData(subSeq.name_)) {
			meta = MetaDataInName(subSeq.name_);
		}
		meta.addMeta("hmmFrom", dom.zeroBasedHmmFrom(), true);
		meta.addMeta("hmmTo", dom.hmmTo_, true);
		meta.addMeta("trimIsRevStrand", njh::boolToStr(region.reverseStrand()), true);
		meta.addMeta("trimStart", region.chromStart_, true);
		meta.addMeta("trimEnd", region.chromEnd_, true);
		meta.addMeta("model", dom.targetName_, true);
		meta.addMeta("ID", dom.targetDesc_, true);
		meta.resetMetaInName(subSeq.name_);
	}
	if(region.reverseStrand()){
		subSeq.reverseComplementRead(false, true);
	}
	return subSeq;
}



void HmmerUtility::sortHmmSeqHitsByEvalue(std::vector<HmmerSeqHitTab> & hits){
	njh::sort(hits,[](const HmmerSeqHitTab & hit1, const HmmerSeqHitTab & hit2){
		if(hit1.modelEvalue_ == hit2.modelEvalue_){
			if(hit1.modelScore_ == hit2.modelScore_){
				return hit1.hmmLen() > hit2.hmmLen();
			}else{
				return hit1.modelScore_ > hit2.modelScore_;
			}
		}else{
			return hit1.modelEvalue_ < hit2.modelEvalue_;
		}
	});
}


void HmmerUtility::sortHmmDomainHitsByEvalue(std::vector<HmmerDomainHitTab> & hits){
	njh::sort(hits,[](const HmmerDomainHitTab & hit1, const HmmerDomainHitTab & hit2){
		if(hit1.domain_i_evalue_ == hit2.domain_i_evalue_){
			if(hit1.domainScore_ == hit2.domainScore_){
				return hit1.hmmLen() > hit2.hmmLen();
			}else{
				return hit1.domainScore_ > hit2.domainScore_;
			}
		}else{
			return hit1.domain_i_evalue_ < hit2.domain_i_evalue_;
		}
	});
}

void HmmerUtility::sortHmmDomainHitsByAccuracy(std::vector<HmmerDomainHitTab> & hits){
	njh::sort(hits,[](const HmmerDomainHitTab & hit1, const HmmerDomainHitTab & hit2){
		if(hit1.acc_ == hit2.acc_){
			if(hit1.domainScore_ == hit2.domainScore_){
				return hit1.hmmLen() > hit2.hmmLen();
			}else{
				return hit1.domainScore_ > hit2.domainScore_;
			}
		}else{
			return hit1.acc_ < hit2.acc_;
		}
	});
}

std::vector<HmmerSeqHitTab> HmmerUtility::getNonoverlapingSeqHitsPostSort(
		const std::vector<HmmerSeqHitTab> &hits,
		const std::function<Bed6RecordCore(const HmmerSeqHitTab)> &getBedLoc,
		uint32_t minOverlap ) {
	std::vector<HmmerSeqHitTab> ret;

	for (const auto &hit : hits) {
		bool overlaps = false;
		for (const auto &otherHit : ret) {
			if (getBedLoc(hit).overlaps(getBedLoc(otherHit), minOverlap)) {
				overlaps = true;
				break;
			}
		}
		if (!overlaps) {
			ret.emplace_back(hit);
		}
	}
	return ret;
}

std::vector<HmmerSeqHitTab> HmmerUtility::getNonoverlapingSeqHitsPostSortEnv(
		const std::vector<HmmerSeqHitTab> &hits, uint32_t minOverlap) {
	auto getBedEnv = [](const HmmerSeqHitTab &hit) {
		return hit.genBed6_env();
	};
	return getNonoverlapingSeqHitsPostSort(hits, getBedEnv, minOverlap);
}

std::vector<HmmerSeqHitTab> HmmerUtility::getNonoverlapingSeqHitsPostSortAln(
		const std::vector<HmmerSeqHitTab> &hits, uint32_t minOverlap) {
	auto getBedEnv = [](const HmmerSeqHitTab &hit) {
		return hit.genBed6_aln();
	};
	return getNonoverlapingSeqHitsPostSort(hits, getBedEnv, minOverlap);
}

std::vector<HmmerDomainHitTab> HmmerUtility::getNonoverlapingDomainHitsPostSort(
		const std::vector<HmmerDomainHitTab> &hits,
		const std::function<Bed6RecordCore(const HmmerDomainHitTab)> &getBedLoc,
		uint32_t minOverlap ) {
	std::vector<HmmerDomainHitTab> ret;

	for (const auto &hit : hits) {
		bool overlaps = false;
		for (const auto &otherHit : ret) {
			if (getBedLoc(hit).overlaps(getBedLoc(otherHit), minOverlap)) {
				overlaps = true;
				break;
			}
		}
		if (!overlaps) {
			ret.emplace_back(hit);
		}
	}
	return ret;
}

std::vector<HmmerDomainHitTab> HmmerUtility::getNonoverlapingDomainHitsPostSortEnv(
		const std::vector<HmmerDomainHitTab> &hits, uint32_t minOverlap) {
	auto getBedEnv = [](const HmmerDomainHitTab &hit) {
		return hit.genBed6_env();
	};
	return getNonoverlapingDomainHitsPostSort(hits, getBedEnv, minOverlap);
}

std::vector<HmmerDomainHitTab> HmmerUtility::getNonoverlapingDomainHitsPostSortAln(
		const std::vector<HmmerDomainHitTab> &hits, uint32_t minOverlap) {
	auto getBedEnv = [](const HmmerDomainHitTab &hit) {
		return hit.genBed6_aln();
	};
	return getNonoverlapingDomainHitsPostSort(hits, getBedEnv, minOverlap);
}




}  // namespace njhseq
