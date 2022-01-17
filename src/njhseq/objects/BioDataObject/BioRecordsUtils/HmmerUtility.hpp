#pragma once

/*
 * HmmerUtility.hpp
 *
 *  Created on: Jan 16, 2022
 *      Author: nick
 */

#include "njhseq/objects/BioDataObject/HmmerObjs.h"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"




namespace njhseq {

class HmmerUtility {
public:

	/**
		 * @fn seqInfo trimSeqToHmmDomainHitEnv(const seqInfo&, const HmmerDomainHitTab&, bool)
	 * @brief trim seq to the sub domain, envelope location
	 *
	 * @param seq the sequence to trim
	 * @param dom the hit
	 * @param mark whether or not to add trimming meta to output
	 * @return trimmed seq
	 */
	static seqInfo trimSeqToHmmDomainHitEnv(const seqInfo &seq,
			const HmmerDomainHitTab &dom, bool mark);

	/**
		 * @fn seqInfo trimSeqToHmmDomainHitAln(const seqInfo&, const HmmerDomainHitTab&, bool)
	 * @brief trim seq to the sub domain, align location
	 *
	 * @param seq the sequence to trim
	 * @param dom the hit
	 * @param mark whether or not to add trimming meta to output
	 * @return trimmed seq
	 */
	static seqInfo trimSeqToHmmDomainHitAln(const seqInfo &seq,
			const HmmerDomainHitTab &dom, bool mark);

	/**
		 * @fn seqInfo trimSeqToHmmDomainHit(const seqInfo&, const Bed6RecordCore&, const HmmerDomainHitTab&, bool)
	 * @brief should only be used by trimSeqToHmmDomainHitEnv or trimSeqToHmmDomainHitAln
	 *
	 * @param seq the sequence to trim
	 * @param region the region to trim to
	 * @param dom the hit
	 * @param mark whether or not to add trimming meta to output
	 * @return trimmed seq
	 */
	static seqInfo trimSeqToHmmDomainHit(const seqInfo &seq,
			const Bed6RecordCore &region,
			const HmmerDomainHitTab &dom, bool mark);


	/**
		 * @fn seqInfo trimSeqToHmmSeqHitEnv(const seqInfo&, const HmmerSeqHitTab&, bool)
	 * @brief assumes the input is what the hit comes from e.g. if hit is reverse complement, trimming will reverse complement seq
	 *
	 * @param seq the sequence to trim
	 * @param dom the hit
	 * @param mark whether or not to add trimming meta to output
	 * @return trimmed seq
	 */
	static seqInfo trimSeqToHmmSeqHitEnv(const seqInfo &seq,
			const HmmerSeqHitTab &dom, bool mark);
	/**
		 * @fn seqInfo trimSeqToHmmSeqHitAlign(const seqInfo&, const HmmerSeqHitTab&, bool)
	 * @brief assumes the input is what the hit comes from e.g. if hit is reverse complement, trimming will reverse complement seq
	 *
	 * @param seq the sequence to trim
	 * @param dom the hit
	 * @param mark whether or not to add trimming meta to output
	 * @return trimmed seq
	 */
	static seqInfo trimSeqToHmmSeqHitAlign(const seqInfo &seq,
			const HmmerSeqHitTab &dom, bool mark);

	/**
		 * @fn seqInfo trimSeqToHmmSeqHit(const seqInfo&, const Bed6RecordCore&, const HmmerSeqHitTab&, bool)
	 * @brief should only be used by trimSeqToHmmSeqHitEnv or trimSeqToHmmSeqHitAlign
	 *
	 * @param seq the sequence to trim
	 * @param region the region to extract from
	 * @param dom the hit
	 * @param mark whether or not to add trimming meta to output
	 * @return trimmed seq
	 */
	static seqInfo trimSeqToHmmSeqHit(const seqInfo &seq,
			const Bed6RecordCore &region,
			const HmmerSeqHitTab &dom,
			bool mark);


	static void sortHmmSeqHitsByEvalue(std::vector<HmmerSeqHitTab> & hits);
	static void sortHmmDomainHitsByEvalue(std::vector<HmmerDomainHitTab> & hits);
	static void sortHmmDomainHitsByAccuracy(std::vector<HmmerDomainHitTab> & hits);

	static std::vector<HmmerSeqHitTab> getNonoverlapingSeqHitsPostSort(
			const std::vector<HmmerSeqHitTab> &hits,
			const std::function<Bed6RecordCore(const HmmerSeqHitTab)> &getBedLoc,
			uint32_t minOverlap = 1);

	static std::vector<HmmerSeqHitTab> getNonoverlapingSeqHitsPostSortEnv(
			const std::vector<HmmerSeqHitTab> &hits, uint32_t minOverlap = 1);

	static std::vector<HmmerSeqHitTab> getNonoverlapingSeqHitsPostSortAln(
			const std::vector<HmmerSeqHitTab> &hits, uint32_t minOverlap = 1);



	static std::vector<HmmerDomainHitTab> getNonoverlapingDomainHitsPostSort(
			const std::vector<HmmerDomainHitTab> &hits,
			const std::function<Bed6RecordCore(const HmmerDomainHitTab)> &getBedLoc,
			uint32_t minOverlap = 1);

	static std::vector<HmmerDomainHitTab> getNonoverlapingDomainHitsPostSortEnv(
			const std::vector<HmmerDomainHitTab> &hits, uint32_t minOverlap = 1);

	static std::vector<HmmerDomainHitTab> getNonoverlapingDomainHitsPostSortAln(
			const std::vector<HmmerDomainHitTab> &hits, uint32_t minOverlap = 1);

};


}  // namespace njhseq


