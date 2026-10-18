#pragma once

/*
 * KmerSetInvertedIndex.hpp
 *
 *  Created on: Oct 17, 2026
 *      Author: nick
 */

#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "njhseq/common.h"

namespace njhseq {

/**@brief One k-mer -> set(s) index over named k-mer sets, replacing a probe of every set per k-mer.
 *
 * Classifying a read against named k-mer sets is normally done by taking each of the read's k-mers
 * and asking every set whether it holds it, which costs one hash-set lookup per set per k-mer. This
 * inverts that: a single lookup per read k-mer names the set, or the sets, that hold it. The work
 * drops by roughly the number of sets, and the counts produced are identical, so the saving grows
 * with how many sets there are.
 *
 * K-mers shared between sets are a first-class case rather than a fallback. A shared k-mer's entry
 * carries SHARED_FLAG and indexes a list of every owning set, so a caller can either increment all
 * of them (`tallyAll`, matching a per-set loop with no break) or take the first (`firstSet`,
 * matching a per-set loop that breaks on the first match). Which of those a caller wants depends on
 * the loop it is replacing, so this makes no choice for it.
 *
 * Sets are keyed by address: `get` caches the index built for a given map and `usable` re-validates
 * it cheaply, so a caller whose sets change between iterations gets a rebuilt index automatically.
 */
class KmerSetInvertedIndex {
public:
	using SetMap = std::unordered_map<std::string, std::unordered_set<uint64_t>>;

	//!@brief marks an entry whose value indexes sharedSets_ rather than naming one set directly
	static constexpr uint32_t SHARED_FLAG = 0x80000000u;
	//!@brief returned by firstSet when the k-mer is in no set
	static constexpr uint32_t NO_SET = 0xFFFFFFFFu;

	std::unordered_map<uint64_t, uint32_t> kmerToSet_;
	std::vector<std::vector<uint32_t>> sharedSets_;
	std::vector<std::string> setNames_;
	const void * builtFrom_{nullptr};
	size_t nSets_{0};
	size_t totalKmers_{0};
	size_t nSharedKmers_{0};   //!< how many k-mers appear in more than one set, for reporting

	/**@brief Cheap check that this index still describes the given sets.
	 *
	 * Catches the common ways sets change without walking them: a different map, a set added or
	 * removed, or k-mers added or removed. It cannot catch a change that leaves both the set count and
	 * the k-mer total untouched, so a caller that mutates sets in place that way must rebuild.
	 */
	bool usable(const SetMap & sets) const {
		if (builtFrom_ != &sets || nSets_ != sets.size() || setNames_.size() != sets.size()) {
			return false;
		}
		size_t total = 0;
		for (const auto & s : sets) {
			total += s.second.size();
		}
		return total == totalKmers_;
	}

	/**@brief Increment every set that holds this k-mer. Matches a per-set loop with no break. */
	template<typename COUNTMAP>
	void tallyAll(uint64_t hash, COUNTMAP & into) const {
		auto search = kmerToSet_.find(hash);
		if (search == kmerToSet_.end()) {
			return;
		}
		if (search->second & SHARED_FLAG) {
			for (const auto idx : sharedSets_[search->second & ~SHARED_FLAG]) {
				++into[setNames_[idx]];
			}
		} else {
			++into[setNames_[search->second]];
		}
	}

	/**@brief The index of the first set holding this k-mer, or NO_SET. Matches a loop that breaks. */
	uint32_t firstSet(uint64_t hash) const {
		auto search = kmerToSet_.find(hash);
		if (search == kmerToSet_.end()) {
			return NO_SET;
		}
		if (search->second & SHARED_FLAG) {
			return sharedSets_[search->second & ~SHARED_FLAG].front();
		}
		return search->second;
	}

	/**@brief Build the index for these sets, or return the one already built for them.
	 *
	 * The hot path takes no lock: each thread remembers the index it last used and re-validates it,
	 * falling back to a mutex-guarded cache only on a miss.
	 */
	static std::shared_ptr<const KmerSetInvertedIndex> get(const SetMap & sets) {
		static thread_local const void * tlKey = nullptr;
		static thread_local std::shared_ptr<const KmerSetInvertedIndex> tlIndex;
		if (tlKey == &sets && nullptr != tlIndex && tlIndex->usable(sets)) {
			return tlIndex;
		}
		static std::mutex mut;
		static std::unordered_map<const void *, std::shared_ptr<const KmerSetInvertedIndex>> cache;
		{
			std::lock_guard<std::mutex> lock(mut);
			auto search = cache.find(&sets);
			if (search != cache.end() && search->second->usable(sets)) {
				tlKey = &sets;
				tlIndex = search->second;
				return tlIndex;
			}
		}
		auto built = std::make_shared<KmerSetInvertedIndex>();
		built->builtFrom_ = &sets;
		built->nSets_ = sets.size();
		size_t total = 0;
		for (const auto & s : sets) {
			total += s.second.size();
		}
		built->totalKmers_ = total;
		built->kmerToSet_.reserve(static_cast<size_t>(static_cast<double>(total) * 1.3));
		for (const auto & s : sets) {
			const uint32_t idx = static_cast<uint32_t>(built->setNames_.size());
			built->setNames_.emplace_back(s.first);
			for (const auto k : s.second) {
				auto res = built->kmerToSet_.emplace(k, idx);
				if (res.second) {
					continue;
				}
				//already present: this k-mer is shared between sets, so hold the full list
				uint32_t & slot = res.first->second;
				if (slot & SHARED_FLAG) {
					auto & lst = built->sharedSets_[slot & ~SHARED_FLAG];
					if (std::find(lst.begin(), lst.end(), idx) == lst.end()) {
						lst.emplace_back(idx);
					}
				} else {
					built->sharedSets_.emplace_back(std::vector<uint32_t>{slot, idx});
					slot = SHARED_FLAG | static_cast<uint32_t>(built->sharedSets_.size() - 1);
					++built->nSharedKmers_;
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock(mut);
			cache[&sets] = built;
		}
		tlKey = &sets;
		tlIndex = built;
		return tlIndex;
	}
};

}  // namespace njhseq
