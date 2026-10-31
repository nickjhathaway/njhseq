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

/**@brief One k-mer -> set(s) index over named k-mer sets.
 *
 * A single lookup per k-mer gives the set, or the sets, that hold it, so the cost of looking up a
 * read's k-mers doesn't grow with the number of sets.
 *
 * A k-mer held by more than one set has an entry flagged with SHARED_FLAG that indexes the list of
 * every set holding it. `tallyAll` increments all of those sets and `firstSet` returns the first.
 *
 * Indexes are keyed by the address of the set map: `get` caches the index built for a map and
 * `usable` checks that it still matches, so an index is rebuilt when its sets change.
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

	/**@brief Whether this index still describes the given sets, checked by map address, set count and total k-mers.
	 *
	 * A change that keeps both the set count and the k-mer total isn't detected.
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

	/**@brief Increment every set that holds this k-mer. */
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

	/**@brief Increment every set that holds this k-mer in allInto, and also in uniqueInto when only one set holds it. */
	template<typename COUNTMAP>
	void tallyAllAndUnique(uint64_t hash, COUNTMAP & allInto, COUNTMAP & uniqueInto) const {
		auto search = kmerToSet_.find(hash);
		if (search == kmerToSet_.end()) {
			return;
		}
		if (search->second & SHARED_FLAG) {
			for (const auto idx : sharedSets_[search->second & ~SHARED_FLAG]) {
				++allInto[setNames_[idx]];
			}
		} else {
			++allInto[setNames_[search->second]];
			++uniqueInto[setNames_[search->second]];
		}
	}

	/**@brief tallyAllAndUnique with counts held by set index, the order of setNames_; both vectors must have nSets_ elements. */
	void tallyAllAndUniqueByIndex(uint64_t hash, std::vector<uint32_t> & allCounts, std::vector<uint32_t> & uniqueCounts) const {
		auto search = kmerToSet_.find(hash);
		if (search == kmerToSet_.end()) {
			return;
		}
		if (search->second & SHARED_FLAG) {
			for (const auto idx : sharedSets_[search->second & ~SHARED_FLAG]) {
				++allCounts[idx];
			}
		} else {
			++allCounts[search->second];
			++uniqueCounts[search->second];
		}
	}

	/**@brief The index of the first set holding this k-mer, or NO_SET. */
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
