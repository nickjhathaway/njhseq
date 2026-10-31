/*
 * KmerCollectionReadMatcher.cpp
 *
 *  Created on: Oct 30, 2026
 *      Author: nick
 */

#include "KmerCollectionReadMatcher.hpp"
#include "njhseq/objects/kmer/KmerGatherer.hpp"

namespace njhseq {

bool KmerCollectionReadMatcher::Match::hasWinner() const {
	return NO_WINNER != winnerIdx;
}

void KmerCollectionReadMatcher::Match::reset() {
	kmers.clear();
	kmersRevComp.clear();
	winnerIdx = NO_WINNER;
	winnerRevComp = false;
	bestFrac = 0;
}

namespace {
size_t totalKmers(const KmerCollectionReadMatcher::Collections & collections) {
	size_t ret = 0;
	for (const auto & collection : collections) {
		ret += collection.second.size();
	}
	return ret;
}
}  // namespace

KmerCollectionReadMatcher::KmerCollectionReadMatcher(const Collections & collections, uint32_t klen, bool useInvertedIndex, bool countUnique) :
		collections_(collections), klen_(klen), useInvertedIndex_(useInvertedIndex), countUnique_(countUnique),
		totalKmers_(totalKmers(collections)) {
	if (useInvertedIndex) {
		index_ = KmerSetInvertedIndex::get(collections_);
		if (nullptr != index_ && !index_->usable(collections_)) {
			index_ = nullptr;
		}
	}
	if (nullptr != index_) {
		names_ = index_->setNames_;
	} else {
		names_ = njh::getVecOfMapKeys(collections_);
		if (countUnique_) {
			uniqueCollections_ = KmerGatherer::filterToKmersUniqueToEachCollection(collections_);
		}
	}
}

const VecStr & KmerCollectionReadMatcher::names() const {
	return names_;
}

bool KmerCollectionReadMatcher::usable(const Collections & collections, uint32_t klen, bool useInvertedIndex, bool countUnique) const {
	return &collections == &collections_ && klen == klen_ && useInvertedIndex == useInvertedIndex_ && countUnique == countUnique_
			&& collections.size() == names_.size() && totalKmers(collections) == totalKmers_;
}

KmerCollectionReadMatcher::Match KmerCollectionReadMatcher::genMatch() const {
	Match ret;
	ret.found.resize(names_.size());
	ret.foundRevComp.resize(names_.size());
	ret.uniqueFound.resize(names_.size());
	ret.uniqueFoundRevComp.resize(names_.size());
	return ret;
}

void KmerCollectionReadMatcher::addKmers(const std::string & seq, bool revComp, Match & match) const {
	if (seq.size() <= klen_) {
		return;
	}
	for (uint32_t pos = 0; pos < seq.size() - klen_ + 1; ++pos) {
		match.kmers.emplace(hasher_.hash(seq, pos, klen_));
		if (revComp) {
			match.kmersRevComp.emplace(hasher_.revCompHash(seq, pos, klen_));
		}
	}
}

void KmerCollectionReadMatcher::countFound(Match & match) const {
	std::fill(match.found.begin(), match.found.end(), 0);
	std::fill(match.foundRevComp.begin(), match.foundRevComp.end(), 0);
	std::fill(match.uniqueFound.begin(), match.uniqueFound.end(), 0);
	std::fill(match.uniqueFoundRevComp.begin(), match.uniqueFoundRevComp.end(), 0);
	if (nullptr != index_) {
		for (const auto & kmer : match.kmers) {
			index_->tallyAllAndUniqueByIndex(kmer, match.found, match.uniqueFound);
		}
		for (const auto & kmer : match.kmersRevComp) {
			index_->tallyAllAndUniqueByIndex(kmer, match.foundRevComp, match.uniqueFoundRevComp);
		}
		return;
	}
	for (uint32_t idx = 0; idx < names_.size(); ++idx) {
		const auto & allKmers = collections_.at(names_[idx]);
		for (const auto & kmer : match.kmers) {
			match.found[idx] += njh::in(kmer, allKmers);
		}
		for (const auto & kmer : match.kmersRevComp) {
			match.foundRevComp[idx] += njh::in(kmer, allKmers);
		}
		if (countUnique_) {
			const auto & uniqueKmers = uniqueCollections_.at(names_[idx]);
			for (const auto & kmer : match.kmers) {
				match.uniqueFound[idx] += njh::in(kmer, uniqueKmers);
			}
			for (const auto & kmer : match.kmersRevComp) {
				match.uniqueFoundRevComp[idx] += njh::in(kmer, uniqueKmers);
			}
		}
	}
}

void KmerCollectionReadMatcher::determineWinner(Match & match) const {
	const auto nKmers = static_cast<double>(match.kmers.size());
	for (uint32_t idx = 0; idx < names_.size(); ++idx) {
		if (static_cast<double>(match.found[idx]) / nKmers > match.bestFrac) {
			match.bestFrac = static_cast<double>(match.found[idx]) / nKmers;
			match.winnerIdx = idx;
			match.winnerRevComp = false;
		}
		if (static_cast<double>(match.foundRevComp[idx]) / nKmers > match.bestFrac) {
			match.bestFrac = static_cast<double>(match.foundRevComp[idx]) / nKmers;
			match.winnerIdx = idx;
			match.winnerRevComp = true;
		}
	}
}

uint32_t KmerCollectionReadMatcher::countCollectionsHit(const Match & match, uint32_t minCount) {
	uint32_t ret = 0;
	for (uint32_t idx = 0; idx < match.uniqueFound.size(); ++idx) {
		if (match.uniqueFound[idx] >= minCount || match.uniqueFoundRevComp[idx] >= minCount) {
			++ret;
		}
	}
	return ret;
}

}  // namespace njhseq
