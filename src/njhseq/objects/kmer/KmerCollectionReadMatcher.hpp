#pragma once

/*
 * KmerCollectionReadMatcher.hpp
 *
 *  Created on: Oct 30, 2026
 *      Author: nick
 */

#include "njhseq/common.h"
#include "njhseq/objects/kmer/SimpleKmerHash.hpp"
#include "njhseq/objects/kmer/KmerSetInvertedIndex.hpp"

namespace njhseq {

/**
 * @brief Matches reads against k-mer collections
 *
 * A read's hashed k-mers are counted against every collection, both all of a collection's k-mers and its unique k-mers
 * (found in no other collection), and the collection holding the highest fraction of the read's k-mers wins.
 * Counts are held by collection index (the order of names()) in a Match, and k-mers are looked up in a KmerSetInvertedIndex.
 */
class KmerCollectionReadMatcher {
public:
	using Collections = std::unordered_map<std::string, std::unordered_set<uint64_t>>;

	/**
	 * @brief The k-mers, counts and winner for one read, can be reused between reads after reset()
	 */
	struct Match {
		static constexpr uint32_t NO_WINNER = std::numeric_limits<uint32_t>::max();

		std::unordered_set<uint64_t> kmers;        ///< the read's distinct hashed k-mers
		std::unordered_set<uint64_t> kmersRevComp; ///< the reverse complement's distinct hashed k-mers
		std::vector<uint32_t> found;               ///< per collection, read k-mers in the collection
		std::vector<uint32_t> foundRevComp;        ///< per collection, reverse complement k-mers in the collection
		std::vector<uint32_t> uniqueFound;         ///< per collection, read k-mers unique to the collection
		std::vector<uint32_t> uniqueFoundRevComp;  ///< per collection, reverse complement k-mers unique to the collection

		uint32_t winnerIdx{NO_WINNER};
		bool winnerRevComp{false};
		double bestFrac{0};

		[[nodiscard]] bool hasWinner() const;
		/**
		 * @brief Clear the k-mers and the winner, ready for the next read
		 */
		void reset();
	};

	/**
	 * @param collections collection name to its hashed k-mers, must outlive the matcher
	 * @param klen the k-mer length the collections were hashed with
	 * @param useInvertedIndex look k-mers up in a k-mer -> collection(s) index, otherwise check each collection
	 * @param countUnique also count uniqueFound and uniqueFoundRevComp when not using the index
	 */
	KmerCollectionReadMatcher(const Collections & collections, uint32_t klen, bool useInvertedIndex = true, bool countUnique = true);

	/**
	 * @brief Collection names in index order
	 */
	[[nodiscard]] const VecStr & names() const;

	/**
	 * @brief Whether this matcher was built for these collections and settings, checked by address, collection count and total k-mers
	 */
	[[nodiscard]] bool usable(const Collections & collections, uint32_t klen, bool useInvertedIndex, bool countUnique) const;

	/**
	 * @brief A Match sized for these collections
	 */
	[[nodiscard]] Match genMatch() const;

	/**
	 * @brief Add a sequence's hashed k-mers to match.kmers and, if revComp, its reverse complement's to match.kmersRevComp
	 *
	 * Sequences no longer than the k-mer length add nothing. Call more than once to pool k-mers, e.g. both mates of a pair.
	 */
	void addKmers(const std::string & seq, bool revComp, Match & match) const;

	/**
	 * @brief Set the found counts from match.kmers and match.kmersRevComp
	 */
	void countFound(Match & match) const;

	/**
	 * @brief Pick the collection with the highest fraction of the read's k-mers in found or foundRevComp
	 *
	 * Fractions are over match.kmers.size(). Collections are checked in index order, forward before reverse complement,
	 * and ties go to the first. Counts below any cut offs should be zeroed before calling.
	 */
	void determineWinner(Match & match) const;

	/**
	 * @brief Number of collections whose unique k-mer count, forward or reverse complement, is at least minCount
	 */
	[[nodiscard]] static uint32_t countCollectionsHit(const Match & match, uint32_t minCount);

private:
	const Collections & collections_;
	uint32_t klen_;
	bool useInvertedIndex_;
	bool countUnique_;
	size_t totalKmers_{0};
	SimpleKmerHash hasher_;
	std::shared_ptr<const KmerSetInvertedIndex> index_;
	VecStr names_;
	//only built when not using the index
	Collections uniqueCollections_;
};

}  // namespace njhseq
