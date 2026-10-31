#pragma once

/*
 * KmerCollectionIO.hpp
 *
 *  Created on: Oct 31, 2026
 *      Author: nick
 */

#include "njhseq/common.h"
#include "njhseq/objects/kmer/SimpleKmerHash.hpp"

namespace njhseq {

/**
 * @brief Reads and writes k-mer collection tables
 *
 * A k-mer collection table is a headerless, tab delimited table of 1) collection name, 2) k-mer, one row per k-mer.
 * Hashed k-mers (SimpleKmerHash) are written out as their k-mer strings and hashed again when read in.
 */
class KmerCollectionIO {
public:
	using HashedCollections = std::unordered_map<std::string, std::unordered_set<uint64_t>>;

	/**
	 * @brief Write one collection's hashed k-mers, one row per k-mer in the container's order
	 * @param out the stream to write to
	 * @param name the collection name for the first column
	 * @param kmers the hashed k-mers
	 */
	template<typename KMERS>
	static void writeHashedCollection(std::ostream & out, const std::string & name, const KMERS & kmers);

	/**
	 * @brief Write every collection's hashed k-mers, collections in the map's order
	 * @param out the stream to write to
	 * @param collections collection name to its hashed k-mers
	 */
	template<typename COLLECTIONS>
	static void writeHashedCollections(std::ostream & out, const COLLECTIONS & collections);

	/**
	 * @brief The k-mer length of a table, from its first row
	 */
	[[nodiscard]] static uint32_t getKmerLength(const bfs::path & tableFnp);

	/**
	 * @brief Read a table into collection name to its hashed k-mers
	 */
	[[nodiscard]] static HashedCollections readHashedCollections(const bfs::path & tableFnp);

	/**
	 * @brief Read a table into collection name to its hashed k-mers
	 * @param kmerLength set to the length of the k-mer in the last row, or 0 if the table is empty
	 */
	[[nodiscard]] static HashedCollections readHashedCollections(const bfs::path & tableFnp, uint32_t & kmerLength);

	/**
	 * @brief Read every hashed k-mer in a table, ignoring which collection it is in
	 */
	[[nodiscard]] static std::unordered_set<uint64_t> readHashedKmers(const bfs::path & tableFnp);
};

template<typename KMERS>
void KmerCollectionIO::writeHashedCollection(std::ostream & out, const std::string & name, const KMERS & kmers) {
	SimpleKmerHash hasher;
	for (const auto & kmer : kmers) {
		out << name << "\t" << hasher.reverseHash(kmer) << "\n";
	}
}

template<typename COLLECTIONS>
void KmerCollectionIO::writeHashedCollections(std::ostream & out, const COLLECTIONS & collections) {
	for (const auto & collection : collections) {
		writeHashedCollection(out, collection.first, collection.second);
	}
}

}  // namespace njhseq
