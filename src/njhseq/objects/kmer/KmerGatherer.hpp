#pragma once

/*
 * KmerGatherer.hpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */


#include <utility>

#include "njhseq/common.h"
#include "njhseq/programUtils/seqSetUp.hpp"
#include "njhseq/objects/kmer/SimpleKmerHash.hpp"


namespace njhseq {

/**
 * @brief Gathers k-mers from fasta files, 2bit files or in-memory sequences
 *
 * A k-mer collection is the distinct k-mers of one named group of sequences, where the group is a whole fasta file,
 * a whole 2bit file or one key of an in-memory map. Collections are returned as a map of collection name to its k-mers.
 * K-mers from both strands are gathered unless KmerGathererPars::noRevComp_ is set.
 */
class KmerGatherer{
public:
	struct KmerGathererPars {
		KmerGathererPars(uint32_t kmerLength,
				bool noRevComp, uint32_t numThreads,
				std::set<char> allowableCharacters);
		KmerGathererPars();
		uint32_t kmerLength_{19}; ///< k-mer length, hashed k-mers can be at most 19
		uint32_t kmerLengthForEntropyCalc_{1}; ///< sub-k-mer length used when computing a k-mer's entropy for the filter
		bool noRevComp_{false}; ///< only gather k-mers from the forward strand
		uint32_t numThreads_{1};
		std::set<char> allowableCharacters_{'A', 'C', 'G', 'T'}; ///< filter, k-mers must only contain these characters
		double entropyFilter_{1.95}; ///< filter, k-mers must have an entropy above this (exclusive), 0 or below turns the entropy filter off
		bool allUpper_ = false; ///< upper case sequences before gathering, otherwise lower case k-mers stay distinct and fail the allowable character filter

		void setOptions(seqSetUp & setUp);
	};

	explicit KmerGatherer(KmerGathererPars pars);
	KmerGathererPars pars_;

	/**
	 * @brief Count every k-mer occurrence across all sequences in a fasta file
	 * @param fastaFnp the fasta file
	 * @return k-mer to number of occurrences
	 */
	[[nodiscard]] std::unordered_map<std::string, uint32_t> countKmers(const bfs::path & fastaFnp) const;

	/**
	 * @brief Gather one collection of string k-mers per fasta file, unfiltered
	 * @param fastaFnps the fasta files
	 * @return fasta file path to its distinct k-mers
	 */
	[[nodiscard]] std::unordered_map<std::string, std::set<std::string>> getKmerCollectionsFromFastas(const std::vector<bfs::path> & fastaFnps) const;

	/**
	 * @brief Gather one collection of hashed k-mers per 2bit file, keeping only k-mers that pass passesFilters()
	 * @param twobitFnps the 2bit files
	 * @return 2bit file path to its distinct hashed k-mers
	 */
	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getFilteredHashedKmerCollectionsFromTwobits(const std::vector<bfs::path> & twobitFnps) const;

	/**
	 * @brief Gather one collection of hashed k-mers per fasta file, keeping only k-mers that pass passesFilters()
	 * @param fastaFnps the fasta files
	 * @return fasta file path to its distinct hashed k-mers
	 */
	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getFilteredHashedKmerCollectionsFromFastas(const std::vector<bfs::path> & fastaFnps) const;

	/**
	 * @brief Gather one collection of hashed k-mers per key of an in-memory map of sequences, unfiltered
	 * @param seqs collection name to its sequences, anything getSeqBase() accepts
	 * @return collection name to its distinct hashed k-mers
	 */
	template<typename T>
	[[nodiscard]] std::unordered_map<std::string, std::unordered_set<uint64_t>> getHashedKmerCollectionsFromSeqs(const std::unordered_map<std::string, std::vector<T>> & seqs) const;

	/**
	 * @brief Whether a k-mer only contains allowable characters and has an entropy above the entropy filter, the entropy check is skipped when the filter is 0 or below
	 * @param kmer the k-mer to check
	 * @return true if the k-mer passes both filters
	 */
	[[nodiscard]] bool passesFilters(const std::string & kmer) const;

	/**
	 * @brief Reduce each collection to the k-mers found in no other collection
	 * @param kmersPerCollection collection name to its distinct k-mers, a std::map or std::unordered_map of std::set or std::unordered_set
	 * @return same type as the input, every input collection is present even if it ends up empty
	 */
	template<typename COLLECTIONS>
	[[nodiscard]] static COLLECTIONS filterToKmersUniqueToEachCollection(const COLLECTIONS & kmersPerCollection);

	/**
	 * @brief The k-mers found in more than one collection, the complement of filterToKmersUniqueToEachCollection()
	 * @param kmersPerCollection collection name to its distinct k-mers
	 * @return the shared k-mers, in the same container type as each input collection
	 */
	template<typename COLLECTIONS>
	[[nodiscard]] static typename COLLECTIONS::mapped_type getKmersSharedBetweenCollections(const COLLECTIONS & kmersPerCollection);

	/**
	 * @brief Merge collections into groups, each group's collection being the union of its members' collections
	 * @param kmersPerMember member name (e.g. file path) to its distinct k-mers
	 * @param membersPerGroup group name to the names of its members, every member must be in kmersPerMember
	 * @return group name to the union of its members' k-mers
	 */
	template<typename COLLECTIONS, typename MEMBERS>
	[[nodiscard]] static std::map<std::string, typename COLLECTIONS::mapped_type> condenseCollections(const COLLECTIONS & kmersPerMember, const MEMBERS & membersPerGroup);

private:
	struct TwobitFnpSeqNamePair{
		TwobitFnpSeqNamePair(bfs::path twoBit, std::string seqName):twoBit_(std::move(twoBit)), seqName_(std::move(seqName)){

		}
		TwobitFnpSeqNamePair()= default;
		bfs::path twoBit_;
		std::string seqName_;
	};

	SimpleKmerHash hasher_;

	/**
	 * @brief Throw if the k-mer length is too long for SimpleKmerHash
	 * @param funcName calling function, for the error message
	 */
	void checkKmerLengthForHashing(const std::string & funcName) const;

	/**
	 * @brief Fasta input options with lower case bases upper cased if pars_.allUpper_ is set
	 */
	[[nodiscard]] SeqIOOptions genFastaInOpts(const bfs::path & fastaFnp) const;

	/**
	 * @brief Every sequence name in every 2bit file
	 */
	[[nodiscard]] static std::vector<TwobitFnpSeqNamePair> getTwobitSeqPairs(const std::vector<bfs::path> & twobitFnps);

	/**
	 * @brief Add the string k-mers of a sequence and, unless pars_.noRevComp_, of its reverse complement
	 */
	void addKmers(const std::string & seq, std::set<std::string> & kmers) const;

	/**
	 * @brief Add the hashed k-mers of a sequence and, unless pars_.noRevComp_, of its reverse complement
	 * @param filter only add k-mers that pass passesFilters(), judged on the forward strand k-mer
	 */
	template<typename KMERS>
	void addHashedKmers(const std::string & seq, bool filter, KMERS & kmers) const;

	/**
	 * @brief Run addWorkKmers over each work item on pars_.numThreads_ threads, each thread fills its own collections that are then merged
	 * @param work the work items, e.g. files or collection names
	 * @param addWorkKmers called as addWorkKmers(item, collections) to add an item's k-mers to its collection
	 */
	template<typename COLLECTION, typename WORK, typename FUNC>
	[[nodiscard]] std::unordered_map<std::string, COLLECTION> gatherCollections(const std::vector<WORK> & work, const FUNC & addWorkKmers) const;

	/**
	 * @brief The number of collections each k-mer is found in
	 */
	template<typename COLLECTIONS>
	[[nodiscard]] static std::unordered_map<typename COLLECTIONS::mapped_type::value_type, uint32_t> countCollectionsPerKmer(const COLLECTIONS & kmersPerCollection);
};


template<typename KMERS>
void KmerGatherer::addHashedKmers(const std::string & seq, bool filter, KMERS & kmers) const {
	if (seq.size() < pars_.kmerLength_) {
		return;
	}
	for (uint32_t pos = 0; pos < seq.size() - pars_.kmerLength_ + 1; ++pos) {
		if (filter && !passesFilters(seq.substr(pos, pars_.kmerLength_))) {
			continue;
		}
		kmers.emplace(hasher_.hash(seq, pos, pars_.kmerLength_));
		if (!pars_.noRevComp_) {
			kmers.emplace(hasher_.revCompHash(seq, pos, pars_.kmerLength_));
		}
	}
}

template<typename COLLECTION, typename WORK, typename FUNC>
std::unordered_map<std::string, COLLECTION> KmerGatherer::gatherCollections(const std::vector<WORK> & work, const FUNC & addWorkKmers) const {
	njh::concurrent::LockableQueue<WORK> workQueue(work);
	std::unordered_map<std::string, COLLECTION> ret;
	std::mutex mut;
	std::function<void()> gather = [&workQueue, &addWorkKmers, &ret, &mut]() {
		WORK item;
		std::unordered_map<std::string, COLLECTION> current;
		while (workQueue.getVal(item)) {
			addWorkKmers(item, current);
		}
		std::lock_guard<std::mutex> lock(mut);
		for (const auto & collection : current) {
			ret[collection.first].insert(collection.second.begin(), collection.second.end());
		}
	};
	njh::concurrent::runVoidFunctionThreaded(gather, pars_.numThreads_);
	return ret;
}

template<typename T>
std::unordered_map<std::string, std::unordered_set<uint64_t>> KmerGatherer::getHashedKmerCollectionsFromSeqs(const std::unordered_map<std::string, std::vector<T>> & seqs) const {
	checkKmerLengthForHashing(__PRETTY_FUNCTION__);
	return gatherCollections<std::unordered_set<uint64_t>>(njh::getVecOfMapKeys(seqs),
		[this, &seqs](const std::string & name, std::unordered_map<std::string, std::unordered_set<uint64_t>> & current) {
			auto & kmers = current[name];
			for (const auto & seq : seqs.at(name)) {
				if (pars_.allUpper_) {
					auto upperSeq = getSeqBase(seq).seq_;
					njh::strToUpper(upperSeq);
					addHashedKmers(upperSeq, false, kmers);
				} else {
					addHashedKmers(getSeqBase(seq).seq_, false, kmers);
				}
			}
		});
}

template<typename COLLECTIONS>
std::unordered_map<typename COLLECTIONS::mapped_type::value_type, uint32_t> KmerGatherer::countCollectionsPerKmer(const COLLECTIONS & kmersPerCollection) {
	std::unordered_map<typename COLLECTIONS::mapped_type::value_type, uint32_t> ret;
	for (const auto & collection : kmersPerCollection) {
		for (const auto & kmer : collection.second) {
			++ret[kmer];
		}
	}
	return ret;
}

template<typename COLLECTIONS>
COLLECTIONS KmerGatherer::filterToKmersUniqueToEachCollection(const COLLECTIONS & kmersPerCollection) {
	auto collectionCounts = countCollectionsPerKmer(kmersPerCollection);
	COLLECTIONS ret;
	for (const auto & collection : kmersPerCollection) {
		auto & uniqueKmers = ret[collection.first];
		for (const auto & kmer : collection.second) {
			if (1 == collectionCounts.at(kmer)) {
				uniqueKmers.emplace(kmer);
			}
		}
	}
	return ret;
}

template<typename COLLECTIONS>
typename COLLECTIONS::mapped_type KmerGatherer::getKmersSharedBetweenCollections(const COLLECTIONS & kmersPerCollection) {
	typename COLLECTIONS::mapped_type ret;
	for (const auto & count : countCollectionsPerKmer(kmersPerCollection)) {
		if (count.second > 1) {
			ret.emplace(count.first);
		}
	}
	return ret;
}

template<typename COLLECTIONS, typename MEMBERS>
std::map<std::string, typename COLLECTIONS::mapped_type> KmerGatherer::condenseCollections(const COLLECTIONS & kmersPerMember, const MEMBERS & membersPerGroup) {
	std::map<std::string, typename COLLECTIONS::mapped_type> ret;
	for (const auto & group : membersPerGroup) {
		auto & kmers = ret[group.first];
		for (const auto & member : group.second) {
			const auto & memberKmers = kmersPerMember.at(member);
			kmers.insert(memberKmers.begin(), memberKmers.end());
		}
	}
	return ret;
}


}  // namespace njhseq
