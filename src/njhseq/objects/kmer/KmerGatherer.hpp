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

class KmerGatherer{
public:
	struct KmerGathererPars {
		KmerGathererPars(uint32_t kmerLength,
				bool noRevComp, uint32_t numThreads,
				std::set<char> allowableCharacters);
		KmerGathererPars();
		uint32_t kmerLength_{19};
		uint32_t kmerLengthForEntropyCalc_{1};
		bool noRevComp_{false};
		uint32_t numThreads_{1};
		std::set<char> allowableCharacters_{'A', 'C', 'G', 'T'};
		double entropyFilter_{1.95};
		bool allUpper_ = false;

		void setOptions(seqSetUp & setUp);
	};

	explicit KmerGatherer(KmerGathererPars pars);
	KmerGathererPars pars_;

	[[nodiscard]] std::unordered_map<std::string, uint32_t> countGenomeKmers(const bfs::path & genomeFnp) const;

	[[nodiscard]] std::unordered_set<std::string> getUniqueKmers(const bfs::path & genomeFnp) const;
	[[nodiscard]] std::set<std::string> getUniqueKmersSet(const bfs::path & genomeFnp) const;

	struct TwobitFnpSeqNamePair{
		TwobitFnpSeqNamePair(bfs::path twoBit, std::string seqName):twoBit_(std::move(twoBit)), seqName_(std::move(seqName)){

		}
		TwobitFnpSeqNamePair()= default;
		bfs::path twoBit_;
		std::string seqName_;
	};
	[[nodiscard]] std::unordered_map<std::string, std::set<std::string>> getUniqueKmersSet(const std::vector<bfs::path> & twobitFnps) const;

	template<typename T>
	std::unordered_map<std::string, std::unordered_set<uint64_t>> getUniqueKmersSetHash(const std::unordered_map<std::string, std::vector<T>> & seqs) const;
	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHash(const std::vector<bfs::path> & twobitFnps) const;
	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHashWithFilters(const std::vector<bfs::path> & twobitFnps) const;
	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHashWithFiltersFromFastas(const std::vector<bfs::path> & fastaTwoBit) const;

};


template<typename T>
std::unordered_map<std::string, std::unordered_set<uint64_t>> KmerGatherer::getUniqueKmersSetHash(const std::unordered_map<std::string, std::vector<T>> & seqs) const {
	if (pars_.kmerLength_ > 19) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "cannot do kmer lengths greater than 19" << "\n";
		throw std::runtime_error { ss.str() };
	}

	VecStr names = njh::getVecOfMapKeys(seqs);
	njh::concurrent::LockableQueue < std::string > namesQueue(names);

	std::unordered_map<std::string, std::unordered_set<uint64_t>> ret;
	std::mutex mut;

	std::function < void() > gatherKmers = [&namesQueue, this, &ret, &mut, &seqs]() {
		SimpleKmerHash hasher;
		std::string name;
		std::unordered_map<std::string, std::set<uint64_t>> current;
		while (namesQueue.getVal(name)) {
			std::set < uint64_t > genomeKmersCurrent;
			for(const auto & seq : seqs.at(name)) {
				for (uint32_t pos = 0; pos < len(getSeqBase(seq).seq_) - pars_.kmerLength_ + 1; ++pos) {
					genomeKmersCurrent.emplace(
							hasher.hash(getSeqBase(seq).seq_.substr(pos, pars_.kmerLength_)));
				}
				if (!pars_.noRevComp_) {
					std::string revComp = seqUtil::reverseComplement(getSeqBase(seq).seq_, "DNA");
					for (uint32_t pos = 0; pos < len(revComp) - pars_.kmerLength_ + 1;
							++pos) {
						genomeKmersCurrent.emplace(
								hasher.hash(revComp.substr(pos, pars_.kmerLength_)));
							}
				}
			}
			current[name].insert(genomeKmersCurrent.begin(),
					genomeKmersCurrent.end());
		}
		{
			std::lock_guard < std::mutex > lock(mut);
			for (const auto &kmerSet : current) {
				ret[kmerSet.first].insert(kmerSet.second.begin(), kmerSet.second.end());
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(gatherKmers, pars_.numThreads_);
	return ret;
}




}  // namespace njhseq




