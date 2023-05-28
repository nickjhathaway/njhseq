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


namespace njhseq {

class KmerGatherer{
public:
	struct KmerGathererPars {
		KmerGathererPars(uint32_t kmerLength,
				bool noRevComp, uint32_t numThreads,
				std::set<char> allowableCharacters);
		KmerGathererPars();
		uint32_t kmerLength_{19};
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


	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHash(const std::vector<bfs::path> & twobitFnps) const;
	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHashWithFilters(const std::vector<bfs::path> & twobitFnps) const;
	[[nodiscard]] std::unordered_map<std::string, std::set<uint64_t>> getUniqueKmersSetHashWithFiltersFromFastas(const std::vector<bfs::path> & fastaTwoBit) const;

};


}  // namespace njhseq




