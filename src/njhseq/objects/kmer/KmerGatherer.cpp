/*
 * KmerGatherer.cpp
 *
 *  Created on: Jul 13, 2021
 *      Author: nick
 */

#include "KmerGatherer.hpp"
#include "njhseq/objects/kmer/KmerUtils.hpp"
#include "njhseq/objects/dataContainers.h"
#include "njhseq/IO/SeqIO.h"
#include <TwoBit.h>

#include <utility>

namespace njhseq {

KmerGatherer::KmerGathererPars::KmerGathererPars(uint32_t kmerLength,
		bool noRevComp, uint32_t numThreads, std::set<char> allowableCharacters) :
		kmerLength_(kmerLength), noRevComp_(noRevComp), numThreads_(numThreads), allowableCharacters_(std::move(
				allowableCharacters)) {

}
KmerGatherer::KmerGathererPars::KmerGathererPars() = default;

void KmerGatherer::KmerGathererPars::setOptions(seqSetUp &setUp) {
	setUp.setOption(kmerLength_, "--kmerLength", "kmer Length");
	setUp.setOption(kmerLengthForEntropyCalc_, "--kmerLengthForEntropyCalc", "kmer Length For Entropy Calculation");
	if(kmerLengthForEntropyCalc_ > kmerLength_){
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr("--kmerLengthForEntropyCalc: ", kmerLengthForEntropyCalc_, " can't be longer than kmer length, --kmerLength, for table, ", kmerLength_));
	}
	setUp.setOption(numThreads_, "--numThreads", "number of threads");
	setUp.setOption(allowableCharacters_, "--allowableCharacters",
			"Only count kmers with these allowable Characters");
	setUp.setOption(noRevComp_, "--noRevComp", "noRevComp");
	setUp.setOption(entropyFilter_, "--entropyFilter", "entropy Filter cut off, exclusive, will only keep kmers above this entropy level");
	setUp.setOption(allUpper_, "--changeSeqsToUpper", "change all sequences to upper case, otherwise they will filter off");
}

KmerGatherer::KmerGatherer(KmerGathererPars pars) :
		pars_(std::move(pars)) {

}

void KmerGatherer::checkKmerLengthForHashing(const std::string & funcName) const {
	if (pars_.kmerLength_ > 19) {
		std::stringstream ss;
		ss << funcName << ", error "
				<< "cannot do kmer lengths greater than 19" << "\n";
		throw std::runtime_error { ss.str() };
	}
}

SeqIOOptions KmerGatherer::genFastaInOpts(const bfs::path & fastaFnp) const {
	auto fastaInOpts = SeqIOOptions::genFastaIn(fastaFnp);
	if (pars_.allUpper_) {
		fastaInOpts.lowerCaseBases_ = "upper";
	}
	return fastaInOpts;
}

std::vector<KmerGatherer::TwobitFnpSeqNamePair> KmerGatherer::getTwobitSeqPairs(const std::vector<bfs::path> & twobitFnps) {
	std::vector<TwobitFnpSeqNamePair> pairs;
	for (const auto & twoBit : twobitFnps) {
		TwoBit::TwoBitFile tReader(twoBit);
		for (const auto & seqName : tReader.sequenceNames()) {
			pairs.emplace_back(twoBit, seqName);
		}
	}
	return pairs;
}

bool KmerGatherer::passesFilters(const std::string & kmer) const {
	if (!std::all_of(kmer.begin(), kmer.end(), [this](char base) { return njh::in(base, pars_.allowableCharacters_); })) {
		return false;
	}
	//an entropy filter of 0 or below turns the entropy check off
	return pars_.entropyFilter_ <= 0
			|| kmerInfo(kmer, pars_.kmerLengthForEntropyCalc_, false).computeKmerEntropy() > pars_.entropyFilter_;
}

void KmerGatherer::addKmers(const std::string & seq, std::set<std::string> & kmers) const {
	if (seq.size() < pars_.kmerLength_) {
		return;
	}
	for (uint32_t pos = 0; pos < seq.size() - pars_.kmerLength_ + 1; ++pos) {
		kmers.emplace(seq.substr(pos, pars_.kmerLength_));
	}
	if (!pars_.noRevComp_) {
		auto revComp = seqUtil::reverseComplement(seq, "DNA");
		for (uint32_t pos = 0; pos < revComp.size() - pars_.kmerLength_ + 1; ++pos) {
			kmers.emplace(revComp.substr(pos, pars_.kmerLength_));
		}
	}
}

std::unordered_map<std::string, uint32_t> KmerGatherer::countKmers(
		const bfs::path &fastaFnp) const {
	std::unordered_map < std::string, uint32_t > ret;
	SeqInput reader(genFastaInOpts(fastaFnp));
	reader.openIn();
	std::mutex countsMut;
	std::function < void() > countKmers =
			[&reader, this, &countsMut, &ret]() {
				seqInfo seq;
				while (reader.readNextReadLock(seq)) {
					std::unordered_map < std::string, uint32_t > countsCurrent;
					if (len(seq) >= pars_.kmerLength_) {
						for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
								++pos) {
							++countsCurrent[seq.seq_.substr(pos, pars_.kmerLength_)];
						}
						if (!pars_.noRevComp_) {
							seq.reverseComplementRead(false, true);
							for (uint32_t pos = 0; pos < len(seq) - pars_.kmerLength_ + 1;
									++pos) {
								++countsCurrent[seq.seq_.substr(pos, pars_.kmerLength_)];
							}
						}
					}
					{
						std::lock_guard < std::mutex > lock(countsMut);
						for (const auto &count : countsCurrent) {
							ret[count.first] += count.second;
						}
					}
				}
			};
	njh::concurrent::runVoidFunctionThreaded(countKmers, pars_.numThreads_);
	return ret;
}

std::unordered_map<std::string, std::set<std::string>> KmerGatherer::getKmerCollectionsFromFastas(
		const std::vector<bfs::path> &fastaFnps) const {
	return gatherCollections<std::set<std::string>>(fastaFnps,
		[this](const bfs::path & fastaFnp, std::unordered_map<std::string, std::set<std::string>> & current) {
			auto & kmers = current[fastaFnp.string()];
			SeqInput reader(genFastaInOpts(fastaFnp));
			reader.openIn();
			seqInfo seq;
			while (reader.readNextRead(seq)) {
				addKmers(seq.seq_, kmers);
			}
		});
}

std::unordered_map<std::string, std::set<uint64_t>> KmerGatherer::getFilteredHashedKmerCollectionsFromFastas(
		const std::vector<bfs::path> &fastaFnps) const {
	checkKmerLengthForHashing(__PRETTY_FUNCTION__);
	return gatherCollections<std::set<uint64_t>>(fastaFnps,
		[this](const bfs::path & fastaFnp, std::unordered_map<std::string, std::set<uint64_t>> & current) {
			auto & kmers = current[fastaFnp.string()];
			SeqInput reader(genFastaInOpts(fastaFnp));
			reader.openIn();
			seqInfo seq;
			while (reader.readNextRead(seq)) {
				addHashedKmers(seq.seq_, true, kmers);
			}
		});
}

std::unordered_map<std::string, std::set<uint64_t>> KmerGatherer::getFilteredHashedKmerCollectionsFromTwobits(
		const std::vector<bfs::path> &twobitFnps) const {
	checkKmerLengthForHashing(__PRETTY_FUNCTION__);
	return gatherCollections<std::set<uint64_t>>(getTwobitSeqPairs(twobitFnps),
		[this](const TwobitFnpSeqNamePair & pair, std::unordered_map<std::string, std::set<uint64_t>> & current) {
			auto & kmers = current[pair.twoBit_.string()];
			TwoBit::TwoBitFile tReader(pair.twoBit_);
			std::string buffer;
			tReader[pair.seqName_]->getSequence(buffer);
			if (pars_.allUpper_) {
				njh::strToUpper(buffer);
			}
			addHashedKmers(buffer, true, kmers);
		});
}

}  // namespace njhseq
