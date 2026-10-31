#include <catch_amalgamated.hpp>
#include "njhseq/objects/kmer/KmerGatherer.hpp"
#include "njhseq/objects/kmer/KmerSetInvertedIndex.hpp"
#include "njhseq/objects/kmer/KmerCollectionReadMatcher.hpp"
#include "njhseq/helpers/seqUtil.hpp"

using namespace njhseq;

namespace {

//reference: repeatedly set_difference sorted collections against every other collection
std::map<std::string, std::set<uint64_t>> uniqueByPairwiseDifference(const std::map<std::string, std::set<uint64_t>> & collections) {
	std::map<std::string, std::set<uint64_t>> ret;
	for (const auto & collection : collections) {
		std::set<uint64_t> uniqueKmers = collection.second;
		for (const auto & other : collections) {
			if (other.first == collection.first) {
				continue;
			}
			std::vector<uint64_t> notShared;
			std::set_difference(uniqueKmers.begin(), uniqueKmers.end(),
					other.second.begin(), other.second.end(),
					std::back_inserter(notShared));
			uniqueKmers = std::set<uint64_t>(notShared.begin(), notShared.end());
		}
		ret[collection.first] = uniqueKmers;
	}
	return ret;
}

template<typename COLLECTIONS>
std::map<std::string, std::set<uint64_t>> toSorted(const COLLECTIONS & collections) {
	std::map<std::string, std::set<uint64_t>> ret;
	for (const auto & collection : collections) {
		ret[collection.first] = std::set<uint64_t>(collection.second.begin(), collection.second.end());
	}
	return ret;
}

template<typename COLLECTIONS>
COLLECTIONS fromSorted(const std::map<std::string, std::set<uint64_t>> & collections) {
	COLLECTIONS ret;
	for (const auto & collection : collections) {
		ret[collection.first].insert(collection.second.begin(), collection.second.end());
	}
	return ret;
}

const std::map<std::string, std::set<uint64_t>> exampleCollections{
		{"a", {1, 2, 3, 4}},
		{"b", {3, 4, 5}},
		{"c", {4, 6, 7}},
		{"d", {}}};

const std::map<std::string, std::set<uint64_t>> exampleUnique{
		{"a", {1, 2}},
		{"b", {5}},
		{"c", {6, 7}},
		{"d", {}}};

}  // namespace

TEMPLATE_TEST_CASE("filterToKmersUniqueToEachCollection keeps k-mers found in only one collection", "[kmer][KmerGatherer]",
		(std::unordered_map<std::string, std::set<uint64_t>>),
		(std::unordered_map<std::string, std::unordered_set<uint64_t>>),
		(std::map<std::string, std::set<uint64_t>>),
		(std::map<std::string, std::unordered_set<uint64_t>>)) {
	auto input = fromSorted<TestType>(exampleCollections);
	TestType unique = KmerGatherer::filterToKmersUniqueToEachCollection(input);
	REQUIRE(toSorted(unique) == exampleUnique);
}

TEST_CASE("filterToKmersUniqueToEachCollection matches pairwise set_difference", "[kmer][KmerGatherer]") {
	std::mt19937 gen(7);
	std::uniform_int_distribution<uint64_t> kmerDist(0, 400);
	std::map<std::string, std::set<uint64_t>> collections;
	for (uint32_t c = 0; c < 6; ++c) {
		auto & kmers = collections[njh::pasteAsStr("c", c)];
		for (uint32_t i = 0; i < 150; ++i) {
			kmers.emplace(kmerDist(gen));
		}
	}
	auto expected = uniqueByPairwiseDifference(collections);
	REQUIRE(KmerGatherer::filterToKmersUniqueToEachCollection(collections) == expected);
	REQUIRE(toSorted(KmerGatherer::filterToKmersUniqueToEachCollection(
			fromSorted<std::unordered_map<std::string, std::unordered_set<uint64_t>>>(collections))) == expected);
}

TEST_CASE("filterToKmersUniqueToEachCollection with one collection keeps all its k-mers", "[kmer][KmerGatherer]") {
	std::map<std::string, std::set<uint64_t>> collections{{"only", {1, 2, 3}}};
	REQUIRE(KmerGatherer::filterToKmersUniqueToEachCollection(collections) == collections);
}

TEST_CASE("getKmersSharedBetweenCollections keeps k-mers found in more than one collection", "[kmer][KmerGatherer]") {
	REQUIRE(KmerGatherer::getKmersSharedBetweenCollections(exampleCollections) == std::set<uint64_t>{3, 4});
	auto unordered = fromSorted<std::unordered_map<std::string, std::unordered_set<uint64_t>>>(exampleCollections);
	REQUIRE(KmerGatherer::getKmersSharedBetweenCollections(unordered) == std::unordered_set<uint64_t>{3, 4});
}

TEST_CASE("KmerSetInvertedIndex tallyAllAndUnique matches tallying all and unique-to-collection k-mers", "[kmer][KmerSetInvertedIndex]") {
	auto all = fromSorted<KmerSetInvertedIndex::SetMap>(exampleCollections);
	auto unique = KmerGatherer::filterToKmersUniqueToEachCollection(all);
	auto index = KmerSetInvertedIndex::get(all);
	std::unordered_map<std::string, uint32_t> allCounts;
	std::unordered_map<std::string, uint32_t> uniqueCounts;
	std::unordered_map<std::string, uint32_t> expectedAll;
	std::unordered_map<std::string, uint32_t> expectedUnique;
	for (uint64_t kmer = 0; kmer < 10; ++kmer) {
		index->tallyAllAndUnique(kmer, allCounts, uniqueCounts);
		for (const auto & collection : all) {
			if (njh::in(kmer, collection.second)) {
				++expectedAll[collection.first];
			}
			if (njh::in(kmer, unique.at(collection.first))) {
				++expectedUnique[collection.first];
			}
		}
	}
	REQUIRE(allCounts == expectedAll);
	REQUIRE(uniqueCounts == expectedUnique);
}

TEST_CASE("condenseCollections unions each group's member collections", "[kmer][KmerGatherer]") {
	std::unordered_map<std::string, std::vector<std::string>> membersPerGroup{
			{"ab", {"a", "b"}},
			{"d", {"d"}}};
	std::map<std::string, std::set<uint64_t>> expected{
			{"ab", {1, 2, 3, 4, 5}},
			{"d", {}}};
	REQUIRE(KmerGatherer::condenseCollections(exampleCollections, membersPerGroup) == expected);
	REQUIRE_THROWS(KmerGatherer::condenseCollections(exampleCollections,
			std::map<std::string, std::set<std::string>>{{"missing", {"z"}}}));
}

TEST_CASE("getHashedKmerCollectionsFromSeqs gathers both strands per collection", "[kmer][KmerGatherer]") {
	KmerGatherer::KmerGathererPars pars;
	pars.kmerLength_ = 5;
	pars.numThreads_ = 2;
	KmerGatherer gatherer(pars);
	std::unordered_map<std::string, std::vector<seqInfo>> seqs{
			{"one", {seqInfo("s1", "ACGTTGCAAG"), seqInfo("s2", "TTTACG")}},
			{"short", {seqInfo("s3", "AC")}}};
	auto collections = gatherer.getHashedKmerCollectionsFromSeqs(seqs);

	SimpleKmerHash hasher;
	std::unordered_set<uint64_t> expected;
	for (const auto & seq : seqs.at("one")) {
		for (const auto & strand : {seq.seq_, seqUtil::reverseComplement(seq.seq_, "DNA")}) {
			for (uint32_t pos = 0; pos + 5 <= strand.size(); ++pos) {
				expected.emplace(hasher.hash(strand.substr(pos, 5)));
			}
		}
	}
	REQUIRE(collections.size() == 2);
	REQUIRE(collections.at("one") == expected);
	REQUIRE(collections.at("short").empty());
}

TEST_CASE("getHashedKmerCollectionsFromSeqs upper cases sequences when allUpper_ is set", "[kmer][KmerGatherer]") {
	KmerGatherer::KmerGathererPars pars;
	pars.kmerLength_ = 5;
	std::unordered_map<std::string, std::vector<seqInfo>> lower{{"c", {seqInfo("s1", "acgtTGCAag")}}};
	std::unordered_map<std::string, std::vector<seqInfo>> upper{{"c", {seqInfo("s1", "ACGTTGCAAG")}}};
	REQUIRE(KmerGatherer(pars).getHashedKmerCollectionsFromSeqs(lower) != KmerGatherer(pars).getHashedKmerCollectionsFromSeqs(upper));
	pars.allUpper_ = true;
	REQUIRE(KmerGatherer(pars).getHashedKmerCollectionsFromSeqs(lower) == KmerGatherer(pars).getHashedKmerCollectionsFromSeqs(upper));
}

TEST_CASE("KmerCollectionReadMatcher counts all and unique k-mers and picks the best collection", "[kmer][KmerCollectionReadMatcher]") {
	//two overlapping targets and one separate one, collections built from the forward strand
	std::string overlapA = "ACGTTGCAAGCTTAGGCATCGATCGGATCCTAGCTAGCAATCG";
	std::string overlapB = overlapA.substr(15) + "TTGACCGGTAACGTAGCTTAGCAGGT";
	std::string separate = "GGATTCCAAGTCGTAGCATGCATTCGACGATCGTAGGCTAATC";
	KmerGatherer::KmerGathererPars pars;
	pars.kmerLength_ = 11;
	pars.noRevComp_ = true;
	std::unordered_map<std::string, std::vector<seqInfo>> refs{
			{"A", {seqInfo("A", overlapA)}},
			{"B", {seqInfo("B", overlapB)}},
			{"C", {seqInfo("C", separate)}}};
	auto collections = KmerGatherer(pars).getHashedKmerCollectionsFromSeqs(refs);
	auto unique = KmerGatherer::filterToKmersUniqueToEachCollection(collections);

	for (bool useIndex : {true, false}) {
		KmerCollectionReadMatcher matcher(collections, 11, useIndex);
		auto match = matcher.genMatch();
		//a read of A, given as its reverse complement
		matcher.addKmers(seqUtil::reverseComplement(overlapA, "DNA"), true, match);
		matcher.countFound(match);
		SimpleKmerHash hasher;
		for (uint32_t idx = 0; idx < matcher.names().size(); ++idx) {
			const auto & name = matcher.names()[idx];
			uint32_t expectedRevComp = 0;
			uint32_t expectedUniqueRevComp = 0;
			for (const auto & kmer : match.kmersRevComp) {
				expectedRevComp += njh::in(kmer, collections.at(name));
				expectedUniqueRevComp += njh::in(kmer, unique.at(name));
			}
			REQUIRE(match.foundRevComp[idx] == expectedRevComp);
			REQUIRE(match.uniqueFoundRevComp[idx] == expectedUniqueRevComp);
		}
		matcher.determineWinner(match);
		REQUIRE(match.hasWinner());
		REQUIRE(matcher.names()[match.winnerIdx] == "A");
		REQUIRE(match.winnerRevComp);
		REQUIRE(match.bestFrac == 1);
		//B shares k-mers with A but has none of A's unique ones, so only A is hit
		REQUIRE(KmerCollectionReadMatcher::countCollectionsHit(match, 1) == 1);

		match.reset();
		REQUIRE_FALSE(match.hasWinner());
		matcher.addKmers("ACGT", true, match);
		REQUIRE(match.kmers.empty());

		REQUIRE(matcher.usable(collections, 11, useIndex, true));
		REQUIRE_FALSE(matcher.usable(collections, 13, useIndex, true));
		auto changed = collections;
		REQUIRE_FALSE(matcher.usable(changed, 11, useIndex, true));
	}
}

TEST_CASE("passesFilters rejects disallowed characters and low entropy", "[kmer][KmerGatherer]") {
	KmerGatherer::KmerGathererPars pars;
	pars.kmerLength_ = 8;
	pars.entropyFilter_ = 1.5;
	KmerGatherer gatherer(pars);
	REQUIRE(gatherer.passesFilters("ACGTACGT"));
	REQUIRE_FALSE(gatherer.passesFilters("ACGNACGT"));
	REQUIRE_FALSE(gatherer.passesFilters("acgtacgt"));
	REQUIRE_FALSE(gatherer.passesFilters("AAAAAAAT"));
}

TEST_CASE("passesFilters skips the entropy check when the entropy filter is 0", "[kmer][KmerGatherer]") {
	KmerGatherer::KmerGathererPars pars;
	pars.kmerLength_ = 8;
	pars.entropyFilter_ = 0;
	KmerGatherer gatherer(pars);
	REQUIRE(gatherer.passesFilters("AAAAAAAA"));
	REQUIRE_FALSE(gatherer.passesFilters("AAAANAAA"));
}
