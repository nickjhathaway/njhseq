#include <catch_amalgamated.hpp>
#include "njhseq/objects/kmer/KmerCollectionIO.hpp"
#include "njhseq/IO/OutputStream.hpp"

using namespace njhseq;

namespace {

struct TempDir {
	TempDir() : dir_(bfs::temp_directory_path() / bfs::unique_path("njhseqTest-%%%%-%%%%")) {
		bfs::create_directories(dir_);
	}
	~TempDir() {
		bfs::remove_all(dir_);
	}
	bfs::path dir_;
};

KmerCollectionIO::HashedCollections exampleCollections() {
	SimpleKmerHash hasher;
	KmerCollectionIO::HashedCollections ret;
	for (const auto & kmer : {"ACGTACGTAC", "TTGACCGGTA", "GGGCATTACA"}) {
		ret["A"].emplace(hasher.hash(kmer));
	}
	for (const auto & kmer : {"CATCATCATG", "ACGTACGTAC"}) {
		ret["B"].emplace(hasher.hash(kmer));
	}
	return ret;
}

}  // namespace

TEST_CASE("KmerCollectionIO writes one row of collection and k-mer per k-mer", "[kmer][KmerCollectionIO]") {
	SimpleKmerHash hasher;
	std::vector<uint64_t> kmers{hasher.hash("ACGTACGTAC"), hasher.hash("TTGACCGGTA")};
	std::stringstream out;
	KmerCollectionIO::writeHashedCollection(out, "set1", kmers);
	REQUIRE(out.str() == "set1\tACGTACGTAC\nset1\tTTGACCGGTA\n");

	std::map<std::string, std::vector<uint64_t>> collections{{"a", {hasher.hash("CATCATCATG")}}, {"b", kmers}};
	std::stringstream outAll;
	KmerCollectionIO::writeHashedCollections(outAll, collections);
	REQUIRE(outAll.str() == "a\tCATCATCATG\nb\tACGTACGTAC\nb\tTTGACCGGTA\n");
}

TEST_CASE("KmerCollectionIO reads back what it writes, plain or gzipped", "[kmer][KmerCollectionIO]") {
	TempDir tmp;
	auto collections = exampleCollections();
	for (const auto & name : {"kmers.tsv", "kmers.tsv.gz"}) {
		auto fnp = tmp.dir_ / name;
		{
			OutputStream out(fnp);
			KmerCollectionIO::writeHashedCollections(out, collections);
		}
		REQUIRE(KmerCollectionIO::readHashedCollections(fnp) == collections);
		uint32_t kmerLength = 0;
		REQUIRE(KmerCollectionIO::readHashedCollections(fnp, kmerLength) == collections);
		REQUIRE(10 == kmerLength);
		REQUIRE(10 == KmerCollectionIO::getKmerLength(fnp));
		std::unordered_set<uint64_t> allKmers;
		for (const auto & collection : collections) {
			allKmers.insert(collection.second.begin(), collection.second.end());
		}
		REQUIRE(KmerCollectionIO::readHashedKmers(fnp) == allKmers);
	}
}
