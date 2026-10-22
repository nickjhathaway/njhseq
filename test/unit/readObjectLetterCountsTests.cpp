#include <catch_amalgamated.hpp>
#include "njhseq/objects/seqObjects/readObject.hpp"
#include "njhseq/objects/seqObjects/Paired/PairedRead.hpp"

#include <thread>

using namespace njhseq;

TEST_CASE("readObject doesn't carry a charCounter", "[readObject][letterCounts]") {
	REQUIRE(sizeof(readObject) < sizeof(charCounter));
}

TEST_CASE("readObject letter counts", "[readObject][letterCounts]") {
	readObject read(seqInfo("read", "AACGTT"));
	auto counts = read.getLetterCounts();
	REQUIRE(counts.chars_['A'] == 2);
	REQUIRE(counts.chars_['C'] == 1);
	REQUIRE(counts.chars_['G'] == 1);
	REQUIRE(counts.chars_['T'] == 2);
	REQUIRE(counts.fractions_['A'] == Catch::Approx(2.0 / 6));
}

TEST_CASE("readObject letter counts are weighted by read count", "[readObject][letterCounts]") {
	seqInfo seq("read", "AACG");
	seq.cnt_ = 3;
	readObject read(seq);
	auto counts = read.getLetterCounts();
	REQUIRE(counts.chars_['A'] == 6);
	REQUIRE(counts.fractions_['A'] == Catch::Approx(0.5));
}

TEST_CASE("readObject letter counts follow the current sequence", "[readObject][letterCounts]") {
	readObject read(seqInfo("read", "AAAA"));
	REQUIRE(read.getLetterCounts().chars_['A'] == 4);
	read.seqBase_.seq_ = "CCCC";
	auto counts = read.getLetterCounts();
	REQUIRE(counts.chars_['C'] == 4);
	REQUIRE(counts.chars_['A'] == 0);
}

TEST_CASE("readObject letter counts with an alphabet", "[readObject][letterCounts]") {
	readObject read(seqInfo("read", "ACGT"));
	auto counts = read.getLetterCounts(std::vector<char>{'A', 'C', 'G', 'T'});
	REQUIRE(counts.alphabet_ == std::vector<char>{'A', 'C', 'G', 'T'});
	REQUIRE(counts.fractions_['G'] == Catch::Approx(0.25));
}

TEST_CASE("readObject letter counts from many threads agree", "[readObject][letterCounts][threads]") {
	const readObject read(seqInfo("read", "AACGTT"));
	std::vector<double> aCounts(8, 0);
	std::vector<std::thread> threads;
	for (uint32_t t = 0; t < aCounts.size(); ++t) {
		threads.emplace_back([&read, &aCounts, t]() {
			aCounts[t] = read.getLetterCounts().chars_['A'];
		});
	}
	for (auto & th : threads) {
		th.join();
	}
	for (const auto count : aCounts) {
		REQUIRE(count == 2);
	}
}

TEST_CASE("readObject gc content", "[readObject][letterCounts]") {
	const readObject read(seqInfo("read", "GGCCAATT"));
	REQUIRE(read.getGCContent() == Catch::Approx(0.5));
}

TEST_CASE("PairedRead letter counts include the mate", "[PairedRead][letterCounts]") {
	SECTION("mate is reverse complemented when counting") {
		PairedRead read(seqInfo("first", "AAAA"), seqInfo("second", "AAAA"));
		auto counts = read.getLetterCounts();
		//mate AAAA is counted as its reverse complement TTTT
		REQUIRE(counts.chars_['A'] == 4);
		REQUIRE(counts.chars_['T'] == 4);
	}
	SECTION("with an alphabet keeps that alphabet") {
		PairedRead read(seqInfo("first", "AAAA"), seqInfo("second", "AAAA"));
		auto counts = read.getLetterCounts(std::vector<char>{'A', 'C', 'G', 'T', 'N'});
		REQUIRE(njh::in('N', counts.alphabet_));
		REQUIRE(counts.chars_['T'] == 4);
	}
}
