#include <catch_amalgamated.hpp>
#include "njhseq/alignment/aligner/aligner.hpp"
#include <functional>

using namespace njhseq;

namespace {
comparison compare(aligner & alignerObj, const seqInfo & ref, const seqInfo & query) {
	alignerObj.alignRegGlobal(ref, query);
	return alignerObj.compareAlignment(ref, query, false);
}
}  // namespace

TEST_CASE("compareAlignment counts an internal indel the same in either direction", "[aligner][compareAlignment]") {
	aligner alignerObj(200, gapScoringParameters(5, 1));
	seqInfo full("full", "ACGTTGCAAGGCTTACCGATAGCTAGGCATCGATCGGATC");
	//same sequence with the C at position 20 removed
	seqInfo deleted("deleted", "ACGTTGCAAGGCTTACCGATAGTAGGCATCGATCGGATC");

	auto fullRef = compare(alignerObj, full, deleted);
	auto deletedRef = compare(alignerObj, deleted, full);

	REQUIRE(fullRef.distances_.alignmentGaps_.size() == 1);
	REQUIRE(deletedRef.distances_.alignmentGaps_.size() == 1);
	REQUIRE(fullRef.distances_.getNumOfEvents(true) == deletedRef.distances_.getNumOfEvents(true));
	REQUIRE(fullRef.distances_.eventBasedIdentity_ == Catch::Approx(deletedRef.distances_.eventBasedIdentity_));
	REQUIRE(fullRef.distances_.eventBasedIdentity_ < 1.0);
}

TEST_CASE("compareAlignment of identical sequences has no events", "[aligner][compareAlignment]") {
	aligner alignerObj(200, gapScoringParameters(5, 1));
	seqInfo first("first", "ACGTTGCAAGGCTTACCGATAGCTAGGCATCGATCGGATC");
	seqInfo second("second", first.seq_);
	auto comp = compare(alignerObj, first, second);
	REQUIRE(comp.distances_.getNumOfEvents(true) == 0);
	REQUIRE(comp.distances_.eventBasedIdentity_ == Catch::Approx(1.0));
}

TEST_CASE("aligning sequences longer than the aligner was made for gives the same alignment as an aligner made large enough", "[aligner]") {
	const std::string piece = "ACGTTGCAAGGCTTACCGATAGCTAGGCATCGATCGGATC";
	std::string longer;
	for (uint32_t i = 0; i < 20; ++i) {
		longer += piece;
	}
	seqInfo ref("ref", longer);
	//same sequence with a base removed, one changed and two added to the end
	seqInfo query("query", longer.substr(0, 300) + longer.substr(301) + "GT");
	query.seq_[500] = 'A' == query.seq_[500] ? 'C' : 'A';

	auto sameFor = [&ref, &query](const std::function<void(aligner &, const seqInfo &, const seqInfo &)> & align) {
		for (const auto & seqs : {std::make_pair(ref, query), std::make_pair(query, ref)}) {
			aligner large(1000, gapScoringParameters(5, 1));
			aligner small(50, gapScoringParameters(5, 1));
			align(large, seqs.first, seqs.second);
			align(small, seqs.first, seqs.second);
			REQUIRE(large.alignObjectA_.seqBase_.seq_ == small.alignObjectA_.seqBase_.seq_);
			REQUIRE(large.alignObjectB_.seqBase_.seq_ == small.alignObjectB_.seqBase_.seq_);
			REQUIRE(large.parts_.score_ == small.parts_.score_);
		}
	};
	SECTION("global") {
		sameFor([](aligner & alignerObj, const seqInfo & first, const seqInfo & second) { alignerObj.alignCacheGlobal(first, second); });
	}
	SECTION("global diagonal") {
		sameFor([](aligner & alignerObj, const seqInfo & first, const seqInfo & second) { alignerObj.alignCacheGlobalDiag(first, second); });
	}
	SECTION("global without internal gaps") {
		sameFor([](aligner & alignerObj, const seqInfo & first, const seqInfo & second) { alignerObj.alignCacheGlobalNoInternalGaps(first, second); });
	}
	SECTION("local") {
		sameFor([](aligner & alignerObj, const seqInfo & first, const seqInfo & second) { alignerObj.alignCacheLocal(first, second); });
	}
}
