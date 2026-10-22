#include <catch_amalgamated.hpp>
#include "njhseq/alignment/aligner/aligner.hpp"

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
