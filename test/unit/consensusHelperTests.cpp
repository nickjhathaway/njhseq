#include <catch_amalgamated.hpp>
#include "njhseq/helpers/consensusHelper.hpp"

using njhseq::aligner;
using njhseq::consensusHelper;
using njhseq::gapScoringParameters;
using njhseq::seqInfo;
using njhseq::substituteMatrix;

namespace {
//the consensus of reads built off of a given starting sequence
std::string consensusFrom(const std::string & start, const std::vector<std::pair<std::string, uint32_t>> & seqCounts) {
	std::vector<seqInfo> reads;
	for (const auto & seqCount : seqCounts) {
		for (uint32_t i = 0; i < seqCount.second; ++i) {
			reads.emplace_back("read", seqCount.first, std::vector<uint8_t>(seqCount.first.size(), 30));
			reads.back().cnt_ = 1;
		}
	}
	seqInfo startSeq("start", start, std::vector<uint8_t>(start.size(), 30));
	startSeq.cnt_ = reads.size();
	aligner alignerObj(100, gapScoringParameters(5, 1, 0, 0, 0, 0), substituteMatrix::createDegenScoreMatrix(2, -2));
	return consensusHelper::buildConsensus(startSeq, reads, [](const seqInfo & seq) -> const seqInfo & { return seq; }, alignerObj).seq_;
}
const std::string shortEnd = "GATGAGTATACAATGCATTCCACTT";
const std::string fullEnd = shortEnd + "T";
}

TEST_CASE("consensus gets the bases the reads have past the end of the starting sequence", "[consensusHelper]") {
	SECTION("one base, majority of reads have it") {
		REQUIRE(fullEnd == consensusFrom(shortEnd, {{shortEnd, 3}, {fullEnd, 7}}));
	}
	SECTION("more than one base") {
		REQUIRE(fullEnd + "GC" == consensusFrom(shortEnd, {{shortEnd, 3}, {fullEnd + "GC", 7}}));
	}
	SECTION("only the bases with a majority are added") {
		REQUIRE(fullEnd == consensusFrom(shortEnd, {{shortEnd, 3}, {fullEnd, 4}, {fullEnd + "G", 3}}));
	}
	SECTION("two reads that both have it") {
		REQUIRE(fullEnd == consensusFrom(shortEnd, {{fullEnd, 2}}));
	}
}

TEST_CASE("consensus doesn't get bases past the end that most reads don't have", "[consensusHelper]") {
	SECTION("a minority of reads") {
		REQUIRE(shortEnd == consensusFrom(shortEnd, {{shortEnd, 7}, {fullEnd, 3}}));
	}
	SECTION("exactly half of the reads") {
		REQUIRE(shortEnd == consensusFrom(shortEnd, {{shortEnd, 5}, {fullEnd, 5}}));
		REQUIRE(shortEnd == consensusFrom(shortEnd, {{shortEnd, 1}, {fullEnd, 1}}));
	}
	SECTION("a majority of reads that don't agree on the base") {
		REQUIRE(shortEnd == consensusFrom(shortEnd, {{shortEnd, 4}, {shortEnd + "T", 3}, {shortEnd + "G", 3}}));
	}
}

TEST_CASE("consensus loses an end base that most reads don't have and keeps one most do", "[consensusHelper]") {
	REQUIRE(shortEnd == consensusFrom(fullEnd, {{shortEnd, 7}, {fullEnd, 3}}));
	REQUIRE(fullEnd == consensusFrom(fullEnd, {{shortEnd, 3}, {fullEnd, 7}}));
}
