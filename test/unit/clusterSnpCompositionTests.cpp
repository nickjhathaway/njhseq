#include <catch_amalgamated.hpp>
#include "njhseq/objects/seqObjects/Clusters/cluster.hpp"

using njhseq::cluster;

TEST_CASE("snp strand composition is not checked when the cluster is one direction", "[cluster][snpBreakout]") {
	SECTION("all forward reads") {
		REQUIRE(cluster::snpPassesStrandComposition(10, 0, 0, 50, 0.2));
	}
	SECTION("all reverse complemented reads") {
		REQUIRE(cluster::snpPassesStrandComposition(10, 10, 50, 50, 0.2));
	}
}

TEST_CASE("snp strand composition uses the snp's own reads when the cluster is mixed", "[cluster][snpBreakout]") {
	SECTION("snp found in both directions passes") {
		REQUIRE(cluster::snpPassesStrandComposition(10, 5, 25, 50, 0.2));
		REQUIRE(cluster::snpPassesStrandComposition(10, 3, 25, 50, 0.2));
		REQUIRE(cluster::snpPassesStrandComposition(10, 7, 25, 50, 0.2));
	}
	SECTION("snp found in only one direction fails") {
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 0, 25, 50, 0.2));
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 10, 25, 50, 0.2));
	}
	SECTION("snp mostly in one direction fails at the cut off") {
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 2, 25, 50, 0.2));
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 8, 25, 50, 0.2));
	}
	SECTION("the cluster wide composition doesn't matter, only the snp's") {
		//cluster is 2% reverse complemented but the snp is balanced
		REQUIRE(cluster::snpPassesStrandComposition(10, 5, 1, 50, 0.2));
		//cluster is balanced but the snp is only forward
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 0, 25, 50, 0.2));
	}
	SECTION("a snp with no reads fails") {
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(0, 0, 25, 50, 0.2));
	}
}
