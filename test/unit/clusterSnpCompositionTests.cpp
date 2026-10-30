#include <catch_amalgamated.hpp>
#include "njhseq/objects/seqObjects/Clusters/cluster.hpp"

using njhseq::cluster;

TEST_CASE("snp strand composition is not checked when the cluster is one direction", "[cluster][snpBreakout]") {
	SECTION("all forward reads") {
		REQUIRE(cluster::snpPassesStrandComposition(10, 0, 0, 50, 0.01));
	}
	SECTION("all reverse complemented reads") {
		REQUIRE(cluster::snpPassesStrandComposition(10, 10, 50, 50, 0.01));
	}
}

TEST_CASE("snp strand composition is tested against the cluster's composition", "[cluster][snpBreakout]") {
	SECTION("balanced cluster, snp in both directions passes") {
		REQUIRE(cluster::snpPassesStrandComposition(10, 5, 25, 50, 0.01));
		//2 of 10 is two sided p ~0.11
		REQUIRE(cluster::snpPassesStrandComposition(10, 2, 25, 50, 0.01));
	}
	SECTION("balanced cluster, snp in only one direction fails") {
		//all 10 forward or all 10 reverse is two sided p ~0.002
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 0, 25, 50, 0.01));
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 10, 25, 50, 0.01));
		//a larger snp mostly in one direction, 30 of 100
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(100, 30, 25, 50, 0.01));
	}
	SECTION("a small snp in one direction isn't enough evidence of bias") {
		//all 4 forward is two sided p 0.125
		REQUIRE(cluster::snpPassesStrandComposition(4, 0, 25, 50, 0.01));
	}
	SECTION("a mostly forward cluster expects mostly forward snps") {
		//cluster is 10% reverse complemented, a snp with no reverse reads is expected
		REQUIRE(cluster::snpPassesStrandComposition(10, 0, 5, 50, 0.01));
		//but one with only reverse reads is not
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(10, 10, 5, 50, 0.01));
	}
	SECTION("a snp with no reads fails") {
		REQUIRE_FALSE(cluster::snpPassesStrandComposition(0, 0, 25, 50, 0.01));
	}
}

TEST_CASE("reverse complemented reads are found by the extractor's _Comp with an optional count or fraction after it", "[cluster][snpBreakout]") {
	REQUIRE(cluster::readIsRevComp("M01:1:2_Comp"));
	//qluster names collapsed identical reads <first read>_t<count>, kluster can give fractional counts
	REQUIRE(cluster::readIsRevComp("M01:1:2_Comp_t5"));
	REQUIRE(cluster::readIsRevComp("M01:1:2_Comp_t2.5"));
	REQUIRE(cluster::readIsRevComp("M01:1:2_Comp_f0.25"));
	REQUIRE(cluster::readIsRevComp("M01:1:2_Comp_t1e+06"));
	REQUIRE_FALSE(cluster::readIsRevComp("M01:1:2_t5"));
	REQUIRE_FALSE(cluster::readIsRevComp("M01:1:2"));
	//_Comp elsewhere in the name, e.g. a sample name
	REQUIRE_FALSE(cluster::readIsRevComp("Pt_Comp2.M01:1:2_t5"));
	REQUIRE_FALSE(cluster::readIsRevComp("M01:1:2[sample=Pt_Comp]_t5"));
	//only one of the known appendages is skipped
	REQUIRE_FALSE(cluster::readIsRevComp("M01:1:2_Comp_t"));
	REQUIRE_FALSE(cluster::readIsRevComp("M01:1:2_Comp_x5"));
	REQUIRE_FALSE(cluster::readIsRevComp("M01:1:2_Comp_t5_t3"));
}

TEST_CASE("snp error model background rate is a quantile counting unobserved mismatches as 0", "[cluster][snpBreakout]") {
	SECTION("mostly unobserved mismatches give the minimum rate") {
		//300 possible mismatches, 5 observed, the 95th percentile is still a 0
		REQUIRE(cluster::snpBackgroundErrorRate({0.01, 0.02, 0.03, 0.2, 0.3}, 300, 0.95, 0.001) == Catch::Approx(0.001));
	}
	SECTION("enough observed mismatches give an observed rate") {
		//10 possible, 9 observed, the 0.5 quantile (position 4 of 0-9) is the 4th observed rate since 1 is a 0
		std::vector<double> rates{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
		REQUIRE(cluster::snpBackgroundErrorRate(rates, 10, 0.5, 0.001) == Catch::Approx(0.04));
	}
	SECTION("no observed mismatches give the minimum rate") {
		REQUIRE(cluster::snpBackgroundErrorRate({}, 300, 0.95, 0.002) == Catch::Approx(0.002));
	}
}

TEST_CASE("snp error model upper tail binomial test", "[cluster][snpBreakout]") {
	//1000 reads with a 1% error rate expect 10 reads with the error
	REQUIRE_FALSE(cluster::snpAboveBackgroundError(12, 1000, 0.01, 0.01));
	REQUIRE(cluster::snpAboveBackgroundError(40, 1000, 0.01, 0.01));
	REQUIRE_FALSE(cluster::snpAboveBackgroundError(0, 1000, 0.01, 0.01));
	REQUIRE_FALSE(cluster::snpAboveBackgroundError(5, 0, 0.01, 0.01));
}

TEST_CASE("snp breakout read count cut offs with and without the error model", "[cluster][snpBreakout]") {
	cluster::snpBreakoutPars pars;
	pars.hardCutOff = 10;
	pars.hardSnpFreqCutOff = 0.01;
	SECTION("count cut offs only") {
		REQUIRE(cluster::snpPassesBreakoutCutOffs(12, 1000, pars, 0.01, 0.01));
		REQUIRE_FALSE(cluster::snpPassesBreakoutCutOffs(8, 1000, pars, 0.001, 0.01));
	}
	SECTION("error model instead of the count cut offs") {
		pars.useErrorModel = true;
		//passes the counts but not above a 1% background
		REQUIRE_FALSE(cluster::snpPassesBreakoutCutOffs(12, 1000, pars, 0.01, 0.01));
		//below the count cut off but well above a 0.1% background
		REQUIRE(cluster::snpPassesBreakoutCutOffs(8, 1000, pars, 0.001, 0.01));
		REQUIRE(cluster::snpPassesBreakoutCutOffs(40, 1000, pars, 0.01, 0.01));
	}
	SECTION("error model as an extra filter needs both") {
		pars.useErrorModel = true;
		pars.errorModelAsFilter = true;
		REQUIRE_FALSE(cluster::snpPassesBreakoutCutOffs(12, 1000, pars, 0.01, 0.01));
		REQUIRE_FALSE(cluster::snpPassesBreakoutCutOffs(8, 1000, pars, 0.001, 0.01));
		REQUIRE(cluster::snpPassesBreakoutCutOffs(40, 1000, pars, 0.01, 0.01));
	}
}
