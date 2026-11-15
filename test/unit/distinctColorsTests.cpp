#include <catch_amalgamated.hpp>
#include "njhseq/utils/kMeans.hpp"
#include "njhseq/utils/distinctColors.hpp"

using njhseq::KMeans;

TEST_CASE("KMeans dist2 and mean", "[kMeans]") {
	REQUIRE(KMeans::dist2({0, 0}, {3, 4}) == 25);
	REQUIRE(KMeans::mean({{0, 0}, {2, 4}}) == KMeans::Point{1, 2});
	REQUIRE(KMeans::mean({}).empty());
}

TEST_CASE("KMeans farthestPoints", "[kMeans]") {
	const std::vector<KMeans::Point> points{{0, 0}, {1, 0}, {10, 0}, {5, 0}};
	//from the mean (4, 0) the farthest is 10, then 0, then 5
	REQUIRE(KMeans::farthestPoints(points, 3, KMeans::mean(points)) == std::vector<uint32_t>{2, 0, 3});
	//all of them orders them
	REQUIRE(KMeans::farthestPoints(points, 4, {0, 0}) == std::vector<uint32_t>{2, 0, 3, 1});
	//ties go to the lowest index
	REQUIRE(KMeans::farthestPoints({{1, 0}, {-1, 0}}, 1, {0, 0}) == std::vector<uint32_t>{0});
	REQUIRE_THROWS(KMeans::farthestPoints(points, 5, {0, 0}));
}

TEST_CASE("KMeans run", "[kMeans]") {
	//two clumps, starting both centers in the first clump
	const std::vector<KMeans::Point> points{{0, 0}, {0, 1}, {1, 0}, {10, 10}, {10, 11}, {11, 10}};
	SECTION("free centers end at the clump means") {
		//starting both in the first clump still finds the two clumps
		auto res = KMeans::run(points, {0, 1});
		REQUIRE(res.converged_);
		REQUIRE(res.assignments_ == std::vector<uint32_t>{0, 0, 0, 1, 1, 1});
		res = KMeans::run(points, {0, 3});
		REQUIRE(res.converged_);
		REQUIRE(res.assignments_ == std::vector<uint32_t>{0, 0, 0, 1, 1, 1});
		REQUIRE(res.centers_[0] == KMeans::Point{1.0 / 3, 1.0 / 3});
		REQUIRE(res.centers_[1] == KMeans::Point{31.0 / 3, 31.0 / 3});
	}
	SECTION("snapped centers are distinct points") {
		KMeans::Pars pars;
		pars.snapToPoints_ = true;
		auto res = KMeans::run(points, {0, 3}, pars);
		REQUIRE(res.converged_);
		REQUIRE(res.centerPoints_ == std::vector<uint32_t>{0, 3});
		//both centers on the same clump can't snap to the same point
		res = KMeans::run({{0, 0}, {0, 0}, {5, 5}}, {0, 1}, pars);
		REQUIRE(res.centerPoints_[0] != res.centerPoints_[1]);
	}
	SECTION("max iterations") {
		KMeans::Pars pars;
		pars.maxIterations_ = 1;
		auto res = KMeans::run(points, {0, 1}, pars);
		REQUIRE(1 == res.iterations_);
		REQUIRE_FALSE(res.converged_);
	}
	REQUIRE(KMeans::run({}, {}).assignments_.empty());
}

TEST_CASE("genDistinctColorHexStrs", "[distinctColors]") {
	REQUIRE(njhseq::genDistinctColorHexStrs(0).empty());
	SECTION("distinct, valid and the same every time") {
		for (const uint32_t count : {1u, 7u, 34u, 120u}) {
			const auto colors = njhseq::genDistinctColorHexStrs(count);
			REQUIRE(colors.size() == count);
			REQUIRE(std::set<std::string>(colors.begin(), colors.end()).size() == count);
			for (const auto & color : colors) {
				REQUIRE(7 == color.size());
				REQUIRE('#' == color.front());
			}
			REQUIRE(colors == njhseq::genDistinctColorHexStrs(count));
		}
	}
	SECTION("the same colors as seqServer utils.js distinctColorHexStrs()") {
		const std::vector<std::pair<uint32_t, std::string>> expected{
		{1, "#b67262"},
		{2, "#9ba454,#ba689e"},
		{3, "#7bab5d,#a76fb8,#c86949"},
		{5, "#5faf66,#a958ab,#ce6552,#6791d7,#b09e40"},
		{10, "#6e6ad5,#c6bc4f,#cf4658,#4fbbbc,#64bf69,#c46d35,#b84dab,#717d39,#d17492,#7b7ec0"},
		{21, "#5cd061,#873ca9,#d85f2d,#48a8d3,#e0819f,#c9bb36,#5c722d,#43b08b,#d53e58,#6b62a5,#bf8535,#d84cb0,#9b442f,#6c9f30,#de886a,#a4a5eb,#af437c,#626dd5,#bc74d8,#c2a962,#7eba69"},
		{34, "#47c25e,#9433a1,#e2712e,#58b5e1,#e793ac,#c9bb36,#1d866d,#6e6e23,#d53e58,#707ae4,#b53c7d,#b4a1e2,#9b442f,#49cbca,#835a9d,#d1b66e,#96c676,#387836,#4a78bb,#8fca40,#d380d4,#e5836a,#f0ac3a,#5bbd8f,#9f6c3b,#a0526c,#c4362b,#d84cb0,#6c9f30,#b88835,#ac64db,#5e4db9,#829951,#e66493"},
		{35, "#47c25e,#9433a1,#e2712e,#58b5e1,#606222,#e793ac,#c9bb36,#1d866d,#d53e58,#707ae4,#d7b478,#b53c7d,#b4a1e2,#9b442f,#49cbca,#835a9d,#5e9222,#4a78bb,#d380d4,#e5836a,#f0ac3a,#8fca40,#80995a,#9f6c3b,#a0526c,#52be98,#8bc876,#387836,#c4362b,#d84cb0,#b88835,#ac64db,#5e4db9,#e66493,#adb14c"},
		{50, "#5cd061,#972ca9,#e2712e,#58b5e1,#606222,#e793ac,#c9bb36,#65558f,#dc3559,#1d866d,#a72d70,#9b442f,#556ede,#cfb778,#d992e2,#49cbca,#b98829,#6c9f30,#de886a,#aaa4e1,#387836,#d84cb0,#c4362b,#a0526c,#e66493,#52be98,#8bc876,#9c6ae4,#8d3f90,#5884c9,#d16063,#7d5719,#3dd188,#e19f1a,#8978c9,#34a54f,#da2d82,#99c840,#769b5a,#ab7950,#b769a7,#b56029,#adb14c,#ecac5c,#266197,#d564d8,#5e4db9,#b02949,#757c24,#e2534a"},
		{100, "#49d261,#972ca9,#e37021,#58b5e1,#606222,#e793ac,#cabb26,#2b72e7,#e22959,#38b094,#65558f,#9b442f,#a72d70,#edaa83,#a8c280,#237b12,#aa7b1b,#b8b2f0,#ca81e6,#ca2d22,#1d7bc3,#f8a82c,#df43b0,#fd6ca0,#a1c54d,#a0526c,#36cbd3,#b769a7,#d1b66e,#33783f,#f982d1,#58bf7c,#5e4db9,#ab7950,#d16063,#e3a3e7,#e78361,#b02949,#8978c9,#8d3f90,#7d5719,#da2d82,#2f5bb1,#48a421,#35846e,#9c6ae4,#ecac5c,#799cee,#adb14c,#a6541d,#3e6914,#c75235,#d88884,#6c9f30,#83cc3f,#d564d8,#8f8546,#8a4777,#e1b338,#3d8bb7,#52dcbc,#8bc876,#80995a,#34a54f,#757c24,#de8a2c,#b29e27,#9a594a,#cc7e41,#7d76e4,#b32831,#266197,#e84c4b,#71a1dc,#a855bc,#569244,#6274b2,#ca6285,#b590cb,#d7529f,#ee663a,#6b4fa0,#7c3eba,#71caa5,#429369,#3d683c,#d34069,#9c4050,#8d69a2,#3dd188,#2dadb8,#4f90f2,#97b723,#bb2c95,#565bd9,#68ce6b,#a93713,#c7799a,#56ac9d,#f49446"},
		};
		for (const auto & exp : expected) {
			INFO(exp.first << " colors");
			REQUIRE(njh::conToStr(njhseq::genDistinctColorHexStrs(exp.first), ",") == exp.second);
		}
	}
	SECTION("more colors than candidates repeat") {
		njhseq::DistinctColorsPars pars;
		pars.lightnessMin_ = 50;
		pars.lightnessMax_ = 50;
		pars.gridStep_ = 50;
		const auto colors = njhseq::genDistinctColorHexStrs(10, pars);
		REQUIRE(colors.size() == 10);
		REQUIRE(std::set<std::string>(colors.begin(), colors.end()).size() < 10);
	}
	SECTION("as colors") {
		const auto colors = njhseq::genDistinctColors(3);
		const auto hexes = njhseq::genDistinctColorHexStrs(3);
		for (const auto pos : {0, 1, 2}) {
			REQUIRE(colors[pos].getHexStr() == hexes[pos]);
		}
	}
}
