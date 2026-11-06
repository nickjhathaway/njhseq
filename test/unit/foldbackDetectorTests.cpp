#include <catch_amalgamated.hpp>
#include "njhseq/seqToolsUtils/FoldbackDetector.hpp"
#include <random>

using namespace njhseq;

namespace {

std::string randomSeq(const uint32_t length, std::mt19937 & gen) {
	std::uniform_int_distribution<int> base(0, 3);
	std::string ret(length, 'A');
	for (auto & letter : ret) {
		letter = "ACGT"[base(gen)];
	}
	return ret;
}

std::string revComp(const std::string & seq) {
	std::string ret(seq.rbegin(), seq.rend());
	for (auto & letter : ret) {
		switch (letter) {
			case 'A': letter = 'T'; break;
			case 'C': letter = 'G'; break;
			case 'G': letter = 'C'; break;
			case 'T': letter = 'A'; break;
			default: break;
		}
	}
	return ret;
}

/** substitutions, insertions and deletions at the given rate */
std::string addErrors(const std::string & seq, const double rate, std::mt19937 & gen) {
	std::uniform_real_distribution<double> roll(0, 1);
	std::uniform_int_distribution<int> base(0, 3);
	std::string ret;
	for (const auto letter : seq) {
		const auto chance = roll(gen);
		if (chance < rate / 2) {
			ret.push_back("ACGT"[base(gen)]);
		} else if (chance < rate * 3 / 4) {
			ret.push_back(letter);
			ret.push_back("ACGT"[base(gen)]);
		} else if (chance >= rate) {
			ret.push_back(letter);
		}
	}
	return ret;
}

}  // namespace

TEST_CASE("a read with no fold isn't called", "[FoldbackDetector]") {
	std::mt19937 gen(1);
	const FoldbackDetector detector;
	for (uint32_t read = 0; read < 50; ++read) {
		CHECK_FALSE(detector.detect(randomSeq(2000, gen)).called_);
	}
}

TEST_CASE("reads too short to have a fold are handled", "[FoldbackDetector]") {
	const FoldbackDetector detector;
	CHECK_FALSE(detector.detect("").hasCandidate());
	CHECK_FALSE(detector.detect("ACGTACGT").hasCandidate());
}

TEST_CASE("a perfect fold gives the exact arms and loop", "[FoldbackDetector]") {
	std::mt19937 gen(2);
	const FoldbackDetector detector;
	const auto arm = randomSeq(300, gen);
	//A can't pair with A so neither the loop's ends nor the read's ends can extend the arms
	const auto loop = "A" + randomSeq(100, gen) + "A";
	const auto fold = detector.detect("C" + arm + loop + revComp(arm) + "C");
	REQUIRE(fold.called_);
	CHECK(1 == fold.leftStart_);
	CHECK(301 == fold.leftEnd_);
	CHECK(loop.size() == fold.loopLen());
	CHECK(fold.rightStart_ == fold.leftEnd_ + loop.size());
	CHECK(fold.rightEnd_ == fold.rightStart_ + 300);
	CHECK(1.0 == fold.identity_);
}

TEST_CASE("the loop doesn't change the call", "[FoldbackDetector]") {
	std::mt19937 gen(3);
	const FoldbackDetector detector;
	for (const auto loopLen : {0u, 1u, 10u, 50u, 200u, 700u}) {
		const auto arm = randomSeq(300, gen);
		const auto fold = detector.detect(arm + randomSeq(loopLen, gen) + revComp(arm));
		CHECK(fold.called_);
		//a few bases of the random loop can pair by chance, on both arms
		CHECK(fold.loopLen() <= loopLen);
		CHECK(fold.loopLen() + 12 >= loopLen);
	}
}

TEST_CASE("folds are found through nanopore levels of error", "[FoldbackDetector]") {
	std::mt19937 gen(4);
	const FoldbackDetector detector;
	uint32_t called = 0;
	std::vector<uint32_t> leftEndOff;
	for (uint32_t read = 0; read < 100; ++read) {
		const auto arm = randomSeq(400, gen);
		const auto left = addErrors(arm, 0.06, gen);
		const auto loop = randomSeq(50, gen);
		const auto right = addErrors(revComp(arm), 0.06, gen);
		const auto fold = detector.detect(left + loop + right);
		if (fold.called_) {
			++called;
			leftEndOff.emplace_back(fold.leftEnd_ > left.size() ? fold.leftEnd_ - left.size() : left.size() - fold.leftEnd_);
		}
	}
	CHECK(called >= 99);
	std::sort(leftEndOff.begin(), leftEndOff.end());
	CHECK(leftEndOff[leftEndOff.size() / 2] <= 3);
}

TEST_CASE("an (AT)n tract can't make a call on its own", "[FoldbackDetector]") {
	std::mt19937 gen(5);
	const FoldbackDetector detector;
	std::string tract;
	for (uint32_t repeat = 0; repeat < 60; ++repeat) {
		tract += "AT";
	}
	CHECK_FALSE(detector.detect(randomSeq(500, gen) + tract + randomSeq(500, gen)).called_);
}

TEST_CASE("arms that pair around an (AT)n tract close the loop", "[FoldbackDetector]") {
	std::mt19937 gen(6);
	const FoldbackDetector detector;
	const auto arm = randomSeq(200, gen);
	const auto fold = detector.detect(arm + std::string("ATATATATATATATATATATATAT") + revComp(arm));
	REQUIRE(fold.called_);
	CHECK(fold.loopLen() <= 2);
}

TEST_CASE("sequence shared by the far ends of a read isn't a fold", "[FoldbackDetector]") {
	std::mt19937 gen(7);
	const FoldbackDetector detector;
	const auto adapter = randomSeq(70, gen);
	const auto fold = detector.detect(adapter + randomSeq(2000, gen) + revComp(adapter));
	CHECK(fold.hasCandidate());
	CHECK_FALSE(fold.called_);
}

TEST_CASE("the case of the read doesn't matter and other letters don't vote", "[FoldbackDetector]") {
	std::mt19937 gen(8);
	const FoldbackDetector detector;
	const auto arm = randomSeq(200, gen);
	const auto read = arm + randomSeq(30, gen) + revComp(arm);
	auto lower = read;
	std::transform(lower.begin(), lower.end(), lower.begin(), [](char letter) { return static_cast<char>(std::tolower(letter)); });
	const auto upperFold = detector.detect(read);
	const auto lowerFold = detector.detect(lower);
	CHECK(upperFold.called_);
	CHECK(lowerFold.called_);
	CHECK(upperFold.leftEnd_ == lowerFold.leftEnd_);
	CHECK(upperFold.rightStart_ == lowerFold.rightStart_);
	CHECK_FALSE(detector.detect(std::string(500, 'N')).hasCandidate());
}
