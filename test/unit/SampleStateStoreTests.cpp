#include <catch_amalgamated.hpp>
#include "njhseq/objects/collapseObjects/SampleCollapseCollectionUtils/SampleStateStore.hpp"

#include <thread>

using njhseq::collapse::SampleState;
using njhseq::collapse::SampleStateStore;

namespace {
SampleState makeState(const std::string & finalText) {
	SampleState state;
	state.final_ = finalText;
	state.finalInputClusters_ = finalText + finalText;
	return state;
}
}  // namespace

TEST_CASE("SampleState size counts only the set texts", "[SampleStateStore]") {
	SampleState state;
	REQUIRE(0 == state.sizeInBytes());
	state.final_ = "abc";
	state.refCompInfos_ = "de";
	REQUIRE(5 == state.sizeInBytes());
}

TEST_CASE("SampleStateStore set, get and erase", "[SampleStateStore]") {
	SampleStateStore store;
	REQUIRE_FALSE(store.has("s1"));
	REQUIRE_FALSE(store.get("s1").has_value());

	store.set("s1", makeState("@r1\nACGT\n+\nIIII\n"));
	REQUIRE(store.has("s1"));
	auto state = store.get("s1");
	REQUIRE(state.has_value());
	REQUIRE(*state->final_ == "@r1\nACGT\n+\nIIII\n");
	REQUIRE_FALSE(state->excluded_.has_value());

	SECTION("setting again replaces the whole state") {
		SampleState replacement;
		replacement.excluded_ = "@r2\nTT\n+\nII\n";
		store.set("s1", replacement);
		auto replaced = store.get("s1");
		REQUIRE_FALSE(replaced->final_.has_value());
		REQUIRE(*replaced->excluded_ == "@r2\nTT\n+\nII\n");
		//a state already handed out is a copy and isn't changed
		REQUIRE(*state->final_ == "@r1\nACGT\n+\nIIII\n");
	}
	SECTION("erase removes the sample") {
		store.erase("s1");
		REQUIRE_FALSE(store.has("s1"));
		REQUIRE(store.samples().empty());
	}
}

TEST_CASE("SampleStateStore samples are sorted", "[SampleStateStore]") {
	SampleStateStore store;
	store.set("s2", makeState("ab"));
	store.set("s1", makeState("c"));
	REQUIRE(store.samples() == njhseq::VecStr{"s1", "s2"});
}

TEST_CASE("SampleStateStore round trips every text and compresses", "[SampleStateStore]") {
	SampleStateStore store;
	SampleState state;
	std::string fastq;
	for (uint32_t pos = 0; pos < 500; ++pos) {
		fastq += "@seq." + std::to_string(pos) + "[sample=s1]_t3\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
	}
	state.final_ = fastq;
	state.finalInputClusters_ = "";
	state.excludedInputClusters_ = fastq + "extra";
	state.refCompInfosExcluded_ = "read\tbestExpected\nseq.0\tref1\n";
	store.set("s1", state);
	auto back = store.get("s1");
	REQUIRE(back.has_value());
	REQUIRE(back->final_ == state.final_);
	REQUIRE(back->finalInputClusters_ == state.finalInputClusters_);
	REQUIRE_FALSE(back->excluded_.has_value());
	REQUIRE(back->excludedInputClusters_ == state.excludedInputClusters_);
	REQUIRE_FALSE(back->refCompInfos_.has_value());
	REQUIRE(back->refCompInfosExcluded_ == state.refCompInfosExcluded_);
	REQUIRE(store.uncompressedSizeInBytes() > state.sizeInBytes());
	REQUIRE(store.sizeInBytes() < state.sizeInBytes() / 4);
}

TEST_CASE("SampleStateStore can be set from many threads", "[SampleStateStore][threads]") {
	SampleStateStore store;
	//Catch2 assertions aren't thread safe so record results and check after joining
	std::vector<uint32_t> foundAfterSet(8, 0);
	std::vector<std::thread> threads;
	for (uint32_t t = 0; t < 8; ++t) {
		threads.emplace_back([&store, &foundAfterSet, t]() {
			for (uint32_t samp = 0; samp < 200; ++samp) {
				auto name = "t" + std::to_string(t) + "_s" + std::to_string(samp);
				store.set(name, makeState(name));
				if (store.has(name)) {
					++foundAfterSet[t];
				}
			}
		});
	}
	for (auto & th : threads) {
		th.join();
	}
	for (const auto found : foundAfterSet) {
		REQUIRE(200 == found);
	}
	REQUIRE(store.samples().size() == 8 * 200);
	REQUIRE(*store.get("t3_s17")->final_ == "t3_s17");
}
