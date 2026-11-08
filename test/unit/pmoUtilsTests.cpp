#include <catch_amalgamated.hpp>
#include "njhseq/objects/BioDataObject/pmo/pmoUtils.hpp"

using njhseq::VecStr;
using njhseq::table;
using namespace njhseq::pmo;

namespace {
PmoGenerationMethod mergeMethod() {
	PmoGenerationMethod ret;
	ret.program_name_ = "tester";
	ret.program_version_ = "1";
	return ret;
}

//a one target PMO, calls are sample, sequence, reads
PortableMicrohaplotypeObject onePmo(const std::string & targetName, const std::vector<std::tuple<std::string, std::string, uint32_t>> & calls, const VecStr & extraSamples = VecStr{}) {
	nlohmann::json j;
	j["pmo_header"]["pmo_version"] = "1.1.0";
	j["panel_info"] = nlohmann::json::array({{{"panel_name", "panelA"}, {"reactions", nlohmann::json::array({{{"reaction_name", "full"}, {"panel_targets", {0}}}})}}});
	j["target_info"] = nlohmann::json::array({{{"target_name", targetName}, {"forward_primer", {{"seq", "ACGT"}}}, {"reverse_primer", {{"seq", "TTGG"}}}}});
	std::vector<std::string> seqs;
	std::vector<std::string> samples;
	for (const auto & call : calls) {
		if (std::find(samples.begin(), samples.end(), std::get<0>(call)) == samples.end()) {
			samples.emplace_back(std::get<0>(call));
		}
		if (std::find(seqs.begin(), seqs.end(), std::get<1>(call)) == seqs.end()) {
			seqs.emplace_back(std::get<1>(call));
		}
	}
	nlohmann::json libraries = nlohmann::json::array();
	nlohmann::json specimens = nlohmann::json::array();
	auto allSamples = samples;
	allSamples.insert(allSamples.end(), extraSamples.begin(), extraSamples.end());
	for (size_t pos = 0; pos < allSamples.size(); ++pos) {
		specimens.push_back({{"specimen_name", allSamples[pos]}});
		libraries.push_back({{"library_sample_name", allSamples[pos]}, {"specimen_id", pos}, {"panel_id", 0}});
	}
	j["specimen_info"] = specimens;
	j["library_sample_info"] = libraries;
	nlohmann::json mhaps = nlohmann::json::array();
	for (size_t pos = 0; pos < seqs.size(); ++pos) {
		mhaps.push_back({{"seq", seqs[pos]}, {"microhaplotype_name", targetName + "." + std::to_string(pos)}});
	}
	j["representative_microhaplotypes"]["targets"] = nlohmann::json::array({{{"target_id", 0}, {"microhaplotypes", mhaps}}});
	nlohmann::json sampleResults = nlohmann::json::array();
	for (size_t pos = 0; pos < samples.size(); ++pos) {
		nlohmann::json calledMhaps = nlohmann::json::array();
		for (const auto & call : calls) {
			if (std::get<0>(call) == samples[pos]) {
				auto seqPos = std::distance(seqs.begin(), std::find(seqs.begin(), seqs.end(), std::get<1>(call)));
				calledMhaps.push_back({{"mhap_id", seqPos}, {"reads", std::get<2>(call)}});
			}
		}
		sampleResults.push_back({{"library_sample_id", pos}, {"target_results", nlohmann::json::array({{{"mhaps_target_id", 0}, {"mhaps", calledMhaps}}})}});
	}
	j["detected_microhaplotypes"] = nlohmann::json::array({{{"library_samples", sampleResults}}});
	auto ret = PortableMicrohaplotypeObject::from_json(j);
	ret.validate();
	return ret;
}

//reads per sample, target and sequence by following the ids in a PMO
std::map<std::tuple<std::string, std::string, std::string>, uint32_t> callsOf(const PortableMicrohaplotypeObject & pmoObj) {
	std::map<std::tuple<std::string, std::string, std::string>, uint32_t> ret;
	const auto & reps = pmoObj.representative_microhaplotypes_.targets_;
	for (const auto & run : pmoObj.detected_microhaplotypes_) {
		for (const auto & sample : run.library_samples_) {
			for (const auto & results : sample.target_results_) {
				const auto & rep = reps.at(results.mhaps_target_id_);
				for (const auto & mhap : results.mhaps_) {
					ret[{pmoObj.library_sample_info_.at(sample.library_sample_id_).library_sample_name_,
							pmoObj.target_info_.at(rep.target_id_).target_name_,
							rep.microhaplotypes_.at(mhap.mhap_id_).seq_}] = mhap.reads_;
				}
			}
		}
	}
	return ret;
}
}  // namespace

TEST_CASE("mergeMinimumPmos merges per target PMOs over the same samples", "[pmoUtils]") {
	auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}, {"s2", "AAAA", 5}, {"s2", "CCCC", 7}});
	auto pmo2 = onePmo("tar2", {{"s2", "GGGG", 3}, {"s3", "TTTT", 4}}, {"s1"});
	auto merged = mergeMinimumPmos({pmo1, pmo2}, mergeMethod());

	REQUIRE(2 == merged.target_info_.size());
	REQUIRE(1 == merged.panel_info_.size());
	REQUIRE(std::vector<uint32_t>{0, 1} == merged.panel_info_[0].reactions_[0].panel_targets_);
	REQUIRE(3 == merged.library_sample_info_.size());
	REQUIRE(3 == merged.specimen_info_.size());
	REQUIRE(1 == merged.detected_microhaplotypes_.size());
	REQUIRE(PMO_SCHEMA_VERSION == merged.pmo_header_.pmo_version_);
	REQUIRE("tester" == merged.pmo_header_.generation_method_.value().program_name_);

	auto expected = callsOf(pmo1);
	auto second = callsOf(pmo2);
	expected.insert(second.begin(), second.end());
	REQUIRE(expected == callsOf(merged));
}

TEST_CASE("mergeMinimumPmos combines the same target over different samples by sequence", "[pmoUtils]") {
	auto merged = mergeMinimumPmos({onePmo("tar1", {{"s1", "AAAA", 10}}), onePmo("tar1", {{"s2", "AAAA", 5}, {"s2", "CCCC", 7}})}, mergeMethod());
	REQUIRE(1 == merged.target_info_.size());
	REQUIRE(1 == merged.representative_microhaplotypes_.targets_.size());
	const auto & mhaps = merged.representative_microhaplotypes_.targets_[0].microhaplotypes_;
	REQUIRE(2 == mhaps.size());
	//tar1.0 was already used for AAAA and tar1.1 is free
	REQUIRE("tar1.0" == mhaps[0].microhaplotype_name_.value());
	REQUIRE("tar1.1" == mhaps[1].microhaplotype_name_.value());
	REQUIRE(3 == callsOf(merged).size());
}

TEST_CASE("mergeMinimumPmos drops haplotype names reused for a different sequence", "[pmoUtils]") {
	auto merged = mergeMinimumPmos({onePmo("tar1", {{"s1", "AAAA", 10}}), onePmo("tar1", {{"s2", "CCCC", 5}})}, mergeMethod());
	const auto & mhaps = merged.representative_microhaplotypes_.targets_[0].microhaplotypes_;
	REQUIRE("tar1.0" == mhaps[0].microhaplotype_name_.value());
	REQUIRE_FALSE(mhaps[1].microhaplotype_name_.has_value());
}

TEST_CASE("mergeMinimumPmos errors on conflicts", "[pmoUtils]") {
	SECTION("the same results twice") {
		auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo1}, mergeMethod()), Catch::Matchers::ContainsSubstring("more than once"));
	}
	SECTION("different primers") {
		auto pmo2 = onePmo("tar1", {{"s2", "AAAA", 10}});
		pmo2.target_info_[0].forward_primer_.seq_ = "GGGG";
		REQUIRE_THROWS_WITH(mergeMinimumPmos({onePmo("tar1", {{"s1", "AAAA", 10}}), pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("different primers"));
	}
	SECTION("different meta for a specimen") {
		auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
		auto pmo2 = onePmo("tar2", {{"s1", "CCCC", 10}});
		pmo1.specimen_info_[0].collection_country_ = "Uganda";
		pmo2.specimen_info_[0].collection_country_ = "Kenya";
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("collection_country"));
	}
	SECTION("sections that can't be merged yet") {
		auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
		pmo1.project_info_ = std::vector<ProjectInfo>{};
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1}, mergeMethod()), Catch::Matchers::ContainsSubstring("only minimum PMOs"));
	}
}

TEST_CASE("MinimumPmoMerger gives the same results added one at a time", "[pmoUtils]") {
	auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
	auto pmo2 = onePmo("tar2", {{"s1", "CCCC", 5}});
	MinimumPmoMerger merger(mergeMethod());
	merger.add(pmo1);
	merger.add(pmo2);
	REQUIRE(merger.finish().to_json() == mergeMinimumPmos({pmo1, pmo2}, mergeMethod()).to_json());
}

TEST_CASE("addSpecimenMeta adds schema and custom fields", "[pmoUtils]") {
	auto pmoObj = onePmo("tar1", {{"s1", "AAAA", 10}, {"s2", "AAAA", 5}});
	table meta(VecStr{"sample", "collection_country", "host_taxon_id", "site"});
	meta.content_.emplace_back(VecStr{"s1", "Uganda", "9606", "A"});
	meta.content_.emplace_back(VecStr{"s2", "NA", "", "B"});
	meta.content_.emplace_back(VecStr{"s3", "Kenya", "9606", "C"});
	auto results = addSpecimenMeta(pmoObj, meta, "sample");

	REQUIRE(2 == results.specimensUpdated_);
	REQUIRE(VecStr{"s3"} == results.metaSamplesNotInPmo_);
	REQUIRE("Uganda" == pmoObj.specimen_info_[0].collection_country_.value());
	REQUIRE(9606 == pmoObj.specimen_info_[0].host_taxon_id_.value());
	REQUIRE("A" == pmoObj.specimen_info_[0].extras_.at("site").get<std::string>());
	REQUIRE_FALSE(pmoObj.specimen_info_[1].collection_country_.has_value());
	REQUIRE("B" == pmoObj.specimen_info_[1].extras_.at("site").get<std::string>());
}

TEST_CASE("addSpecimenMeta errors on bad values and duplicate samples", "[pmoUtils]") {
	auto pmoObj = onePmo("tar1", {{"s1", "AAAA", 10}});
	SECTION("a number field that isn't a number") {
		table meta(VecStr{"sample", "host_taxon_id"});
		meta.content_.emplace_back(VecStr{"s1", "human"});
		REQUIRE_THROWS_WITH(addSpecimenMeta(pmoObj, meta, "sample"), Catch::Matchers::ContainsSubstring("host_taxon_id"));
	}
	SECTION("duplicate samples") {
		table meta(VecStr{"sample", "site"});
		meta.content_.emplace_back(VecStr{"s1", "A"});
		meta.content_.emplace_back(VecStr{"s1", "B"});
		REQUIRE_THROWS_WITH(addSpecimenMeta(pmoObj, meta, "sample"), Catch::Matchers::ContainsSubstring("more than once"));
	}
}
