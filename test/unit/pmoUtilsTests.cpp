#include <catch_amalgamated.hpp>
#include "njhseq/objects/BioDataObject/pmo/pmoUtils.hpp"
#include "njhseq/GenomeUtils/TargetGenomeAssignment.hpp"
#include <fstream>

using njhseq::VecStr;
using njhseq::table;
using namespace njhseq::pmo;
namespace bfs = njh::bfs;

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

//read counts by stage for every library sample of a PMO, one count per target, the umis go in the stage extras
ReadCountsByStage readCountsFor(const PortableMicrohaplotypeObject & pmoObj, const std::vector<uint32_t> & targetIds, uint32_t reads) {
	ReadCountsByStage ret;
	for (size_t libraryId = 0; libraryId < pmoObj.library_sample_info_.size(); ++libraryId) {
		ReadCountsByStageForLibrarySample libraryCounts;
		libraryCounts.library_sample_id_ = libraryId;
		libraryCounts.total_raw_count_ = 1000;
		libraryCounts.read_counts_for_targets_ = std::vector<ReadCountsByStageForTarget>{};
		for (const auto targetId : targetIds) {
			ReadCountsByStageForTarget targetCounts;
			targetCounts.target_id_ = targetId;
			StageReadCounts stage;
			stage.stage_ = "extracted";
			stage.reads_ = reads;
			stage.extras_["umis"] = reads / 10;
			targetCounts.stages_.emplace_back(stage);
			libraryCounts.read_counts_for_targets_->emplace_back(targetCounts);
		}
		ret.read_counts_by_library_sample_by_stage_.emplace_back(libraryCounts);
	}
	return ret;
}

void addBioinformatics(PortableMicrohaplotypeObject & pmoObj, const std::string & runName) {
	BioMethod method;
	method.program_ = "MIPWrangler";
	method.program_version_ = "1";
	BioinformaticsMethodInfo methods;
	methods.methods_.emplace_back(method);
	pmoObj.bioinformatics_methods_info_ = std::vector<BioinformaticsMethodInfo>{methods};
	BioinformaticsRunInfo run;
	run.bioinformatics_methods_id_ = 0;
	run.bioinformatics_run_name_ = runName;
	pmoObj.bioinformatics_run_info_ = std::vector<BioinformaticsRunInfo>{run};
}

//reads and umis per sample, target and stage by following the ids in a PMO
std::map<std::tuple<std::string, std::string, std::string>, std::pair<uint32_t, uint32_t>> readCountsOf(const PortableMicrohaplotypeObject & pmoObj) {
	std::map<std::tuple<std::string, std::string, std::string>, std::pair<uint32_t, uint32_t>> ret;
	for (const auto & counts : pmoObj.read_counts_by_stage_.value()) {
		for (const auto & libraryCounts : counts.read_counts_by_library_sample_by_stage_) {
			for (const auto & targetCounts : libraryCounts.read_counts_for_targets_.value()) {
				for (const auto & stage : targetCounts.stages_) {
					ret[{pmoObj.library_sample_info_.at(libraryCounts.library_sample_id_).library_sample_name_,
							pmoObj.target_info_.at(targetCounts.target_id_).target_name_, stage.stage_}] = {stage.reads_, stage.extras_.at("umis").get<uint32_t>()};
				}
			}
		}
	}
	return ret;
}

GenomeInfo genome(const std::string & name) {
	GenomeInfo ret;
	ret.name_ = name;
	ret.genome_version_ = "v1";
	ret.url_ = "NA";
	return ret;
}

GenomicLocation location(uint32_t genomeId, const std::string & chrom, uint32_t start, const std::string & alt) {
	GenomicLocation ret;
	ret.genome_id_ = genomeId;
	ret.chrom_ = chrom;
	ret.start_ = start;
	ret.end_ = start + 1;
	ret.alt_seq_ = alt;
	return ret;
}

ProteinVariant proteinVariant(uint32_t genomeId, uint32_t aaPos, const std::string & aa) {
	ProteinVariant ret;
	ret.protein_location_ = location(genomeId, "transcript1", aaPos, aa);
	ret.codon_genomic_location_ = location(genomeId, "chr1", 100 + aaPos * 3, "AAA");
	return ret;
}
}  // namespace

TEST_CASE("mergeMinimumPmos carries over genomes and haplotype variants", "[pmoUtils]") {
	auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
	pmo1.targeted_genomes_ = std::vector<GenomeInfo>{genome("g1")};
	pmo1.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].associated_seq_variants_ = std::vector<GenomicLocation>{location(0, "chr1", 5, "T")};
	pmo1.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].associated_protein_variants_ = std::vector<ProteinVariant>{proteinVariant(0, 7, "K")};
	//genome g1 is the second genome in this one
	auto pmo2 = onePmo("tar2", {{"s1", "CCCC", 5}});
	pmo2.targeted_genomes_ = std::vector<GenomeInfo>{genome("g2"), genome("g1")};
	pmo2.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].associated_seq_variants_ = std::vector<GenomicLocation>{location(1, "chr1", 9, "G"), location(0, "chr9", 3, "A")};
	pmo2.target_info_[0].insert_location_ = location(1, "chr1", 0, "A");
	auto merged = mergeMinimumPmos({pmo1, pmo2}, mergeMethod());

	REQUIRE(2 == merged.targeted_genomes_.value().size());
	REQUIRE("g1" == merged.targeted_genomes_.value()[0].name_);
	REQUIRE("g2" == merged.targeted_genomes_.value()[1].name_);
	const auto & reps = merged.representative_microhaplotypes_.targets_;
	REQUIRE(0 == reps[0].microhaplotypes_[0].associated_seq_variants_.value()[0].genome_id_);
	REQUIRE(0 == reps[0].microhaplotypes_[0].associated_protein_variants_.value()[0].protein_location_.genome_id_);
	REQUIRE(0 == reps[0].microhaplotypes_[0].associated_protein_variants_.value()[0].codon_genomic_location_.value().genome_id_);
	const auto & tar2Variants = reps[1].microhaplotypes_[0].associated_seq_variants_.value();
	REQUIRE(0 == tar2Variants[0].genome_id_);
	REQUIRE("chr1" == tar2Variants[0].chrom_);
	REQUIRE(1 == tar2Variants[1].genome_id_);
	REQUIRE("chr9" == tar2Variants[1].chrom_);
	REQUIRE(0 == merged.target_info_[1].insert_location_.value().genome_id_);
}

TEST_CASE("mergeMinimumPmos keeps the variants of a haplotype from every PMO once", "[pmoUtils]") {
	auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
	pmo1.targeted_genomes_ = std::vector<GenomeInfo>{genome("g1")};
	pmo1.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].associated_protein_variants_ = std::vector<ProteinVariant>{proteinVariant(0, 7, "K")};
	auto pmo2 = onePmo("tar1", {{"s2", "AAAA", 5}});
	pmo2.targeted_genomes_ = std::vector<GenomeInfo>{genome("g1")};
	pmo2.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].associated_protein_variants_ = std::vector<ProteinVariant>{proteinVariant(0, 7, "K"), proteinVariant(0, 9, "N")};
	auto merged = mergeMinimumPmos({pmo1, pmo2}, mergeMethod());
	REQUIRE(1 == merged.targeted_genomes_.value().size());
	const auto & variants = merged.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].associated_protein_variants_.value();
	REQUIRE(2 == variants.size());
	REQUIRE(7 == variants[0].protein_location_.start_);
	REQUIRE(9 == variants[1].protein_location_.start_);
}

TEST_CASE("mergeMinimumPmos errors on a genome id with no genome", "[pmoUtils]") {
	auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
	pmo1.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].associated_seq_variants_ = std::vector<GenomicLocation>{location(0, "chr1", 5, "T")};
	REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1}, mergeMethod()), Catch::Matchers::ContainsSubstring("refers to genome 0"));
}

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

TEST_CASE("mergeMinimumPmos keeps read counts by stage", "[pmoUtils]") {
	//the samples are in a different order in each PMO
	auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}, {"s2", "AAAA", 20}});
	auto pmo2 = onePmo("tar2", {{"s2", "CCCC", 5}, {"s1", "CCCC", 7}});
	//tar3 has read counts but no haplotypes and is first in target_info, so target ids and representative_microhaplotypes positions differ
	TargetInfo tar3;
	tar3.target_name_ = "tar3";
	tar3.forward_primer_.seq_ = "ACGT";
	tar3.reverse_primer_.seq_ = "TTGG";
	pmo2.target_info_.insert(pmo2.target_info_.begin(), tar3);
	pmo2.representative_microhaplotypes_.targets_[0].target_id_ = 1;
	pmo2.panel_info_[0].reactions_[0].panel_targets_ = {0, 1};
	pmo1.read_counts_by_stage_ = std::vector<ReadCountsByStage>{readCountsFor(pmo1, {0}, 100)};
	pmo2.read_counts_by_stage_ = std::vector<ReadCountsByStage>{readCountsFor(pmo2, {1, 0}, 50)};
	pmo1.read_counts_by_stage_->front().read_counts_by_library_sample_by_stage_[0].extras_["batch"] = "b1";
	pmo2.read_counts_by_stage_->front().read_counts_by_library_sample_by_stage_[1].extras_["lane"] = 2;
	addBioinformatics(pmo1, "run1");
	addBioinformatics(pmo2, "run1");
	pmo1.read_counts_by_stage_->front().bioinformatics_run_id_ = 0;
	pmo2.detected_microhaplotypes_[0].bioinformatics_run_id_ = 0;
	auto merged = mergeMinimumPmos({pmo1, pmo2}, mergeMethod());

	REQUIRE(1 == merged.bioinformatics_run_info_.value().size());
	REQUIRE(1 == merged.bioinformatics_methods_info_.value().size());
	REQUIRE(0 == merged.detected_microhaplotypes_[0].bioinformatics_run_id_.value());
	REQUIRE(1 == merged.read_counts_by_stage_.value().size());
	const auto & counts = merged.read_counts_by_stage_->front();
	REQUIRE(0 == counts.bioinformatics_run_id_.value());
	REQUIRE(2 == counts.read_counts_by_library_sample_by_stage_.size());
	REQUIRE(0 == counts.read_counts_by_library_sample_by_stage_[0].library_sample_id_);
	REQUIRE(1 == counts.read_counts_by_library_sample_by_stage_[1].library_sample_id_);
	for (const auto & libraryCounts : counts.read_counts_by_library_sample_by_stage_) {
		REQUIRE(1000 == libraryCounts.total_raw_count_);
		REQUIRE(3 == libraryCounts.read_counts_for_targets_.value().size());
	}
	//s1 is the first library in pmo1 and the second in pmo2
	REQUIRE(nlohmann::json("b1") == counts.read_counts_by_library_sample_by_stage_[0].extras_.at("batch"));
	REQUIRE(nlohmann::json(2) == counts.read_counts_by_library_sample_by_stage_[0].extras_.at("lane"));
	REQUIRE(0 == counts.read_counts_by_library_sample_by_stage_[1].extras_.count("lane"));
	std::map<std::tuple<std::string, std::string, std::string>, std::pair<uint32_t, uint32_t>> expected{
		{{"s1", "tar1", "extracted"}, {100, 10}}, {{"s2", "tar1", "extracted"}, {100, 10}},
		{{"s1", "tar2", "extracted"}, {50, 5}}, {{"s2", "tar2", "extracted"}, {50, 5}},
		{{"s1", "tar3", "extracted"}, {50, 5}}, {{"s2", "tar3", "extracted"}, {50, 5}}};
	REQUIRE(expected == readCountsOf(merged));

	//through a file the way validatePmo reads it
	auto fnp = bfs::temp_directory_path() / bfs::unique_path("njhseqPmoTest-%%%%-%%%%_PMO.json.gz");
	writePmo(merged, njhseq::OutOptions(fnp, ".gz"));
	auto reread = readPmo(fnp);
	bfs::remove(fnp);
	REQUIRE(expected == readCountsOf(reread));
	REQUIRE(merged.to_json() == reread.to_json());
}

TEST_CASE("mergeMinimumPmos leaves out sections none of the PMOs have", "[pmoUtils]") {
	auto merged = mergeMinimumPmos({onePmo("tar1", {{"s1", "AAAA", 10}}), onePmo("tar2", {{"s1", "CCCC", 5}})}, mergeMethod());
	REQUIRE_FALSE(merged.read_counts_by_stage_.has_value());
	REQUIRE_FALSE(merged.bioinformatics_run_info_.has_value());
	REQUIRE_FALSE(merged.bioinformatics_methods_info_.has_value());
	REQUIRE_FALSE(merged.detected_microhaplotypes_[0].bioinformatics_run_id_.has_value());
	const auto j = merged.to_json();
	REQUIRE_FALSE(j.contains("read_counts_by_stage"));
	REQUIRE_FALSE(j["detected_microhaplotypes"][0].contains("bioinformatics_run_id"));
}

TEST_CASE("mergeMinimumPmos errors on read count conflicts", "[pmoUtils]") {
	auto pmo1 = onePmo("tar1", {{"s1", "AAAA", 10}});
	auto pmo2 = onePmo("tar2", {{"s1", "CCCC", 5}});
	pmo1.read_counts_by_stage_ = std::vector<ReadCountsByStage>{readCountsFor(pmo1, {0}, 100)};
	pmo2.read_counts_by_stage_ = std::vector<ReadCountsByStage>{readCountsFor(pmo2, {0}, 50)};
	SECTION("different total raw counts") {
		pmo2.read_counts_by_stage_->front().read_counts_by_library_sample_by_stage_[0].total_raw_count_ = 999;
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("total_raw_count: 1000 and 999"));
	}
	SECTION("the same library and target twice") {
		pmo1.read_counts_by_stage_->emplace_back(readCountsFor(pmo1, {0}, 100));
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("read counts for target tar1 more than once"));
	}
	SECTION("different library extras") {
		pmo1.read_counts_by_stage_->front().read_counts_by_library_sample_by_stage_[0].extras_["batch"] = "b1";
		pmo2.read_counts_by_stage_->front().read_counts_by_library_sample_by_stage_[0].extras_["batch"] = "b2";
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("conflicting values for batch"));
	}
	SECTION("a target id with no target") {
		pmo2.read_counts_by_stage_ = std::vector<ReadCountsByStage>{readCountsFor(pmo2, {1}, 50)};
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("refer to target 1"));
	}
	SECTION("different bioinformatics runs") {
		addBioinformatics(pmo1, "run1");
		addBioinformatics(pmo2, "run2");
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("different bioinformatics"));
	}
	SECTION("bioinformatics sections in only one PMO") {
		addBioinformatics(pmo2, "run1");
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("different bioinformatics"));
	}
	SECTION("more than one bioinformatics run") {
		addBioinformatics(pmo1, "run1");
		pmo1.bioinformatics_run_info_->emplace_back(pmo1.bioinformatics_run_info_->front());
		pmo1.bioinformatics_run_info_->back().bioinformatics_run_name_ = "run2";
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("at most one bioinformatics run"));
	}
	SECTION("a run id with no run") {
		pmo1.read_counts_by_stage_->front().bioinformatics_run_id_ = 0;
		REQUIRE_THROWS_WITH(mergeMinimumPmos({pmo1, pmo2}, mergeMethod()), Catch::Matchers::ContainsSubstring("refers to bioinformatics run 0"));
	}
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

namespace {
//genTargetInfoFromGenomes style locations and a gff, removed at the end of the test
struct TargetLocationsFixture {
	bfs::path dir_;
	TargetLocationsFixture() :
			dir_(bfs::temp_directory_path() / bfs::unique_path("njhseqPmoTargetGenomes-%%%%-%%%%")) {
		bfs::create_directories(dir_ / "locations");
		bfs::create_directories(dir_ / "gff");
		const std::string header = "#chrom\tfullStart\tfullStop\tname\tlength\tstrand\tgenome\ttarget\tfPrimerStart\tfPrimerStop\terrorsInFPrimer\tinsertStart\tinsertStop\trPrimerStart\trPrimerStop\terrorsInRPrimer\n";
		{
			std::ofstream out((dir_ / "locations" / "GenA_infos.tab.txt").string());
			out << header;
			out << "chrA1\t100\t200\tx\t100\t+\tGenA\ttarA\t100\t120\t0\t120\t180\t180\t200\t0\n";
			out << "chrA1\t500\t600\tx\t100\t-\tGenA\ttarB\t580\t600\t1\t520\t580\t500\t520\t1\n";
			out << "chrA2\t10\t110\tx\t100\t+\tGenA\ttarC\t10\t30\t0\t30\t90\t90\t110\t0\n";
		}
		{
			std::ofstream out((dir_ / "locations" / "GenB_infos.tab.txt").string());
			out << header;
			out << "chrB1\t1000\t1100\tx\t100\t-\tGenB\ttarB\t1080\t1100\t0\t1020\t1080\t1000\t1020\t0\n";
			out << "chrB1\t2000\t2100\tx\t100\t+\tGenB\ttarC\t2000\t2020\t0\t2020\t2080\t2080\t2100\t0\n";
		}
		{
			//tarB's mates, on the reverse strand R1 is at the end of the insert and R2 is on the other strand
			std::ofstream out((dir_ / "locations" / "GenB_inner.bed").string());
			out << "chrB1\t1020\t1080\ttarB\t60\t-\n";
			out << "chrB1\t1020\t1050\ttarB-R2\t30\t+\n";
			out << "chrB1\t1050\t1080\ttarB-R1\t30\t-\n";
			out << "chrB1\t2020\t2080\ttarC\t60\t+\n";
		}
		//a genome with no hits writes an empty table
		std::ofstream((dir_ / "locations" / "GenC_infos.tab.txt").string());
		{
			std::ofstream out((dir_ / "gff" / "GenB.gff").string());
			out << "##gff-version 3\n";
			out << "chrB1\tsrc\tprotein_coding_gene\t1001\t1050\t.\t-\t.\tID=GENE1;Name=kelch13\n";
			out << "chrB1\tsrc\tmRNA\t1001\t1050\t.\t-\t.\tID=GENE1.1;Parent=GENE1\n";
			out << "chrB1\tsrc\tncRNA_gene\t1060\t1200\t.\t+\t.\tID=GENE2\n";
			out << "chrB1\tsrc\tprotein_coding_gene\t1081\t1090\t.\t+\t.\tID=GENE3;Name=outside\n";
		}
	}
	~TargetLocationsFixture() {
		bfs::remove_all(dir_);
	}
};

PortableMicrohaplotypeObject multiTargetPmo() {
	return mergeMinimumPmos({
		onePmo("tarA", {{"s1", "AAAA", 10}}),
		onePmo("tarB", {{"s1", "CCCC", 10}}),
		onePmo("tarB-R1", {{"s1", "GGGG", 10}}),
		onePmo("tarB-R2", {{"s1", "GGCC", 10}}),
		onePmo("tarC", {{"s1", "TTTT", 10}}),
		onePmo("tarMissing", {{"s1", "ACAC", 10}})}, mergeMethod());
}

const TargetInfo & targetNamed(const PortableMicrohaplotypeObject & pmoObj, const std::string & name) {
	return *std::find_if(pmoObj.target_info_.begin(), pmoObj.target_info_.end(), [&name](const TargetInfo & target) {
		return target.target_name_ == name;
	});
}
}  // namespace

TEST_CASE("addTargetGenomeInfo places targets in their genomes with gene names", "[pmoUtils]") {
	TargetLocationsFixture fixture;
	auto pmoObj = multiTargetPmo();
	//a genome already there from variant calling is reused
	GenomeInfo existing;
	existing.name_ = "GenB";
	existing.genome_version_ = "v2";
	existing.url_ = "NA";
	pmoObj.targeted_genomes_ = std::vector<GenomeInfo>{existing};

	AddTargetGenomeInfoPars pars;
	pars.locationsDir_ = fixture.dir_ / "locations";
	pars.gffDir_ = fixture.dir_ / "gff";
	pars.preferredGenomes_ = VecStr{"GenA"};
	auto results = addTargetGenomeInfo(pmoObj, pars);

	REQUIRE(VecStr{"tarMissing"} == results.targetsNotFound_);
	//everything found in GenA goes there
	REQUIRE((std::map<std::string, uint32_t>{{"GenA", 5}}) == results.targetsPerGenome_);
	//no gff for GenA
	REQUIRE(0 == results.targetsWithGenes_);
	REQUIRE(2 == pmoObj.targeted_genomes_->size());
	REQUIRE("GenB" == pmoObj.targeted_genomes_->at(0).name_);
	REQUIRE("v2" == pmoObj.targeted_genomes_->at(0).genome_version_);
	REQUIRE("GenA" == pmoObj.targeted_genomes_->at(1).name_);

	//only in GenA
	const auto & tarA = targetNamed(pmoObj, "tarA");
	REQUIRE(1 == tarA.insert_location_->genome_id_);
	REQUIRE("chrA1" == tarA.insert_location_->chrom_);
	REQUIRE(120 == tarA.insert_location_->start_);
	REQUIRE(180 == tarA.insert_location_->end_);
	REQUIRE(100 == tarA.forward_primer_.location_->start_);
	REQUIRE("-" == tarA.reverse_primer_.location_->strand_.value());
	REQUIRE_FALSE(tarA.gene_name_.has_value());

	//in both, GenA is preferred so it's used even with primer errors
	REQUIRE(1 == targetNamed(pmoObj, "tarB").insert_location_->genome_id_);
	REQUIRE(1 == targetNamed(pmoObj, "tarC").insert_location_->genome_id_);
}

TEST_CASE("addTargetGenomeInfo picks the genome with the fewest primer errors and names overlapping genes", "[pmoUtils]") {
	TargetLocationsFixture fixture;
	auto pmoObj = multiTargetPmo();
	AddTargetGenomeInfoPars pars;
	pars.locationsDir_ = fixture.dir_ / "locations";
	pars.gffDir_ = fixture.dir_ / "gff";
	addTargetGenomeInfo(pmoObj, pars);

	const auto & tarB = targetNamed(pmoObj, "tarB");
	REQUIRE("GenB" == pmoObj.targeted_genomes_->at(tarB.insert_location_->genome_id_).name_);
	REQUIRE("chrB1" == tarB.insert_location_->chrom_);
	REQUIRE(1020 == tarB.insert_location_->start_);
	REQUIRE(1080 == tarB.insert_location_->end_);
	REQUIRE("-" == tarB.insert_location_->strand_.value());
	REQUIRE("-" == tarB.forward_primer_.location_->strand_.value());
	REQUIRE("+" == tarB.reverse_primer_.location_->strand_.value());
	//genes in position order, NA for no Name, the mRNA isn't a gene and GENE3 starts at the insert's end
	REQUIRE("GENE1:kelch13,GENE2:NA" == tarB.gene_name_.value());
	//mates get their window of the target, R2 on the other strand, and the genes of their window
	const auto & tarBR1 = targetNamed(pmoObj, "tarB-R1");
	REQUIRE(1050 == tarBR1.insert_location_->start_);
	REQUIRE(1080 == tarBR1.insert_location_->end_);
	REQUIRE("-" == tarBR1.insert_location_->strand_.value());
	REQUIRE("GENE2:NA" == tarBR1.gene_name_.value());
	const auto & tarBR2 = targetNamed(pmoObj, "tarB-R2");
	REQUIRE(1020 == tarBR2.insert_location_->start_);
	REQUIRE(1050 == tarBR2.insert_location_->end_);
	REQUIRE("+" == tarBR2.insert_location_->strand_.value());
	REQUIRE("GENE1:kelch13" == tarBR2.gene_name_.value());
	REQUIRE(tarB.forward_primer_.location_->to_json() == tarBR2.forward_primer_.location_->to_json());

	//same primer errors in both, by genome name
	REQUIRE("GenA" == pmoObj.targeted_genomes_->at(targetNamed(pmoObj, "tarC").insert_location_->genome_id_).name_);
	//still valid through json
	REQUIRE_NOTHROW(PortableMicrohaplotypeObject::from_json(pmoObj.to_json()).validate());
}

TEST_CASE("chooseTargetGenomes can be limited to the preferred genomes", "[pmoUtils]") {
	TargetLocationsFixture fixture;
	const auto locations = njhseq::TargetGenomeLocations::read(fixture.dir_ / "locations");
	REQUIRE("tarB" == locations.baseTarget("tarB-R1"));
	REQUIRE("tarA" == locations.baseTarget("tarA"));

	njhseq::TargetGenomeChoicePars pars;
	pars.preferredGenomes_ = VecStr{"GenB"};
	SECTION("the others go to their best genome") {
		const auto chosen = njhseq::chooseTargetGenomes(locations, pars);
		REQUIRE(3 == chosen.size());
		REQUIRE("GenA" == chosen.at("tarA").genome_);
		REQUIRE("GenB" == chosen.at("tarB").genome_);
		REQUIRE("GenB" == chosen.at("tarC").genome_);
	}
	SECTION("only the preferred") {
		pars.onlyPreferred_ = true;
		const auto chosen = njhseq::chooseTargetGenomes(locations, pars);
		REQUIRE(2 == chosen.size());
		REQUIRE_FALSE(njh::in(std::string("tarA"), chosen));
		const auto region = locations.innerRegion("tarB-R2", chosen.at("tarB"));
		REQUIRE(nullptr != region);
		REQUIRE('+' == region->strand_);
		REQUIRE(nullptr == locations.innerRegion("tarB-R2", chosen.at("tarC")));
	}
}
