#include <catch_amalgamated.hpp>
#include "njhseq/objects/collapseObjects/SampleCollapseCollection.hpp"
#include "njhseq/objects/BioDataObject/pmo.h"
#include "njhseq/IO/OutputStream.hpp"

using njhseq::VecStr;
using njhseq::collapse::SampleCollapseCollection;
namespace bfs = njh::bfs;

namespace {
//a fresh output directory with a selectedClustersInfo table, removed at the end of the test
struct PmoFixture {
	bfs::path dir_;
	explicit PmoFixture(const VecStr & rows) :
			dir_(bfs::temp_directory_path() / bfs::unique_path("njhseqPmoTest-%%%%-%%%%")) {
		bfs::create_directories(dir_);
		njhseq::OutOptions outOpts(dir_ / "selectedClustersInfo.tab.txt.gz", ".gz");
		njhseq::OutputStream out(outOpts);
		out << "s_Sample\tp_name\th_popUID\th_Consensus\tc_ReadCnt\n";
		for (const auto & row : rows) {
			out << row << "\n";
		}
	}
	~PmoFixture() {
		bfs::remove_all(dir_);
	}
	SampleCollapseCollection collection(const VecStr & samples) const {
		return SampleCollapseCollection(njhseq::SeqIOOptions(), dir_, dir_,
				njhseq::PopNamesInfo("tar1", samples, VecStr{}),
				SampleCollapseCollection::PreFilteringCutOffs());
	}
};

SampleCollapseCollection::MinimumPmoPars pmoPars() {
	SampleCollapseCollection::MinimumPmoPars pars;
	pars.forwardPrimer_ = "ACGT";
	pars.reversePrimer_ = "TTGG";
	pars.programName_ = "SeekDeep processClusters";
	pars.programVersion_ = "3.0.1";
	return pars;
}
}  // namespace

TEST_CASE("minimum PMO holds the population haplotypes and reads per sample", "[minimumPMO]") {
	PmoFixture fixture({
		"s1\ttar1\ttar1.01\tTTTT\t10.4",
		"s1\ttar1\ttar1.00\tAAAA\t90.6",
		"s2\ttar1\ttar1.00\tAAAA\t5"});
	auto pmoObj = fixture.collection({"s1", "s2"}).genMinimumPMO(pmoPars());

	REQUIRE(1 == pmoObj.target_info_.size());
	REQUIRE("tar1" == pmoObj.target_info_[0].target_name_);
	REQUIRE("ACGT" == pmoObj.target_info_[0].forward_primer_.seq_);
	REQUIRE("TTGG" == pmoObj.target_info_[0].reverse_primer_.seq_);
	REQUIRE("tar1" == pmoObj.panel_info_[0].panel_name_);
	REQUIRE(VecStr{"s1", "s2"} == VecStr{pmoObj.library_sample_info_[0].library_sample_name_, pmoObj.library_sample_info_[1].library_sample_name_});
	REQUIRE(1 == pmoObj.library_sample_info_[1].specimen_id_);
	REQUIRE("s2" == pmoObj.specimen_info_[1].specimen_name_);

	//representatives in population order
	const auto & mhaps = pmoObj.representative_microhaplotypes_.targets_[0].microhaplotypes_;
	REQUIRE(2 == mhaps.size());
	REQUIRE("AAAA" == mhaps[0].seq_);
	REQUIRE("TTTT" == mhaps[1].seq_);
	REQUIRE("tar1.00" == mhaps[0].microhaplotype_name_.value());
	REQUIRE("tar1.01" == mhaps[1].microhaplotype_name_.value());
	REQUIRE(njhseq::pmo::PMO_SCHEMA_VERSION == pmoObj.pmo_header_.pmo_version_);

	const auto & samples = pmoObj.detected_microhaplotypes_[0].library_samples_;
	REQUIRE(2 == samples.size());
	const auto & s1Mhaps = samples[0].target_results_[0].mhaps_;
	REQUIRE(2 == s1Mhaps.size());
	REQUIRE(0 == s1Mhaps[0].mhap_id_);
	REQUIRE(91 == s1Mhaps[0].reads_);
	REQUIRE(1 == s1Mhaps[1].mhap_id_);
	REQUIRE(10 == s1Mhaps[1].reads_);
	REQUIRE(1 == samples[1].library_sample_id_);
}

TEST_CASE("minimum PMO lists samples without final haplotypes but not as detected", "[minimumPMO]") {
	PmoFixture fixture({"s2\ttar1\ttar1.00\tAAAA\t5"});
	auto pmoObj = fixture.collection({"s1", "s2"}).genMinimumPMO(pmoPars());
	REQUIRE(2 == pmoObj.library_sample_info_.size());
	REQUIRE(2 == pmoObj.specimen_info_.size());
	REQUIRE(1 == pmoObj.detected_microhaplotypes_[0].library_samples_.size());
	REQUIRE(1 == pmoObj.detected_microhaplotypes_[0].library_samples_[0].library_sample_id_);
}

TEST_CASE("minimum PMO combines population haplotypes with the same sequence", "[minimumPMO]") {
	PmoFixture fixture({
		"s1\ttar1\ttar1.00\tAAAA\t10.25",
		"s1\ttar1\ttar1.01\tAAAA\t10.25"});
	auto pmoObj = fixture.collection({"s1"}).genMinimumPMO(pmoPars());
	REQUIRE(1 == pmoObj.representative_microhaplotypes_.targets_[0].microhaplotypes_.size());
	const auto & s1Mhaps = pmoObj.detected_microhaplotypes_[0].library_samples_[0].target_results_[0].mhaps_;
	REQUIRE(1 == s1Mhaps.size());
	REQUIRE(21 == s1Mhaps[0].reads_);
}

TEST_CASE("minimum PMO with no final haplotypes is still valid", "[minimumPMO]") {
	PmoFixture fixture(VecStr{});
	auto pmoObj = fixture.collection({"s1"}).genMinimumPMO(pmoPars());
	REQUIRE(1 == pmoObj.library_sample_info_.size());
	REQUIRE(pmoObj.representative_microhaplotypes_.targets_[0].microhaplotypes_.empty());
	REQUIRE(pmoObj.detected_microhaplotypes_[0].library_samples_.empty());
}

TEST_CASE("minimum PMO panel name can be set", "[minimumPMO]") {
	PmoFixture fixture({"s1\ttar1\ttar1.00\tAAAA\t5"});
	auto pars = pmoPars();
	pars.panelName_ = "panelA";
	auto pmoObj = fixture.collection({"s1"}).genMinimumPMO(pars);
	REQUIRE("panelA" == pmoObj.panel_info_[0].panel_name_);
}

TEST_CASE("minimum PMO fails validation on names the schema doesn't allow", "[minimumPMO]") {
	PmoFixture fixture({"bad#name\ttar1\ttar1.00\tAAAA\t5"});
	REQUIRE_THROWS(fixture.collection({"bad#name"}).genMinimumPMO(pmoPars()));
}

TEST_CASE("minimum PMO marks fake primers in the target attributes", "[minimumPMO]") {
	PmoFixture fixture({"s1\ttar1\ttar1.00\tAAAA\t5"});
	auto pars = pmoPars();
	REQUIRE_FALSE(fixture.collection({"s1"}).genMinimumPMO(pars).target_info_[0].target_attributes_.has_value());
	pars.fakePrimers_ = true;
	auto pmoObj = fixture.collection({"s1"}).genMinimumPMO(pars);
	REQUIRE(VecStr{njhseq::pmo::PMO_FAKE_PRIMERS_ATTRIBUTE} == pmoObj.target_info_[0].target_attributes_.value());
}

TEST_CASE("minimum PMO round trips through a gzipped file", "[minimumPMO]") {
	PmoFixture fixture({"s1\ttar1\ttar1.00\tAAAA\t5"});
	auto fnp = fixture.dir_ / "PMO.json.gz";
	auto pars = pmoPars();
	fixture.collection({"s1"}).writeMinimumPMO(fnp, pars);
	njhseq::InputStream in{njhseq::InOptions(fnp)};
	auto pmoObj = njhseq::pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	pmoObj.validate();
	REQUIRE("AAAA" == pmoObj.representative_microhaplotypes_.targets_[0].microhaplotypes_[0].seq_);
	REQUIRE(5 == pmoObj.detected_microhaplotypes_[0].library_samples_[0].target_results_[0].mhaps_[0].reads_);
}
