//
//  PmoRunner.cpp
//
//  Programs for working with whole PMO files, add as a ring to give a program these operations
//

#include "PmoRunner.hpp"
#include "njhseq/objects/BioDataObject/pmo/pmoUtils.hpp"

namespace njhseq {

PmoRunner::PmoRunner() :
		njh::progutils::ProgramRunner(
				{ addFunc("mergeMinimumPmos", mergeMinimumPmos, false),
					addFunc("addSpecimenMetaToPmo", addSpecimenMetaToPmo, false),
					addFunc("validatePmo", validatePmo, false),
				}, "Pmo", "1") {
}

int PmoRunner::mergeMinimumPmos(const njh::progutils::CmdArgs & inputCommands) {
	VecStr pmoFnps;
	bfs::path directory = "./";
	std::string contains;
	bool recursive = false;
	std::string programVersion = "1";
	OutOptions outOpts(bfs::path(""), ".json.gz");
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Merge minimum PMOs, for example PMOs per target over the same samples, into one PMO, PMOs with project, sequencing, bioinformatics, genome or read count sections are refused";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --pmos target1_PMO.json.gz,target2_PMO.json.gz --out merged");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --contains _PMO.json.gz --out merged");
	setUp.setOption(pmoFnps, "--pmos", "A list of PMO files to merge", false);
	setUp.setOption(contains, "--contains", "Merge the PMO files in --directory with this in their name", pmoFnps.empty());
	setUp.setOption(directory, "--directory", "Directory to search for PMO files");
	setUp.setOption(recursive, "--recursive", "Search --directory recursively");
	setUp.setOption(programVersion, "--programVersion", "Version of this program to record in the merged PMO");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	if (pmoFnps.empty()) {
		auto files = njh::files::listAllFiles(directory, recursive, VecStr{contains});
		for (const auto & file : files) {
			if (!file.second) {
				pmoFnps.emplace_back(file.first.string());
			}
		}
		if (pmoFnps.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, found no files containing " << contains << " in " << directory << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	njh::sort(pmoFnps);
	pmo::PmoGenerationMethod generationMethod;
	generationMethod.program_name_ = njh::pasteAsStr(setUp.commands_.masterProgram_, " mergeMinimumPmos");
	generationMethod.program_version_ = programVersion;
	pmo::MinimumPmoMerger merger(generationMethod);
	for (const auto & pmoFnp : pmoFnps) {
		merger.add(pmo::readPmo(pmoFnp));
	}
	pmo::writePmo(merger.finish(), outOpts);
	return 0;
}

int PmoRunner::addSpecimenMetaToPmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmoFnp;
	bfs::path metaFnp;
	std::string sampleCol = "sample";
	OutOptions outOpts(bfs::path(""), ".json.gz");
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Add the columns of a table as specimen fields in a PMO, columns that aren't in the PMO schema are added as custom fields";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --pmo merged_PMO.json.gz --meta meta.tsv --out withMeta");
	setUp.setOption(pmoFnp, "--pmo", "The PMO file", true);
	setUp.setOption(metaFnp, "--meta", "A tab delimited table with a header", true);
	setUp.setOption(sampleCol, "--sampleCol", "The column with the library sample or specimen name");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	auto pmoObj = pmo::readPmo(pmoFnp);
	table meta(metaFnp.string(), "\t", true);
	auto results = pmo::addSpecimenMeta(pmoObj, meta, sampleCol);
	std::cout << "specimens updated: " << results.specimensUpdated_ << std::endl;
	if (!results.metaSamplesNotInPmo_.empty()) {
		std::cerr << results.metaSamplesNotInPmo_.size() << " samples in " << metaFnp << " aren't in the PMO: " << njh::conToStr(results.metaSamplesNotInPmo_, ",") << std::endl;
	}
	pmo::writePmo(pmoObj, outOpts);
	return 0;
}

int PmoRunner::validatePmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmoFnp;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Check a PMO against the format and print what's in it";
	setUp.setOption(pmoFnp, "--pmo", "The PMO file", true);
	setUp.finishSetUp(std::cout);

	auto pmoObj = pmo::readPmo(pmoFnp);
	uint32_t librariesWithResults = 0;
	for (const auto & run : pmoObj.detected_microhaplotypes_) {
		librariesWithResults += run.library_samples_.size();
	}
	std::cout << "pmo_version\t" << pmoObj.pmo_header_.pmo_version_ << std::endl;
	std::cout << "panels\t" << pmoObj.panel_info_.size() << std::endl;
	std::cout << "targets\t" << pmoObj.target_info_.size() << std::endl;
	std::cout << "specimens\t" << pmoObj.specimen_info_.size() << std::endl;
	std::cout << "library_samples\t" << pmoObj.library_sample_info_.size() << std::endl;
	std::cout << "library_samples_with_results\t" << librariesWithResults << std::endl;
	return 0;
}

}  // namespace njhseq
