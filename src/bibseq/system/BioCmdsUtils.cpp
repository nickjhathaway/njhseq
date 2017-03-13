/*
 * BioCmdsUtils.cpp
 *
 *  Created on: Mar 13, 2017
 *      Author: nick
 */

#include "BioCmdsUtils.hpp"

namespace bibseq {

BioCmdsUtils::BioCmdsUtils() {
}
BioCmdsUtils::BioCmdsUtils(bool verbose) :
		verbose_(verbose) {
}
bool verbose_ = false;

bib::sys::RunOutput BioCmdsUtils::RunBowtie2Index(const bfs::path & genomeFnp) {
	bib::sys::requireExternalProgramThrow("bowtie2-build");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	std::string templateCmd = "bowtie2-build -q " + genomeFnp.string() + " "
			+ prefix.string();
	auto bowtie2CheckFile = prefix.string() + ".1.bt2";
	return runCmdCheck(templateCmd, genomeFnp, bowtie2CheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunBwaIndex(const bfs::path & genomeFnp) {
	bib::sys::requireExternalProgramThrow("bwa");
	std::string templateCmd = "bwa index " + genomeFnp.string();
	auto bwaCheckFile = genomeFnp.string() + ".bwt";
	return runCmdCheck(templateCmd, genomeFnp, bwaCheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunSamtoolsFastaIndex(const bfs::path & genomeFnp) {
	bib::sys::requireExternalProgramThrow("samtools");
	std::string templateCmd = "samtools faidx " + genomeFnp.string();
	auto samCheckFile = genomeFnp.string() + ".fai";
	return runCmdCheck(templateCmd, genomeFnp, samCheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunPicardFastaSeqDict(const bfs::path & genomeFnp) {
	bib::sys::requireExternalProgramThrow("picard");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".dict";
	std::string templateCmd = "picard CreateSequenceDictionary R="
			+ genomeFnp.string() + " O=" + outFile;
	return runCmdCheck(templateCmd, genomeFnp, outFile);
}

bib::sys::RunOutput BioCmdsUtils::RunFaToTwoBit(const bfs::path & genomeFnp) {
	bib::sys::requireExternalProgramThrow("TwoBit");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".2bit";
	std::string templateCmd = "TwoBit faToTwoBit --in " + genomeFnp.string()
			+ " --out " + outFile + " --overWrite";
	return runCmdCheck(templateCmd, genomeFnp, outFile);
}

bib::sys::RunOutput BioCmdsUtils::runCmdCheck(const std::string & cmd,
		const bfs::path & input, const bfs::path & check) {
	if (!bfs::exists(check) || bib::files::firstFileIsOlder(check, input)) {
		if (verbose_) {
			std::cout << bib::bashCT::bold << bib::bashCT::blue << "Running: "
					<< bib::bashCT::green << cmd << std::endl;
		}
		if(bfs::exists(check)){
			if (verbose_) {
				std::cout << bib::bashCT::bold << bib::bashCT::blue << "\tRemoving out of date file: "
										<< bib::bashCT::green << check << std::endl;
			}
			bfs::remove(check);
		}
		auto runOut = bib::sys::run( { cmd });
		checkRunOutThrow(runOut, __PRETTY_FUNCTION__);
		return runOut;
	} else {
		if (verbose_) {
			std::cout << bib::bashCT::bold << bib::bashCT::blue
					<< "No need to run: " << bib::bashCT::green << cmd << std::endl;
		}
	}
	return bib::sys::RunOutput();
}

void BioCmdsUtils::checkRunOutThrow(const bib::sys::RunOutput & runOut,
		const std::string & funcName) {
	if (!runOut.success_) {
		std::stringstream ss;
		ss << funcName << " error in running " << runOut.cmd_ << "\n";
		if ("" != runOut.stdOut_) {
			ss << runOut.stdOut_ << "\n";
		}
		if ("" != runOut.stdErr_) {
			ss << runOut.stdErr_ << "\n";
		}
		throw std::runtime_error { ss.str() };
	}
}
} /* namespace bibseq */
