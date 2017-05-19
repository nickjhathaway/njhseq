/*
 * BioCmdsUtils.cpp
 *
 *  Created on: Mar 13, 2017
 *      Author: nick
 */

#include "BioCmdsUtils.hpp"
#include <unordered_map>

namespace bibseq {

BioCmdsUtils::BioCmdsUtils() {
}
BioCmdsUtils::BioCmdsUtils(bool verbose) :
		verbose_(verbose) {
}

void checkGenomeFnpExistsThrow(const bfs::path & genomeFnp,
		const std::string & funcName) {
	if (!bfs::exists(genomeFnp)) {
		std::stringstream ss;
		ss << funcName << ", error, genome file: " << bib::bashCT::blue << genomeFnp
				<< bib::bashCT::reset << bib::bashCT::red << ", doesn't exists "
				<< bib::bashCT::reset << "\n";
		throw std::runtime_error { ss.str() };
	}
}

bib::sys::RunOutput BioCmdsUtils::RunBowtie2Index(const bfs::path & genomeFnp) {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("bowtie2-build");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	std::string templateCmd = "bowtie2-build -q " + genomeFnp.string() + " "
			+ prefix.string();
	auto bowtie2CheckFile = prefix.string() + ".1.bt2";
	return runCmdCheck(templateCmd, genomeFnp, bowtie2CheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunBwaIndex(const bfs::path & genomeFnp) {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("bwa");
	std::string templateCmd = "bwa index " + genomeFnp.string();
	auto bwaCheckFile = genomeFnp.string() + ".bwt";
	return runCmdCheck(templateCmd, genomeFnp, bwaCheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunSamtoolsFastaIndex(const bfs::path & genomeFnp) {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("samtools");
	std::string templateCmd = "samtools faidx " + genomeFnp.string();
	auto samCheckFile = genomeFnp.string() + ".fai";
	return runCmdCheck(templateCmd, genomeFnp, samCheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunPicardFastaSeqDict(const bfs::path & genomeFnp) {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("picard");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".dict";
	std::string templateCmd = "picard CreateSequenceDictionary R="
			+ genomeFnp.string() + " O=" + outFile;
	return runCmdCheck(templateCmd, genomeFnp, outFile);
}

bib::sys::RunOutput BioCmdsUtils::RunFaToTwoBit(const bfs::path & genomeFnp) {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("TwoBit");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".2bit";
	std::string templateCmd = "TwoBit faToTwoBit --in " + genomeFnp.string()
			+ " --out " + outFile + " --overWrite";
	return runCmdCheck(templateCmd, genomeFnp, outFile);
}


std::unordered_map<std::string, bib::sys::RunOutput> BioCmdsUtils::runAllPossibleIndexes(const bfs::path & genomeFnp){
	std::unordered_map<std::string, bib::sys::RunOutput> outputs;
	if(bib::sys::hasSysCommand("bowtie2")){
		outputs.emplace("bowtie2", RunBowtie2Index(genomeFnp));
	}
	if(bib::sys::hasSysCommand("bwa")){
		outputs.emplace("bwa", RunBwaIndex(genomeFnp));
	}
	if(bib::sys::hasSysCommand("samtools")){
		outputs.emplace("samtools", RunSamtoolsFastaIndex(genomeFnp));
	}
	if(bib::sys::hasSysCommand("picard")){
		outputs.emplace("picard", RunPicardFastaSeqDict(genomeFnp));
	}
	if(bib::sys::hasSysCommand("TwoBit")){
		outputs.emplace("TwoBit", RunFaToTwoBit(genomeFnp));
	}
	return outputs;
}

bib::sys::RunOutput BioCmdsUtils::runCmdCheck(const std::string & cmd,
		const bfs::path & input, const bfs::path & check) {
	if (!bfs::exists(check) || bib::files::firstFileIsOlder(check, input)) {
		if (verbose_) {
			std::cout << bib::bashCT::bold << bib::bashCT::blue << "Running: "
					<< bib::bashCT::green << cmd << bib::bashCT::reset << std::endl;
		}
		if(bfs::exists(check)){
			if (verbose_) {
				std::cout << bib::bashCT::bold << bib::bashCT::blue << "\tRemoving out of date file: "
										<< bib::bashCT::green << check  << bib::bashCT::reset << std::endl;
			}
			bfs::remove(check);
		}
		auto runOut = bib::sys::run( { cmd });
		checkRunOutThrow(runOut, __PRETTY_FUNCTION__);
		return runOut;
	} else {
		if (verbose_) {
			std::cout << bib::bashCT::bold << bib::bashCT::blue
					<< "No need to run: " << bib::bashCT::green << cmd << bib::bashCT::reset << std::endl;
		}
	}
	return bib::sys::RunOutput();
}

bib::sys::RunOutput BioCmdsUtils::bowtie2Align(const SeqIOOptions & opts,
		const bfs::path & genomeFnp, std::string additionalBowtie2Args) {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bfs::path outputFnp = bib::appendAsNeededRet(opts.out_.outFilename_.string(),
			".sorted.bam");
	if (bfs::exists(outputFnp) && !opts.out_.overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, " << outputFnp
				<< " already exists, use --overWrite to over write it" << "\n";
		throw std::runtime_error { ss.str() };
	}
	RunBowtie2Index(genomeFnp);
	auto genomePrefix = genomeFnp;
	genomePrefix.replace_extension("");
	std::stringstream templateCmd;
	if(SeqIOOptions::inFormats::FASTA == opts.inFormat_){
		additionalBowtie2Args += " -f ";
	}
	templateCmd << "bowtie2 -U " << opts.firstName_ << " -x " << genomePrefix
			<< " " <<  additionalBowtie2Args << " | samtools view - -b | samtools sort - -o " << outputFnp
			<< " " << "&& samtools index " << outputFnp;
	if (verbose_) {
		std::cout << "Running: " << bib::bashCT::green << templateCmd.str()
				<< bib::bashCT::reset << std::endl;
	}
	auto ret = bib::sys::run( { templateCmd.str() });
	BioCmdsUtils::checkRunOutThrow(ret, __PRETTY_FUNCTION__);
	return ret;
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
