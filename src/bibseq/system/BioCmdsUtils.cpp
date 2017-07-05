/*
 * BioCmdsUtils.cpp
 *
 *  Created on: Mar 13, 2017
 *      Author: nick
 */

#include "BioCmdsUtils.hpp"
#include <unordered_map>

#include "bibseq/IO/SeqIO.h"
#include "bibseq/IO/OutputStream.hpp"

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

bib::sys::RunOutput BioCmdsUtils::RunBowtie2Index(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("bowtie2-build");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	std::string templateCmd = "bowtie2-build -q " + genomeFnp.string() + " "
			+ prefix.string();
	auto bowtie2CheckFile = prefix.string() + ".1.bt2";
	return runCmdCheck(templateCmd, genomeFnp, bowtie2CheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunBwaIndex(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("bwa");
	std::string templateCmd = "bwa index " + genomeFnp.string();
	auto bwaCheckFile = genomeFnp.string() + ".bwt";
	return runCmdCheck(templateCmd, genomeFnp, bwaCheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunSamtoolsFastaIndex(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("samtools");
	std::string templateCmd = "samtools faidx " + genomeFnp.string();
	auto samCheckFile = genomeFnp.string() + ".fai";
	return runCmdCheck(templateCmd, genomeFnp, samCheckFile);
}

bib::sys::RunOutput BioCmdsUtils::RunPicardFastaSeqDict(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("picard");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".dict";
	std::string templateCmd = "picard CreateSequenceDictionary R="
			+ genomeFnp.string() + " O=" + outFile;
	return runCmdCheck(templateCmd, genomeFnp, outFile);
}

bib::sys::RunOutput BioCmdsUtils::RunFaToTwoBit(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bib::sys::requireExternalProgramThrow("TwoBit");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".2bit";
	std::string templateCmd = "TwoBit faToTwoBit --in " + genomeFnp.string()
			+ " --out " + outFile + " --overWrite";
	return runCmdCheck(templateCmd, genomeFnp, outFile);
}


std::unordered_map<std::string, bib::sys::RunOutput> BioCmdsUtils::runAllPossibleIndexes(const bfs::path & genomeFnp) const {
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
		const bfs::path & input, const bfs::path & check) const {
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
		const bfs::path & genomeFnp, std::string additionalBowtie2Args) const {

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

bib::sys::RunOutput BioCmdsUtils::lastzAlign(const SeqIOOptions & opts, const LastZPars & pars) const {
	checkGenomeFnpExistsThrow(pars.genomeFnp, __PRETTY_FUNCTION__);
	bfs::path outputFnp = bib::appendAsNeededRet(opts.out_.outFilename_.string(),
			".sorted.bam");
	if (bfs::exists(outputFnp) && !opts.out_.overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, " << outputFnp
				<< " already exists, use --overWrite to over write it" << "\n";
		throw std::runtime_error { ss.str() };
	}
	std::stringstream lastzCmd;
	lastzCmd << "lastz " << pars.genomeFnp << "[multiple] "
			<< opts.firstName_ << " --format=" << pars.outFormat
			<< " --coverage=" << pars.coverage << " --identity=" << pars.identity << " "
			<< pars.extraLastzArgs << " | samtools view - -b | samtools sort - -o "
			<< outputFnp << " " << "&& samtools index "
			<< outputFnp;
	if (verbose_) {
		std::cout << "Running: " << bib::bashCT::green << lastzCmd.str()
				<< bib::bashCT::reset << std::endl;
	}
	auto ret = bib::sys::run( { lastzCmd.str() });
	BioCmdsUtils::checkRunOutThrow(ret, __PRETTY_FUNCTION__);
	return ret;
}


BioCmdsUtils::FastqDumpResults BioCmdsUtils::runFastqDump(const FastqDumpPars & pars) const{
	BioCmdsUtils::FastqDumpResults ret;
	if(!bfs::exists(pars.sraFnp_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << pars.sraFnp_ << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}
	bib::sys::requireExternalProgramThrow(fastqDumpCmd_);
	std::stringstream cmd;
	cmd <<  fastqDumpCmd_ << " --log-level 0 -X 1 -Z --split-spot " << pars.sraFnp_;
	auto cmdOutput = bib::sys::run({cmd.str()});
	BioCmdsUtils::checkRunOutThrow(cmdOutput, __PRETTY_FUNCTION__);
	//bib::sys::run trim end white space so have to add one;
	uint32_t newLines = countOccurences(cmdOutput.stdOut_, "\n") + 1;

	std::string extraSraArgs = pars.extraSraOptions_;
	if (pars.exportBarCode_) {
		extraSraArgs = bib::replaceString(extraSraArgs, "--skip-technical", "");
	} else if (!bib::containsSubString(pars.extraSraOptions_,
			"--skip-technical")) {
		extraSraArgs += " --skip-technical ";
	}
	bfs::path checkFile1        = bib::files::replaceExtension(pars.sraFnp_, "_1.fastq");
	bfs::path checkFile2        = bib::files::replaceExtension(pars.sraFnp_, "_2.fastq");
	bfs::path checkFileBarcodes = bib::files::replaceExtension(pars.sraFnp_, "_barcodes.fastq");
	if(pars.gzip_){
		checkFile1 = checkFile1.string() + ".gz";
		checkFile2 = checkFile2.string() + ".gz";
		checkFileBarcodes = checkFileBarcodes.string() + ".gz";
		ret.isGzipped_ = true;
	}
	bool needToRun = true;
	ret.firstMateFnp_ = checkFile1;
	if (4 == newLines) {
		ret.isPairedEnd_ = false;
		if(bfs::exists(checkFile1) &&
				bib::files::firstFileIsOlder(pars.sraFnp_, checkFile1)){
			needToRun = false;
		}
	} else if (8 == newLines) {
		ret.secondMateFnp_ = checkFile2;
		ret.isPairedEnd_ = true;
		if(bfs::exists(checkFile1) &&
				bfs::exists(checkFile2) &&
				bib::files::firstFileIsOlder(pars.sraFnp_, checkFile1)&&
				bib::files::firstFileIsOlder(pars.sraFnp_, checkFile2)){
			needToRun = false;
		}
	} else if (12 == newLines) {
		ret.secondMateFnp_ = checkFile2;
		if(pars.exportBarCode_){
			ret.barcodeFnp_ = checkFileBarcodes;
			if(bfs::exists(checkFile1) &&
					bfs::exists(checkFile2) &&
					bfs::exists(checkFileBarcodes) &&
					bib::files::firstFileIsOlder(pars.sraFnp_, checkFile1)&&
					bib::files::firstFileIsOlder(pars.sraFnp_, checkFile2)&&
					bib::files::firstFileIsOlder(pars.sraFnp_, checkFileBarcodes)){
				needToRun = false;
			}
		}else{
			if(bfs::exists(checkFile1) &&
					bfs::exists(checkFile2) &&
					bib::files::firstFileIsOlder(pars.sraFnp_, checkFile1)&&
					bib::files::firstFileIsOlder(pars.sraFnp_, checkFile2)){
				needToRun = false;
			}
		}
		ret.isPairedEnd_ = true;
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error with " << cmd.str()
				<< " was expecting 4, 8, or 12 lines but got " << newLines << " instead "
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	if(needToRun){
		std::stringstream ss;
		if(12 != newLines){
			extraSraArgs += " --gzip ";
		}
		ss << fastqDumpCmd_ << " --split-files --defline-seq '@$sn/$ri' " << extraSraArgs << " " << pars.sraFnp_;
		auto runOutput = bib::sys::run({ss.str()});
		checkRunOutThrow(runOutput, __PRETTY_FUNCTION__);
		if(12 == newLines){
			//now to fix the crazy that is the SRA way of dump things with three files
			if(pars.exportBarCode_){
				bfs::path currentBarcodeFnp = bib::files::replaceExtension(pars.sraFnp_, "_2.fastq");
				if(pars.gzip_){
					IoOptions barIoOpts{InOptions(currentBarcodeFnp), OutOptions(checkFileBarcodes)};
					barIoOpts.out_.overWriteFile_ = true;
					gzZipFile(barIoOpts);
				}else{
					if(bfs::exists(checkFileBarcodes)){
						bfs::remove(checkFileBarcodes);
					}
					bfs::copy(currentBarcodeFnp, checkFileBarcodes);
				}
			}
			if(pars.gzip_){
				auto currentFirstMateFnp = bib::files::replaceExtension(pars.sraFnp_, "_1.fastq");
				IoOptions firstMateIoOpts{InOptions(currentFirstMateFnp), OutOptions(checkFile1)};
				firstMateIoOpts.out_.overWriteFile_ = true;
				gzZipFile(firstMateIoOpts);
			}
			bfs::path currentSecondMateFnp = bib::files::replaceExtension(pars.sraFnp_, "_3.fastq");
			auto secondMateIn = SeqIOOptions::genFastqIn(currentSecondMateFnp);
			SeqInputExp reader{secondMateIn};
			seqInfo seq;
			OutOptions secondMateOutOpts(checkFile2);
			OutputStream secondMateOut(secondMateOutOpts);
			reader.openIn();
			while(reader.readNextRead(seq)){
				seq.name_.back() = '2';
				seq.outPutFastq(secondMateOut);
			}
		}
	}

	return ret;
}

bool BioCmdsUtils::isSRAPairedEnd(const bfs::path & sraFnp) const{
	if(!bfs::exists(sraFnp)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << sraFnp << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}
	bib::sys::requireExternalProgramThrow(fastqDumpCmd_);
	std::stringstream cmd;
	cmd <<  fastqDumpCmd_ << " --log-level 0 -X 1 -Z --split-spot " << sraFnp;
	auto cmdOutput = bib::sys::run({cmd.str()});
	BioCmdsUtils::checkRunOutThrow(cmdOutput, __PRETTY_FUNCTION__);
	//bib::sys::run trim end white space so have to add one;
	uint32_t newLines = countOccurences(cmdOutput.stdOut_, "\n") + 1;
	if(8 == newLines || 12 == newLines){
		//check for 8 or 12, 12 means triple file _1 forward mate, _2 barcode, _3 reverse mate
		//could also be solved by checking for only 8 and using --skip-technical in the command to skip the barcode file
		return true;
	}
	return false;
}

void BioCmdsUtils::checkRunOutThrow(const bib::sys::RunOutput & runOut,
		const std::string & funcName)  {
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
