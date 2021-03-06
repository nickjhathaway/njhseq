/*
 * BioCmdsUtils.cpp
 *
 *  Created on: Mar 13, 2017
 *      Author: nick
 */
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "BioCmdsUtils.hpp"
#include <unordered_map>

#include "njhseq/IO/SeqIO.h"
#include "njhseq/IO/OutputStream.hpp"

#include <TwoBit.h>

namespace njhseq {

BioCmdsUtils::BioCmdsUtils() {
}
BioCmdsUtils::BioCmdsUtils(bool verbose) :
		verbose_(verbose) {
}

void checkGenomeFnpExistsThrow(const bfs::path & genomeFnp,
		const std::string & funcName) {
	if (!bfs::exists(genomeFnp)) {
		std::stringstream ss;
		ss << funcName << ", error, genome file: " << njh::bashCT::blue << genomeFnp
				<< njh::bashCT::reset << njh::bashCT::red << ", doesn't exists "
				<< njh::bashCT::reset << "\n";
		throw std::runtime_error { ss.str() };
	}
}

njh::sys::RunOutput BioCmdsUtils::RunBowtie2Index(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	njh::sys::requireExternalProgramThrow("bowtie2-build");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	std::string templateCmd = "bowtie2-build -q " + genomeFnp.string() + " "
			+ prefix.string();
	auto bowtie2CheckFile = prefix.string() + ".1.bt2";
	return runCmdCheck(templateCmd, genomeFnp, bowtie2CheckFile);
}

njh::sys::RunOutput BioCmdsUtils::RunBwaIndex(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	njh::sys::requireExternalProgramThrow("bwa");
	std::string templateCmd = "bwa index " + genomeFnp.string();
	auto bwaCheckFile = genomeFnp.string() + ".bwt";
	return runCmdCheck(templateCmd, genomeFnp, bwaCheckFile);
}

njh::sys::RunOutput BioCmdsUtils::RunSamtoolsFastaIndex(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	njh::sys::requireExternalProgramThrow("samtools");
	std::string templateCmd = "samtools faidx " + genomeFnp.string();
	auto samCheckFile = genomeFnp.string() + ".fai";
	return runCmdCheck(templateCmd, genomeFnp, samCheckFile);
}

njh::sys::RunOutput BioCmdsUtils::RunPicardFastaSeqDict(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	njh::sys::requireExternalProgramThrow("picard");
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".dict";
	std::string templateCmd = "picard CreateSequenceDictionary R="
			+ genomeFnp.string() + " O=" + outFile;
	return runCmdCheck(templateCmd, genomeFnp, outFile);
}

njh::sys::RunOutput BioCmdsUtils::RunFaToTwoBit(const bfs::path & genomeFnp) const {
	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	TwoBit::faToTwoBitPars pars;
	pars.inputFilename = genomeFnp.string();
	auto prefix = bfs::path(genomeFnp).replace_extension("");
	auto outFile = prefix.string() + ".2bit";
	pars.outFilename = outFile;
	if(!bfs::exists(pars.outFilename) || njh::files::firstFileIsOlder(pars.outFilename, genomeFnp)){
		pars.overWrite = true;
		TwoBit::fastasToTwoBit(pars);
	}
	njh::sys::RunOutput fakeOut;
	fakeOut.success_ = true;
	return fakeOut;
}


std::unordered_map<std::string, njh::sys::RunOutput> BioCmdsUtils::runAllPossibleIndexes(
		const bfs::path & genomeFnp) const {
	/**@todo add a force option to force the indexing even if it isn't needed*/
	std::unordered_map<std::string, njh::sys::RunOutput> outputs;

	if (njh::sys::hasSysCommand("bowtie2")) {
		outputs.emplace("bowtie2", RunBowtie2Index(genomeFnp));
	}else	if(verbose_){
		std::cerr << "Couldn't find " << "bowtie2" << " skipping bowtie2 indexing" << std::endl;
	}

	if (njh::sys::hasSysCommand("bwa")) {
		outputs.emplace("bwa", RunBwaIndex(genomeFnp));
	}else	if(verbose_){
		std::cerr << "Couldn't find " << "bwa" << " skipping bwa indexing" << std::endl;
	}

	if (njh::sys::hasSysCommand("samtools")) {
		outputs.emplace("samtools", RunSamtoolsFastaIndex(genomeFnp));
	}else	if(verbose_){
		std::cerr << "Couldn't find " << "samtools" << " skipping samtools faidx" << std::endl;
	}

	if (njh::sys::hasSysCommand("picard")) {
		outputs.emplace("picard", RunPicardFastaSeqDict(genomeFnp));
	}else	if(verbose_){
		std::cerr << "Couldn't find " << "picard" << " skipping picard CreateSequenceDictionary" << std::endl;
	}
	outputs.emplace("TwoBit", RunFaToTwoBit(genomeFnp));
//	if (njh::sys::hasSysCommand("TwoBit")) {
//		outputs.emplace("TwoBit", RunFaToTwoBit(genomeFnp));
//	}else	if(verbose_){
//		std::cerr << "Couldn't find " << "TwoBit" << " skipping 2bit file creation" << std::endl;
//	}
	return outputs;
}

njh::sys::RunOutput BioCmdsUtils::runCmdCheck(const std::string & cmd,
		const bfs::path & input, const bfs::path & check) const {
	if (!bfs::exists(check) || njh::files::firstFileIsOlder(check, input)) {
		if (verbose_) {
			std::cout << njh::bashCT::bold << njh::bashCT::blue << "Running: "
					<< njh::bashCT::green << cmd << njh::bashCT::reset << std::endl;
		}
		if(bfs::exists(check)){
			if (verbose_) {
				std::cout << njh::bashCT::bold << njh::bashCT::blue << "\tRemoving out of date file: "
										<< njh::bashCT::green << check  << njh::bashCT::reset << std::endl;
			}
			bfs::remove(check);
		}
		auto runOut = njh::sys::run( { cmd });
		checkRunOutThrow(runOut, __PRETTY_FUNCTION__);
		return runOut;
	} else {
		if (verbose_) {
			std::cout << njh::bashCT::bold << njh::bashCT::blue
					<< "No need to run: " << njh::bashCT::green << cmd << njh::bashCT::reset << std::endl;
		}
	}
	return njh::sys::RunOutput();
}


njh::sys::RunOutput BioCmdsUtils::bowtie2AlignNoSort(const SeqIOOptions & opts,
		const bfs::path & genomeFnp, std::string additionalBowtie2Args) const {

	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bfs::path outputFnp = njh::appendAsNeededRet(opts.out_.outFilename_.string(),
			".bam");
	auto outputFnpTempSam = njh::appendAsNeededRet(opts.out_.outFilename_.string(),
			".sam");
	if (bfs::exists(outputFnpTempSam) && !opts.out_.overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, " << outputFnp
				<< " already exists, use --overWrite to over write it" << "\n";
		throw std::runtime_error { ss.str() };
	}
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

	if(opts.isPairedIn()){
		templateCmd << "bowtie2 -1 " << opts.firstName_ << " -2 " << opts.secondName_ << " -x " << genomePrefix
				<< " " <<  additionalBowtie2Args << " > " << outputFnpTempSam;
	}else{
		templateCmd << "bowtie2 -U " << opts.firstName_ << " -x " << genomePrefix
				<< " " <<  additionalBowtie2Args << " > " << outputFnpTempSam;
	}
	if (verbose_) {
		std::cout << "Running: " << njh::bashCT::green << templateCmd.str()
				<< njh::bashCT::reset << std::endl;
	}
	auto ret = njh::sys::run( { templateCmd.str() });
	BioCmdsUtils::checkRunOutThrow(ret, __PRETTY_FUNCTION__);
	std::stringstream samtoolsCmds;
	samtoolsCmds << "samtools view " << outputFnpTempSam << " -o " << outputFnp;
	if (verbose_) {
		std::cout << "Running: " << njh::bashCT::green << samtoolsCmds.str()
				<< njh::bashCT::reset << std::endl;
	}
	auto samtoolsRunOut = njh::sys::run( { samtoolsCmds.str() });
	BioCmdsUtils::checkRunOutThrow(samtoolsRunOut, __PRETTY_FUNCTION__);
	bfs::remove(outputFnpTempSam);
	return ret;
}

njh::sys::RunOutput BioCmdsUtils::bowtie2Align(const SeqIOOptions & opts,
		const bfs::path & genomeFnp, std::string additionalBowtie2Args) const {

	checkGenomeFnpExistsThrow(genomeFnp, __PRETTY_FUNCTION__);
	bfs::path outputFnp = njh::appendAsNeededRet(opts.out_.outFilename_.string(),
			".sorted.bam");
	auto outputFnpTempSam = njh::appendAsNeededRet(opts.out_.outFilename_.string(),
			".sam");
	if (bfs::exists(outputFnpTempSam) && !opts.out_.overWriteFile_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, " << outputFnp
				<< " already exists, use --overWrite to over write it" << "\n";
		throw std::runtime_error { ss.str() };
	}
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

//	if(opts.isPairedIn()){
//		templateCmd << "bowtie2 -1 " << opts.firstName_ << " -2 " << opts.secondName_ << " -x " << genomePrefix
//				<< " " <<  additionalBowtie2Args << " | samtools view - -b | samtools sort - -o " << outputFnp
//				<< " " << "&& samtools index " << outputFnp;
//	}else{
//		templateCmd << "bowtie2 -U " << opts.firstName_ << " -x " << genomePrefix
//				<< " " <<  additionalBowtie2Args << " | samtools view - -b | samtools sort - -o " << outputFnp
//				<< " " << "&& samtools index " << outputFnp;
//	}
	if(opts.isPairedIn()){
		templateCmd << "bowtie2 -1 " << opts.firstName_ << " -2 " << opts.secondName_ << " -x " << genomePrefix
				<< " " <<  additionalBowtie2Args << " > " << outputFnpTempSam;
	}else{
		templateCmd << "bowtie2 -U " << opts.firstName_ << " -x " << genomePrefix
				<< " " <<  additionalBowtie2Args << " > " << outputFnpTempSam;
	}
	if (verbose_) {
		std::cout << "Running: " << njh::bashCT::green << templateCmd.str()
				<< njh::bashCT::reset << std::endl;
	}
	auto ret = njh::sys::run( { templateCmd.str() });

	BioCmdsUtils::checkRunOutThrow(ret, __PRETTY_FUNCTION__);
	std::stringstream samtoolsCmds;
	samtoolsCmds << "samtools sort " << outputFnpTempSam <<" -o " << outputFnp
					<< " " << "&& samtools index " << outputFnp;
	if (verbose_) {
		std::cout << "Running: " << njh::bashCT::green << samtoolsCmds.str()
				<< njh::bashCT::reset << std::endl;
	}
	auto samtoolsRunOut = njh::sys::run( { samtoolsCmds.str() });
	BioCmdsUtils::checkRunOutThrow(samtoolsRunOut, __PRETTY_FUNCTION__);
	bfs::remove(outputFnpTempSam);
	return ret;
}

njh::sys::RunOutput BioCmdsUtils::lastzAlign(const SeqIOOptions & opts, const LastZPars & pars) const {
	checkGenomeFnpExistsThrow(pars.genomeFnp, __PRETTY_FUNCTION__);
	bfs::path outputFnp = njh::appendAsNeededRet(opts.out_.outFilename_.string(),
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
			<< " --coverage=" << pars.coverage << " --identity=" << pars.identity << " --ambiguous=iupac "
			<< pars.extraLastzArgs << " | samtools view - -b | samtools sort - -o "
			<< outputFnp << " " << "&& samtools index "
			<< outputFnp;
	if (verbose_) {
		std::cout << "Running: " << njh::bashCT::green << lastzCmd.str() << njh::bashCT::reset << std::endl;
	}
	auto ret = njh::sys::run( { lastzCmd.str() });
	BioCmdsUtils::checkRunOutThrow(ret, __PRETTY_FUNCTION__);
	return ret;
}

Json::Value BioCmdsUtils::FastqDumpResults::toJson() const{
	Json::Value ret;

	ret["class"] = njh::typeStr(*this);
	ret["firstMateFnp_"] = njh::json::toJson(firstMateFnp_);
	ret["barcodeFnp_"] = njh::json::toJson(barcodeFnp_);
	ret["secondMateFnp_"] = njh::json::toJson(secondMateFnp_);
	ret["isGzipped_"] = njh::json::toJson(isGzipped_);
	ret["isPairedEnd_"] = njh::json::toJson(isPairedEnd_);
	ret["output_"] = njh::json::toJson(output_);


	return ret;
}
BioCmdsUtils::FastqDumpResults BioCmdsUtils::runFastqDump(const FastqDumpPars & pars) const{
	BioCmdsUtils::FastqDumpResults ret;
	if(!bfs::exists(pars.sraFnp_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << pars.sraFnp_ << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}
	njh::sys::requireExternalProgramThrow(fastqDumpCmd_);
	std::stringstream cmd;
	cmd << fastqDumpCmd_ << " --log-level 0 -X 1 -Z --split-spot " << pars.sraFnp_;
	auto cmdOutput = njh::sys::run({cmd.str()});
	BioCmdsUtils::checkRunOutThrow(cmdOutput, __PRETTY_FUNCTION__);
	//njh::sys::run trim end white space so have to add one;
	uint32_t newLines = countOccurences(cmdOutput.stdOut_, "\n") + 1;

	std::string extraSraArgs = pars.extraSraOptions_;
	if (pars.exportBarCode_) {
		extraSraArgs = njh::replaceString(extraSraArgs, "--skip-technical", "");
	} else if (!njh::containsSubString(pars.extraSraOptions_,
			"--skip-technical")) {
		extraSraArgs += " --skip-technical ";
		//
	}
	if(pars.gzip_){
		extraSraArgs = njh::replaceString(extraSraArgs, "--gzip", "");
	}
	auto outputStub = njh::files::make_path(pars.outputDir_, pars.sraFnp_.filename());

	bfs::path checkFile1        = njh::files::replaceExtension(outputStub, "_1.fastq");
	bfs::path checkFile2        = njh::files::replaceExtension(outputStub, "_2.fastq");
	bfs::path checkFileBarcodes = njh::files::replaceExtension(outputStub, "_barcodes.fastq");
	if(pars.gzip_){
		checkFile1 = checkFile1.string() + ".gz";
		checkFile2 = checkFile2.string() + ".gz";
		checkFileBarcodes = checkFileBarcodes.string() + ".gz";
		ret.isGzipped_ = true;
	}
	bfs::path orgFirstMateFnp  = njh::files::replaceExtension(outputStub, "_1.fastq");
	bfs::path orgBarcodesFnp   = "";
	bfs::path orgSecondMateFnp = "";

	bool needToRun = true;
	ret.firstMateFnp_ = checkFile1;
	//std::cout << "newLines: " << newLines << std::endl;
	if (4 == newLines) {
		ret.isPairedEnd_ = false;
		if(bfs::exists(checkFile1) &&
				njh::files::firstFileIsOlder(pars.sraFnp_, checkFile1)){
			needToRun = false;
		}
	} else if (8 == newLines) {
		ret.secondMateFnp_ = checkFile2;
		ret.isPairedEnd_ = true;
		if(bfs::exists(checkFile1) &&
				bfs::exists(checkFile2) &&
				njh::files::firstFileIsOlder(pars.sraFnp_, checkFile1)&&
				njh::files::firstFileIsOlder(pars.sraFnp_, checkFile2)){
			needToRun = false;
		}
		orgSecondMateFnp = njh::files::replaceExtension(outputStub, "_2.fastq");
	} else if (12 == newLines) {
		ret.secondMateFnp_ = checkFile2;
		if(pars.exportBarCode_){
			ret.barcodeFnp_ = checkFileBarcodes;
			if(bfs::exists(checkFile1) &&
					bfs::exists(checkFile2) &&
					bfs::exists(checkFileBarcodes) &&
					njh::files::firstFileIsOlder(pars.sraFnp_, checkFile1)&&
					njh::files::firstFileIsOlder(pars.sraFnp_, checkFile2)&&
					njh::files::firstFileIsOlder(pars.sraFnp_, checkFileBarcodes)){
				needToRun = false;
			}
		}else{
			if(bfs::exists(checkFile1) &&
					bfs::exists(checkFile2) &&
					njh::files::firstFileIsOlder(pars.sraFnp_, checkFile1)&&
					njh::files::firstFileIsOlder(pars.sraFnp_, checkFile2)){
				needToRun = false;
			}
		}
		orgBarcodesFnp   = njh::files::replaceExtension(outputStub, "_2.fastq");
		orgSecondMateFnp = njh::files::replaceExtension(outputStub, "_3.fastq");
		ret.isPairedEnd_ = true;
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error with " << cmd.str()
				<< " was expecting 4, 8, or 12 lines but got " << newLines << " instead "
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	if(pars.force_){
		needToRun = true;
	}
	if(needToRun){
		auto zipFile = [](IoOptions opts){
			gzZipFile(opts); //zip file from in_ to file in out_
			bfs::remove(opts.in_.inFilename_); //remove original file
		};
		std::stringstream ss;
		ss << fastqDumpCmd_ << " --split-files --defline-seq '@$sn/$ri' --outdir " << pars.outputDir_
				<< " " << extraSraArgs << " " << pars.sraFnp_;
		auto runOutput = njh::sys::run({ss.str()});
		checkRunOutThrow(runOutput, __PRETTY_FUNCTION__);
		if(12 == newLines){
			std::unique_ptr<std::thread> gzBarcodeTh;
			//now to fix the crazy that is the SRA way of dump things with three files
			if(pars.exportBarCode_){
				bfs::path currentBarcodeFnp = njh::files::replaceExtension(outputStub, "_2.fastq");
				if(pars.gzip_){
					IoOptions barIoOpts{InOptions(currentBarcodeFnp), OutOptions(checkFileBarcodes)};
					gzBarcodeTh = std::make_unique<std::thread>(zipFile, barIoOpts);
				}else{
					if(bfs::exists(checkFileBarcodes)){
						bfs::remove(checkFileBarcodes);
					}
					bfs::copy_file(currentBarcodeFnp, checkFileBarcodes);
				}
			}
			IoOptions firstMateIoOpts { InOptions(orgFirstMateFnp), OutOptions(
					checkFile1) };
			firstMateIoOpts.out_.overWriteFile_ = true;
			std::unique_ptr<std::thread> gzFirstMateTh;
			if(pars.gzip_){
				gzFirstMateTh = std::make_unique<std::thread>(zipFile, firstMateIoOpts);
			}
			bfs::path currentSecondMateFnp = njh::files::replaceExtension(outputStub, "_3.fastq");
			auto secondMateIn = SeqIOOptions::genFastqIn(currentSecondMateFnp);
			SeqInput reader{secondMateIn};
			seqInfo seq;
			OutOptions secondMateOutOpts(checkFile2);
			secondMateOutOpts.overWriteFile_ = true;
			OutputStream secondMateOut(secondMateOutOpts);
			reader.openIn();
			while(reader.readNextRead(seq)){
				seq.name_.back() = '2';
				seq.outPutFastq(secondMateOut);
			}
			reader.closeIn();
			bfs::remove(currentSecondMateFnp);
			if(pars.gzip_){
				gzFirstMateTh->join();
				if(pars.exportBarCode_){
					gzBarcodeTh->join();
				}
			}
		} else if(4 == newLines){
			if(pars.gzip_){
				IoOptions firstMateIoOpts { InOptions(orgFirstMateFnp), OutOptions(
						checkFile1) };
				firstMateIoOpts.out_.overWriteFile_ = true;
				zipFile(firstMateIoOpts);
			}
		} else if(8 == newLines){
			//so a couple of samples were found to have a _1 and _3 file for when there were 8 lines, no barcode _2 file present
			bfs::path properSecondMate = njh::files::replaceExtension(outputStub, "_2.fastq");
			bfs::path possibleSecondMate = njh::files::replaceExtension(outputStub, "_3.fastq");
			if(bfs::exists(possibleSecondMate) && !bfs::exists(properSecondMate)){
				IoOptions firstMateIoOpts { InOptions(orgFirstMateFnp), OutOptions(
						checkFile1) };
				firstMateIoOpts.out_.overWriteFile_ = true;
				std::unique_ptr<std::thread> gzFirstMateTh;
				if (pars.gzip_) {
					gzFirstMateTh = std::make_unique<std::thread>(zipFile,
							firstMateIoOpts);
				}
				bfs::path currentSecondMateFnp = njh::files::replaceExtension(
						outputStub, "_3.fastq");
				auto secondMateIn = SeqIOOptions::genFastqIn(currentSecondMateFnp);
				SeqInput reader { secondMateIn };
				seqInfo seq;
				OutOptions secondMateOutOpts(checkFile2);
				secondMateOutOpts.overWriteFile_ = true;
				OutputStream secondMateOut(secondMateOutOpts);
				reader.openIn();
				while (reader.readNextRead(seq)) {
					seq.name_.back() = '2';
					seq.outPutFastq(secondMateOut);
				}
				reader.closeIn();
				bfs::remove(currentSecondMateFnp);
				if (pars.gzip_) {
					gzFirstMateTh->join();
				}
			} else {
				if(pars.gzip_){
					IoOptions firstMateIoOpts { InOptions(orgFirstMateFnp), OutOptions(
							checkFile1) };
					firstMateIoOpts.out_.overWriteFile_ = true;
					IoOptions secondMateIoOpts { InOptions(orgSecondMateFnp), OutOptions(
							checkFile2) };
					secondMateIoOpts.out_.overWriteFile_ = true;
					std::thread firstMateTh(zipFile, firstMateIoOpts);
					std::thread secondMateTh(zipFile, secondMateIoOpts);
					firstMateTh.join();
					secondMateTh.join();
				}
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
	njh::sys::requireExternalProgramThrow(fastqDumpCmd_);
	std::stringstream cmd;
	cmd <<  fastqDumpCmd_ << " --log-level 0 -X 1 -Z --split-spot " << sraFnp;
	auto cmdOutput = njh::sys::run({cmd.str()});
	BioCmdsUtils::checkRunOutThrow(cmdOutput, __PRETTY_FUNCTION__);
	//njh::sys::run trim end white space so have to add one;
	uint32_t newLines = countOccurences(cmdOutput.stdOut_, "\n") + 1;
	if(8 == newLines || 12 == newLines){
		//check for 8 or 12, 12 means triple file _1 forward mate, _2 barcode, _3 reverse mate
		//could also be solved by checking for only 8 and using --skip-technical in the command to skip the barcode file
		return true;
	}
	return false;
}

void BioCmdsUtils::checkRunOutThrow(const njh::sys::RunOutput & runOut,
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
} /* namespace njhseq */
