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
#include <atomic>

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

njh::sys::RunOutput BioCmdsUtils::lastzAlignNoSort(const SeqIOOptions & opts, const LastZPars & pars) const {
	checkGenomeFnpExistsThrow(pars.genomeFnp, __PRETTY_FUNCTION__);
	bfs::path outputFnp = njh::appendAsNeededRet(opts.out_.outFilename_.string(),
			".bam");
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
			<< pars.extraLastzArgs << " | samtools view - -b -o " << outputFnp;
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

Json::Value BioCmdsUtils::FasterqDumpResults::toJson() const{
	Json::Value ret;

	ret["class"] = njh::typeStr(*this);
	ret["firstMateFnp_"] = njh::json::toJson(firstMateFnp_);
	ret["secondMateFnp_"] = njh::json::toJson(secondMateFnp_);
	ret["isPairedEnd_"] = njh::json::toJson(isPairedEnd_);
	ret["output_"] = njh::json::toJson(output_);


	return ret;
}




BioCmdsUtils::FasterqDumpResults BioCmdsUtils::runFasterqDump(const FasterqDumpPars & inPars) const{
	auto pars = inPars;

	BioCmdsUtils::FasterqDumpResults ret;
	if(!bfs::exists(pars.sraFnp_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error SRA file " << pars.sraFnp_ << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(!njh::endsWith(pars.sraFnp_.string(), ".sra")){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error SRA file " << pars.sraFnp_ << " doesn't end in .sra" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(!bfs::exists(pars.outputDir_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error output directory: " << pars.outputDir_ << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(!bfs::exists(pars.tempDir_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error temp directory: " << pars.tempDir_ << " doesn't exist" << "\n";
		throw std::runtime_error{ss.str()};
	}
	pars.outputDir_ = njh::files::normalize(inPars.outputDir_);
	pars.sraFnp_ = njh::files::normalize(inPars.sraFnp_);
	pars.tempDir_ = njh::files::normalize(inPars.tempDir_);

	njh::sys::requireExternalProgramThrow(fasterqDumpCmd_);
	std::stringstream cmd;
	//apparently fastq-dump might go away in the future which in that case this below won't work anymore
	cmd  << "cd " << pars.sraFnp_.parent_path() << " && " << fastqDumpCmd_ << " --log-level 0 -X 1 -Z --split-spot " << pars.sraFnp_;
	auto cmdOutput = njh::sys::run({cmd.str()});
	BioCmdsUtils::checkRunOutThrow(cmdOutput, __PRETTY_FUNCTION__);
	//njh::sys::run trim end white space so have to add 1;
	uint32_t newLines = countOccurences(cmdOutput.stdOut_, "\n") + 1;
	std::string extraSraArgs = pars.extraSraOptions_;
	if (!njh::containsSubString(pars.extraSraOptions_, "--skip-technical")) {
		extraSraArgs += " --skip-technical ";
	}
	if (!njh::containsSubString(pars.extraSraOptions_, "--force")) {
		extraSraArgs += " --force ";
	}


	auto outputStub = njh::files::make_path(pars.outputDir_, pars.sraFnp_.filename());

	bfs::path checkFile1        = njh::files::replaceExtension(outputStub, "_1.fastq.gz");
	bfs::path checkFile2        = njh::files::replaceExtension(outputStub, "_2.fastq.gz");

	bfs::path orgFirstMateFnp   = njh::files::replaceExtension(outputStub, "_1.fastq");
	bfs::path orgSecondMateFnp  = "";

	bool needToRun = true;
	ret.firstMateFnp_ = checkFile1;
	//std::cout << "newLines: " << newLines << std::endl;
	if (4 == newLines) {
		checkFile1        = njh::files::replaceExtension(outputStub, ".fastq.gz");
		orgFirstMateFnp   = njh::files::replaceExtension(outputStub, ".fastq");

		ret.isPairedEnd_ = false;
		if(bfs::exists(checkFile1) &&
				njh::files::firstFileIsOlder(pars.sraFnp_, checkFile1)){
			needToRun = false;
		}
	} else if (8 == newLines || 12 == newLines) {
		ret.secondMateFnp_ = checkFile2;
		ret.isPairedEnd_ = true;
		if(bfs::exists(checkFile1) &&
				bfs::exists(checkFile2) &&
				njh::files::firstFileIsOlder(pars.sraFnp_, checkFile1)&&
				njh::files::firstFileIsOlder(pars.sraFnp_, checkFile2)){
			needToRun = false;
		}
		orgSecondMateFnp = njh::files::replaceExtension(outputStub, "_2.fastq");
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

		std::stringstream ss;
		//fasterq-dump  --outfile ERR012220  --temp ./
		ss  << "cd " << pars.sraFnp_.parent_path() << " && " << fasterqDumpCmd_ << " --temp " << pars.tempDir_ << " ";
		if(ret.isPairedEnd_){
			ss << " --outfile " << njh::files::replaceExtension(outputStub, "") << " ";
		}else{
			ss << " --outfile " << njh::files::replaceExtension(outputStub, ".fastq") << " ";
		}
		ss << " --threads " << pars.numThreads_ << " ";
		ss << extraSraArgs << " " << bfs::absolute(pars.sraFnp_);
		std::string fasterDumpCmd = ss.str();
		//std::cout << fasterDumpCmd << std::endl;
		auto runOutput = njh::sys::run({fasterDumpCmd});

		checkRunOutThrow(runOutput, __PRETTY_FUNCTION__);
		std::function<void()> readInSeqs;
		struct SeqCache{
			std::vector<seqInfo> seqs_;

			std::atomic<uint32_t> count_{0};

			std::mutex mut_;
			void addSeq(const seqInfo & seq){
				std::lock_guard<std::mutex> lock(mut_);
				seqs_.emplace_back(seq);
				count_++;
			}
			void addSeqs(const std::vector<seqInfo> & seqs){
				std::lock_guard<std::mutex> lock(mut_);
				std::move(seqs.begin(), seqs.end(), std::back_inserter(seqs_));
				count_+= seqs.size();
			}

			std::vector<seqInfo> getSeqs(){
				std::vector<seqInfo> writeSeqs;
				//lock the seqs_ vector and move it to the writeSews
				std::lock_guard<std::mutex> lock(mut_);
				//std::cout << "r1_cache.seqs_: " << r1_cache.seqs_.size() << std::endl;
				writeSeqs = std::move(seqs_);
				seqs_.clear();
				count_ = 0;
				return writeSeqs;
			}
			std::atomic<bool> done_{false};

			uint32_t maxCacheSize_{1000000}; //! maxCacheSize_ needs to be greater than writeCacheSize_ will there will be a forever loop
			uint32_t writeCacheSize_{50};
			uint32_t subCacheReadSize_{100};
			std::chrono::microseconds writeCheckWaitTime{std::chrono::milliseconds(1)};

		};

		if(ret.isPairedEnd_){
			auto r1_inputSeqOpts = SeqIOOptions::genFastqIn(orgFirstMateFnp);
			r1_inputSeqOpts.includeWhiteSpaceInName_ = false;
			auto r1_outputSeqOpts = SeqIOOptions::genFastqOutGz(checkFile1);
			r1_outputSeqOpts.out_.overWriteFile_ = true;

			auto r2_inputSeqOpts = SeqIOOptions::genFastqIn(orgSecondMateFnp);
			r2_inputSeqOpts.includeWhiteSpaceInName_ = false;
			auto r2_outputSeqOpts = SeqIOOptions::genFastqOutGz(checkFile2);
			r2_outputSeqOpts.out_.overWriteFile_ = true;

			if (pars.numThreads_ > 1) {
				SeqCache r1_cache;
				r1_cache.writeCacheSize_ = pars.writeCacheSize;
				r1_cache.writeCheckWaitTime = pars.writeWaitTime;
				r1_cache.subCacheReadSize_ = pars.readSubCacheSize_;

				SeqCache r2_cache;
				r2_cache.writeCacheSize_ = pars.writeCacheSize;
				r2_cache.writeCheckWaitTime = pars.writeWaitTime;
				r2_cache.subCacheReadSize_ = pars.readSubCacheSize_;

				auto r1_readInSeqsTest= [&r1_inputSeqOpts,&r1_cache](){
					SeqInput reader(r1_inputSeqOpts);
					reader.openIn();
					seqInfo seq;
					std::vector<seqInfo> subCache;
					while(reader.readNextRead(seq)){
						while(r1_cache.count_.load() > r1_cache.maxCacheSize_){
							std::this_thread::sleep_for(r1_cache.writeCheckWaitTime);
						}
						subCache.emplace_back(seq);
						if(subCache.size() > r1_cache.subCacheReadSize_){
							r1_cache.addSeqs(subCache);
							subCache.clear();
						}
					}
					if(!subCache.empty()){
						r1_cache.addSeqs(subCache);
						subCache.clear();
					}
					r1_cache.done_ = true;
				};

				auto r1_writeSeqsTest= [&r1_outputSeqOpts,&r1_cache](){
					SeqOutput writer(r1_outputSeqOpts);
					writer.openOut();
					//keep attempting to write if the cache is either not done being filled or it has left over sequence
					while(!r1_cache.done_.load() || r1_cache.count_.load() > 0){
						std::this_thread::sleep_for(r1_cache.writeCheckWaitTime);
						//write up cache reaches write limit or it's done being filled
						if(r1_cache.count_.load() > r1_cache.writeCacheSize_ || r1_cache.done_.load()){
							std::vector<seqInfo> writeSeqs;
							{
								//lock the seqs_ vector and move it to the writeSews
								std::lock_guard<std::mutex> lock(r1_cache.mut_);
								//std::cout << "r1_cache.seqs_: " << r1_cache.seqs_.size() << std::endl;
								writeSeqs = std::move(r1_cache.seqs_);
								r1_cache.seqs_.clear();
								r1_cache.count_ = 0;
							}
							writer.write(writeSeqs);
						}
					}
				};

				auto r2_readInSeqsTest= [&r2_inputSeqOpts,&r2_cache](){
					SeqInput reader(r2_inputSeqOpts);
					reader.openIn();
					seqInfo seq;
					std::vector<seqInfo> subCache;
					while(reader.readNextRead(seq)){
						while(r2_cache.count_.load() > r2_cache.maxCacheSize_){
							std::this_thread::sleep_for(r2_cache.writeCheckWaitTime);
						}
						subCache.emplace_back(seq);
						if(subCache.size() > r2_cache.subCacheReadSize_){
							r2_cache.addSeqs(subCache);
							subCache.clear();
						}
					}
					if(!subCache.empty()){
						r2_cache.addSeqs(subCache);
						subCache.clear();
					}
					r2_cache.done_ = true;
				};

				auto r2_writeSeqsTest= [&r2_outputSeqOpts,&r2_cache](){
					SeqOutput writer(r2_outputSeqOpts);
					writer.openOut();
					//keep attempting to write if the cache is either not done being filled or it has left over sequence
					while(!r2_cache.done_.load() || r2_cache.count_.load() > 0){
						std::this_thread::sleep_for(r2_cache.writeCheckWaitTime);
						//write up cache reaches write limit or it's done being filled
						if(r2_cache.count_.load() > r2_cache.writeCacheSize_ || r2_cache.done_.load()){
							std::vector<seqInfo> writeSeqs;
							{
								//lock the seqs_ vector and move it to the writeSews
								std::lock_guard<std::mutex> lock(r2_cache.mut_);
								writeSeqs = std::move(r2_cache.seqs_);
								r2_cache.seqs_.clear();
								r2_cache.count_ = 0;
							}
							writer.write(writeSeqs);
						}
					}
				};
				if (pars.numThreads_ > 3) {
					std::thread r1_readingThread(r1_readInSeqsTest);
					std::thread r1_writingThread(r1_writeSeqsTest);
					std::thread r2_readingThread(r2_readInSeqsTest);
					std::thread r2_writingThread(r2_writeSeqsTest);
					r1_readingThread.join();
					r1_writingThread.join();
					r2_readingThread.join();
					r2_writingThread.join();
				} else {
					std::thread r1_readingThread(r1_readInSeqsTest);
					std::thread r1_writingThread(r1_writeSeqsTest);
					r1_readingThread.join();
					r1_writingThread.join();
					std::thread r2_readingThread(r2_readInSeqsTest);
					std::thread r2_writingThread(r2_writeSeqsTest);
					r2_readingThread.join();
					r2_writingThread.join();
				}
			} else {
				{
					SeqInput reader(r1_inputSeqOpts);
					SeqOutput writer(r1_outputSeqOpts);
					reader.openIn();
					writer.openOut();
					seqInfo seq;
					while(reader.readNextRead(seq)){
						writer.write(seq);
					}
				}
				{
					SeqInput reader(r2_inputSeqOpts);
					SeqOutput writer(r2_outputSeqOpts);
					reader.openIn();
					writer.openOut();
					seqInfo seq;
					while(reader.readNextRead(seq)){
						writer.write(seq);
					}
				}
			}
			bfs::remove(orgFirstMateFnp);
			bfs::remove(orgSecondMateFnp);



		} else {
			auto inputSeqOpts = SeqIOOptions::genFastqIn(orgFirstMateFnp);
			inputSeqOpts.includeWhiteSpaceInName_ = false;
			auto outputSeqOpts = SeqIOOptions::genFastqOutGz(checkFile1);
			outputSeqOpts.out_.overWriteFile_ = true;
			if(pars.numThreads_ > 1){


				SeqCache cache;
				cache.writeCacheSize_ = pars.writeCacheSize;
				cache.writeCheckWaitTime = pars.writeWaitTime;
				cache.subCacheReadSize_ = pars.readSubCacheSize_;

				auto readInSeqsTest= [&inputSeqOpts,&cache](){
					SeqInput reader(inputSeqOpts);
					reader.openIn();
					seqInfo seq;
					std::vector<seqInfo> subCache;
					while(reader.readNextRead(seq)){
						while(cache.count_.load() > cache.maxCacheSize_){
							std::this_thread::sleep_for(cache.writeCheckWaitTime);
						}
						subCache.emplace_back(seq);
						if(subCache.size() > cache.subCacheReadSize_){
							cache.addSeqs(subCache);
							subCache.clear();
						}
					}
					if(!subCache.empty()){
						cache.addSeqs(subCache);
						subCache.clear();
					}
					cache.done_ = true;
				};

				auto writeSeqsTest= [&outputSeqOpts,&cache](){
					SeqOutput writer(outputSeqOpts);
					writer.openOut();
					//keep attempting to write if the cache is either not done being filled or it has left over sequence
					while(!cache.done_.load() || cache.count_.load() > 0){
						std::this_thread::sleep_for(cache.writeCheckWaitTime);
						//write up cache reaches write limit or it's done being filled
						//std::cout << cache.count_.load() << std::endl;
						if(cache.count_.load() > cache.writeCacheSize_ || cache.done_.load()){
							std::vector<seqInfo> writeSeqs;
							{
								//lock the seqs_ vector and move it to the writeSews
								std::lock_guard<std::mutex> lock(cache.mut_);
								writeSeqs = std::move(cache.seqs_);

								cache.seqs_.clear();
								cache.count_ = 0;
							}
							writer.write(writeSeqs);
						}
					}
				};
				std::thread readingThread(readInSeqsTest);
				std::thread writingThread(writeSeqsTest);
				readingThread.join();
				writingThread.join();

			} else {
				SeqInput reader(inputSeqOpts);
				SeqOutput writer(outputSeqOpts);
				reader.openIn();
				writer.openOut();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					writer.write(seq);
				}
			}
			bfs::remove(orgFirstMateFnp);
		}
	}
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
	cmd  << "fastq-dump --log-level 0 -X 1 -Z --split-spot " << pars.sraFnp_;
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
