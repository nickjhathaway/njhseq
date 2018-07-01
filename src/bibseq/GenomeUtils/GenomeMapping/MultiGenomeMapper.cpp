/*
 * MultiGenomeMapper.cpp
 *
 *  Created on: Mar 31, 2017
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "MultiGenomeMapper.hpp"
#include "bibseq/BamToolsUtils.h"




namespace bibseq {

MultiGenomeMapper::MultiGenomeMapper(const inputParameters & pars):pars_(pars){

}

MultiGenomeMapper::MultiGenomeMapper(const bfs::path & genomeDir,
		const std::string & primaryGenome) :
		pars_ { genomeDir, primaryGenome } {

}

MultiGenomeMapper::Genome::Genome(const std::string & name,
		const bfs::path & fnp) :
		name_(name), fnp_(fnp) {
	bib::files::checkExistenceThrow(fnp_, __PRETTY_FUNCTION__);
	fnpTwoBit_ = fnp_;
	fnpTwoBit_.replace_extension(".2bit");
}

void MultiGenomeMapper::Genome::createTwoBit() {
	TwoBit::faToTwoBitPars pars;
	pars.inputFilename = fnp_.string();
	pars.outFilename = fnpTwoBit_.string();
	bool buildTwoBit = false;
	if (!bfs::exists(fnpTwoBit_)) {
		buildTwoBit = true;
	} else if (bib::files::firstFileIsOlder(fnpTwoBit_, fnp_)) {
		buildTwoBit = true;
		pars.overWrite = true;
	}
	if (buildTwoBit) {
		TwoBit::fastasToTwoBit(pars);
	}
	fnpTwoBit_ = pars.outFilename;
	TwoBit::TwoBitFile genFile(fnpTwoBit_);
	chromosomeLengths_ = genFile.getSeqLens();
}

void MultiGenomeMapper::Genome::buildBowtie2Index() const {
	BioCmdsUtils bioRunner;
	bioRunner.RunBowtie2Index(fnp_);
}

Json::Value MultiGenomeMapper::Genome::chromosomeLengths() const {
	Json::Value ret;

	TwoBit::TwoBitFile genFile(fnpTwoBit_);
	auto lens = genFile.getSeqLens();
	auto lenKeys = bib::getVecOfMapKeys(lens);
	bib::sort(lenKeys);
	for (const auto & lenKey : lenKeys) {
		Json::Value lenObj;
		lenObj["name"] = lenKey;
		lenObj["len"] = lens[lenKey];
		ret.append(lenObj);
	}
	return ret;
}


Json::Value MultiGenomeMapper::Genome::toJson() const{
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["name_"] = bib::json::toJson(name_);
	ret["fnp_"] = bib::json::toJson(fnp_);
	ret["fnpTwoBit_"] = bib::json::toJson(fnpTwoBit_);
	ret["gffFnp_"] = bib::json::toJson(gffFnp_);
	ret["chromosomeLengths_"] = bib::json::toJson(chromosomeLengths_);
	return ret;
}

Json::Value MultiGenomeMapper::inputParameters::toJson() const{
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["genomeDir_"] = bib::json::toJson(genomeDir_);
	ret["primaryGenome_"] = bib::json::toJson(primaryGenome_);
	ret["selectedGenomes_"] = bib::json::toJson(selectedGenomes_);
	ret["numThreads_"] = bib::json::toJson(numThreads_);
	ret["workingDirectory_"] = bib::json::toJson(workingDirectory_);
	return ret;
}


Json::Value MultiGenomeMapper::toJson() const{
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["pars_"] = bib::json::toJson(pars_);
	ret["genomes_"] = bib::json::toJson(genomes_);
	return ret;
}


void MultiGenomeMapper::loadGffFnps(){
	loadGffFnps(pars_.gffDir_);
}

void MultiGenomeMapper::loadGffFnps(const bfs::path & gffDir) {
	for (auto & genome : genomes_) {
		bfs::path gffFnp = bib::files::make_path(pars_.gffDir_,
				genome.second->fnp_.filename()).replace_extension("gff");
		if (bfs::exists(gffFnp)) {
			genome.second->gffFnp_ = gffFnp;
		} else {
			gffFnp = gffFnp.string() + "3";
			if (bfs::exists(gffFnp)) {
				genome.second->gffFnp_ = gffFnp;
			}
		}
	}
}


void MultiGenomeMapper::loadInGenomes() {
	std::lock_guard<std::mutex> lock(mut_);
	auto fastaFiles = bib::files::gatherFiles(pars_.genomeDir_, ".fasta", false);
	if (fastaFiles.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": Error, there has to be at least one genome in " << pars_.genomeDir_
				<< "\n";
		ss << "Found none ending with .fasta" << "\n";
		throw std::runtime_error { ss.str() };
	}
	for (const auto & f : fastaFiles) {
		std::string genomeName = bfs::basename(f);
		if (std::string::npos != genomeName.rfind(".")) {
			genomeName = genomeName.substr(0, genomeName.rfind("."));
		}
		if (pars_.selectedGenomes_.empty() || bib::in(genomeName, pars_.selectedGenomes_)) {
			genomes_[genomeName] = std::make_unique<Genome>(genomeName, f);
		}
	}
}

void MultiGenomeMapper::setUpGenomes() {
	std::lock_guard<std::mutex> lock(mut_);
	bib::concurrent::LockableQueue<std::string> queueGenome(getVectorOfMapKeys(genomes_));
	auto setUpGenome = [&queueGenome,this](){
		std::string genome = "";
		BioCmdsUtils bioRunner(pars_.verbose_);
		while(queueGenome.getVal(genome)){
			genomes_.at(genome)->createTwoBit();
			genomes_.at(genome)->buildBowtie2Index();
			//bioRunner.runAllPossibleIndexes(genomes_.at(genome)->fnp_);
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars_.numThreads_; ++t){
		threads.emplace_back(setUpGenome);
	}
	for(auto & t : threads){
		t.join();
	}
}

void MultiGenomeMapper::bioIndexAllGenomes() {
	std::lock_guard<std::mutex> lock(mut_);
	bib::concurrent::LockableQueue<std::string> queueGenome(getVectorOfMapKeys(genomes_));
	auto setUpGenome = [&queueGenome,this](){
		std::string genome = "";
		BioCmdsUtils bioRunner(pars_.verbose_);
		while(queueGenome.getVal(genome)){
			bioRunner.runAllPossibleIndexes(genomes_.at(genome)->fnp_);
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars_.numThreads_; ++t){
		threads.emplace_back(setUpGenome);
	}
	for(auto & t : threads){
		t.join();
	}
}



void MultiGenomeMapper::init() {
	loadInGenomes();
	if ("" != pars_.primaryGenome_ && !bib::in(pars_.primaryGenome_, genomes_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error primary genome: " << pars_.primaryGenome_
				<< " wasn't found in " << pars_.genomeDir_ << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(genomes_), ", ")
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	setUpGenomes();
	loadGffFnps();
}


void MultiGenomeMapper::checkForGenomeThrow(const std::string & genome,
		const std::string & funcName) const{
	if(!hasGenome(genome)){
		std::stringstream ss;
		ss << funcName << ", error don't have genome: " << genome << "\n";
		ss << "Options are: " << bib::conToStr(getVectorOfMapKeys(genomes_), ",") << "\n";
		throw std::runtime_error{ss.str()};
	}
}

bool MultiGenomeMapper::hasGenome(const std::string & genome) const{
	return bib::in(genome, genomes_);
}


std::vector<bfs::path> MultiGenomeMapper::getGenomeFnps() const {
	std::vector<bfs::path> ret;
	for (const auto & genome : genomes_) {
		ret.emplace_back(genome.second->fnp_);
	}
	return ret;
}

void MultiGenomeMapper::setSelectedGenomes(
		const std::string & genomesStr) {
	if("" != genomesStr){
		auto selectedGenomesVec = bib::tokenizeString(genomesStr, ",");
		setSelectedGenomes(selectedGenomesVec);
	}
}

void MultiGenomeMapper::setSelectedGenomes(
		const std::set<std::string> & genomes) {
	pars_.selectedGenomes_ = genomes;
	if(!bib::in(pars_.primaryGenome_, pars_.selectedGenomes_)){
		pars_.selectedGenomes_.emplace(pars_.primaryGenome_);
	}
}

void MultiGenomeMapper::setSelectedGenomes(const VecStr & genomes) {
	setSelectedGenomes(std::set<std::string> { genomes.begin(), genomes.end() });
}

std::unordered_map<std::string, MultiGenomeMapper::AlignCmdOutput> MultiGenomeMapper::alignToGenomes(
		const SeqIOOptions & inputOpts, const bfs::path & outputPrefix) const{
	std::unordered_map<std::string, AlignCmdOutput> ret;
	BioCmdsUtils cmdRunner;
	for (const auto & genome : genomes_) {
		SeqIOOptions opts = inputOpts;
		opts.out_.outFilename_ = outputPrefix.string() + genome.first
				+ "_aligned.sorted.bam";
		AlignCmdOutput output(opts.out_.outFilename_);
		if (bfs::exists(opts.out_.outFilename_)
				&& bib::files::firstFileIsOlder(opts.firstName_, opts.out_.outFilename_)
				&& bib::files::firstFileIsOlder(genome.second->fnp_,
						opts.out_.outFilename_)) {
		} else {
			output.rOutput_ = cmdRunner.bowtie2Align(opts, genome.second->fnp_);
			//std::cout << bib::json::toJson(output.rOutput_) << std::endl;
		}
		ret.emplace(genome.first, output);
	}
	return ret;
}

std::unordered_map<std::string, MultiGenomeMapper::AlignCmdOutput> MultiGenomeMapper::alignToGenomesLastz(
		const SeqIOOptions & inputOpts,
		const bfs::path & outputPrefix,
		const BioCmdsUtils::LastZPars & pars) const {

	bib::sys::requireExternalProgramThrow("samtools");
	bib::sys::requireExternalProgramThrow("lastz");

	std::unordered_map<std::string, AlignCmdOutput> ret;

	std::mutex retMut;
	auto genomeKeys = getVectorOfMapKeys(genomes_);
	bib::concurrent::LockableQueue<std::string> queue(genomeKeys);
	auto alignGenome = [this,&queue, &retMut,&inputOpts,&outputPrefix,&ret,&pars](){
		std::string genome = "";
		BioCmdsUtils cmdRunner;
		while(queue.getVal(genome)){
			SeqIOOptions opts = inputOpts;
			opts.out_.outFilename_ = outputPrefix.string() + genome
					+ "_aligned.sorted.bam";
			AlignCmdOutput output(opts.out_.outFilename_);
			if (bfs::exists(opts.out_.outFilename_)
					&& bib::files::firstFileIsOlder(opts.firstName_, opts.out_.outFilename_)
					&& bib::files::firstFileIsOlder(genomes_.at(genome)->fnp_,
							opts.out_.outFilename_)) {
			} else {
				BioCmdsUtils::LastZPars parsCopy = pars;
				parsCopy.genomeFnp = genomes_.at(genome)->fnpTwoBit_;
				output.rOutput_ = cmdRunner.lastzAlign(opts, parsCopy);
				//output.rOutput_ = cmdRunner.bowtie2Align(opts, genome.second->fnp_);
			}
			{
				std::lock_guard<std::mutex> lock(retMut);
				ret.emplace(genome, output);
			}
		}
	};

	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars_.numThreads_; ++t){
		threads.emplace_back(std::thread(alignGenome));
	}
	for(auto & t : threads){
		t.join();
	}

	return ret;
}



std::unordered_map<std::string, std::vector<GenomicRegion>> MultiGenomeMapper::getRegionsFromBams(
		const std::unordered_map<std::string, bfs::path> & bamFnps) const {
	std::unordered_map<std::string, std::vector<GenomicRegion>> ret;
	for (const auto & bamFnp : bamFnps) {
		if (bib::in(bamFnp.first, genomes_)) {
			BamTools::BamReader bReader;
			bReader.Open(bamFnp.second.string());
			checkBamOpenThrow(bReader, bamFnp.second.string());
			BamTools::BamAlignment bAln;
			std::vector<BamTools::BamAlignment> bAlns;
			auto refIds = bReader.GetReferenceData();
			while (bReader.GetNextAlignment(bAln)) {
				if (bAln.IsMapped()) {
					bAlns.emplace_back(bAln);
				}
			}
			GenomicRegionCounter gCounter;
			gCounter.increaseCounts(bAlns, refIds);
			ret[bamFnp.first] = gCounter.getRegionsLargestOnTop();
		}
	}
	return ret;
}


seqInfo MultiGenomeMapper::extractGenomeRegion(const std::string & genome, const GenomicRegion & region) const{
	checkForGenomeThrow(genome,__PRETTY_FUNCTION__	);
	TwoBit::TwoBitFile refReader(genomes_.at(genome)->fnpTwoBit_);
	std::string refSeq = "";
	refReader[region.chrom_]->getSequence(refSeq, region.start_,
			region.end_, region.reverseSrand_);
	return seqInfo{genome, refSeq};
}


std::vector<seqInfo> MultiGenomeMapper::extractRegions(
		const std::unordered_map<std::string, GenomicRegion> & regions) const {

	std::vector<seqInfo> refSeqs;
	for (const auto & genome : regions) {
		refSeqs.emplace_back(extractGenomeRegion(genome.first, genome.second));
	}
	std::vector<seqInfo> ret;
	for (const auto & ref : refSeqs) {
		bool found = false;
		for (auto & otherRef : ret) {
			if (otherRef.seq_ == ref.seq_) {
				otherRef.name_ += "-" + ref.name_;
				found = true;
				break;
			}
		}
		if (!found) {
			ret.emplace_back(ref);
		}
	}
	return ret;
}

std::unordered_map<std::string, std::vector<seqInfo>> MultiGenomeMapper::getRefSeqsWithPrimaryGenomeAll(const GenomicRegion & region,
		const bfs::path & alignmentsDir,
		const BioCmdsUtils::LastZPars & lzPars) const{
	std::unordered_map<std::string, GenomicRegion> primaryRegion { {
		pars_.primaryGenome_, region } };
	auto primaryRef = extractRegions(primaryRegion);
	auto refFnp = bib::files::make_path(alignmentsDir, "primaryRefSeq.fasta");
	SeqOutput::write(primaryRef, SeqIOOptions::genFastaOut(refFnp));
	/*
	 * 	auto alignOutputs = alignToGenomes(
	 SeqIOOptions::genFastaIn(refFnp),
	 bib::appendAsNeededRet(refAlignsDir.string(), "/"));
	 */
	auto alignOutputs = alignToGenomesLastz(SeqIOOptions::genFastaIn(refFnp),
			bib::appendAsNeededRet(alignmentsDir.string(), "/"),
			lzPars);
	auto bamFnps = MultiGenomeMapper::getBamFnps(alignOutputs);
	auto allRegions = getRegionsFromBams(bamFnps);

	std::unordered_map<std::string, std::vector<seqInfo>>  ret;
	std::mutex retMut;
	bib::concurrent::LockableQueue<std::string> genomesQueue(getVectorOfMapKeys(allRegions));

	auto extractBestGenomeSeq = [this,&genomesQueue,&ret,&retMut, &allRegions,&alignmentsDir](){
		std::string genome = "";
		while(genomesQueue.getVal(genome)){
			if (pars_.primaryGenome_  != genome) {
				TwoBit::TwoBitFile refReader(genomes_.at(genome)->fnpTwoBit_);
				std::vector<seqInfo> genomeSeqs;
				std::string refSeq = "";
				uint32_t extractionCount = 0;
				for(const auto reg : allRegions.at(genome)){
					refReader[reg.chrom_]->getSequence(refSeq, reg.start_,
							reg.end_, reg.reverseSrand_);
					genomeSeqs.emplace_back(genome, refSeq);
					MetaDataInName meta;
					meta.addMeta("genome", genome);
					meta.addMeta("extractionCount", extractionCount);
					meta.addMeta("chrom", reg.chrom_);
					meta.addMeta("end", reg.end_);
					meta.addMeta("start", reg.start_);
					meta.addMeta("regionUid", reg.uid_);
					meta.resetMetaInName(genomeSeqs.back().name_);
					++extractionCount;
				}
				{
					std::lock_guard<std::mutex> lock(retMut);
					for(const auto & genomeSeqPos : iter::range(genomeSeqs.size())){
						ret[genome].emplace_back(genomeSeqs[genomeSeqPos] );
					}
				}
			}
			std::ofstream outBedFile;
			OutOptions outBedOpts(
					bib::files::make_path(alignmentsDir, genome + "_regions.bed"));
			outBedOpts.openFile(outBedFile);
			for (const auto & region : allRegions.at(genome)) {
				outBedFile << region.genBedRecordCore().toDelimStr() << std::endl;
			}
		}
	};

	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars_.numThreads_; ++t){
		threads.emplace_back(std::thread(extractBestGenomeSeq));
	}
	for(auto & t : threads){
		t.join();
	}
	return ret;
}

std::vector<seqInfo> MultiGenomeMapper::getRefSeqsWithPrimaryGenome(const GenomicRegion & region,
		const bfs::path & alignmentsDir,
		const BioCmdsUtils::LastZPars & lzPars,
		bool keepBestOnly) const{
	bool extendAndTrim = false;
	uint32_t extendAndTrimLen = 10;
	aligner alignerObj;
	return getRefSeqsWithPrimaryGenome(region, alignmentsDir, lzPars, keepBestOnly,extendAndTrim, extendAndTrimLen, alignerObj );
}


std::vector<seqInfo> MultiGenomeMapper::getRefSeqsWithPrimaryGenome(
		const GenomicRegion & region, const bfs::path & refAlignsDir,
		const BioCmdsUtils::LastZPars & lzPars,
		bool keepBestOnly,
		bool extendAndTrim,
		uint32_t extendAndTrimLen,
		aligner & orgAlignerObj) const {
	std::unordered_map<std::string, GenomicRegion> primaryRegion { {
		pars_.primaryGenome_, region } };
	auto primaryRef = extractRegions(primaryRegion);
	auto refFnp = bib::files::make_path(refAlignsDir, "primaryRefSeq.fasta");
	SeqOutput::write(primaryRef, SeqIOOptions::genFastaOut(refFnp));
	/*
	 * 	auto alignOutputs = alignToGenomes(
	 SeqIOOptions::genFastaIn(refFnp),
	 bib::appendAsNeededRet(refAlignsDir.string(), "/"));
	 */
	auto alignOutputs = alignToGenomesLastz(SeqIOOptions::genFastaIn(refFnp),
			bib::appendAsNeededRet(refAlignsDir.string(), "/"),
			lzPars);
	auto bamFnps = MultiGenomeMapper::getBamFnps(alignOutputs);
	auto allRegions = getRegionsFromBams(bamFnps);
	std::vector<seqInfo> refSeqs;
	std::vector<aligner> aligners;
	if(pars_.numThreads_ > 1){
		for(uint32_t t = 0; t < pars_.numThreads_; ++t){
			aligners.emplace_back(orgAlignerObj);
		}
	}
	seqInfo primaryRefInfo;
	{
		TwoBit::TwoBitFile refReader(genomes_.at(pars_.primaryGenome_)->fnpTwoBit_);
		std::string refSeq = "";
		refReader[region.chrom_]->getSequence(refSeq, region.start_,
				region.end_, region.reverseSrand_);
		primaryRefInfo = seqInfo(pars_.primaryGenome_, refSeq);
	}
	std::mutex refSeqsMut;
	bib::concurrent::LockableQueue<std::string> genomesQueue(getVectorOfMapKeys(allRegions));
	struct GenExtracRes{
		uint32_t forwardHits_{0};
		uint32_t reverseHits_{0};
		uint32_t extractCounts_{0};
	};
	std::unordered_map<std::string, GenExtracRes> genomeExtractionsResults;
	for(const auto & genome : genomes_){
		genomeExtractionsResults[genome.first] = GenExtracRes{};
	}
	auto extractBestGenomeSeq = [this,&genomesQueue,&refSeqs,&refSeqsMut,&primaryRefInfo,
															 &allRegions, &refAlignsDir,&keepBestOnly,&genomeExtractionsResults,
															 &aligners, &orgAlignerObj,
															 &extendAndTrim, &extendAndTrimLen](uint32_t threadNumber){

		std::string genome = "";
		while(genomesQueue.getVal(genome)){
			auto & regions = allRegions.at(genome);
			if(regions.size() == 0){
				continue;
			}
			if(extendAndTrim){
				for(auto & reg : regions){
					auto extenedRegion = reg;
					extenedRegion.start_ = reg.start_ <= extendAndTrimLen ? 0 : reg.start_ - extendAndTrimLen;
					extenedRegion.end_ = reg.end_ + extendAndTrimLen < genomes_.at(genome)->chromosomeLengths_.at(reg.chrom_) ? reg.end_ + extendAndTrimLen : genomes_.at(genome)->chromosomeLengths_.at(reg.chrom_);
					TwoBit::TwoBitFile tReader(genomes_.at(genome)->fnpTwoBit_);
					auto extractedSeq = extenedRegion.extractSeq(tReader);
					auto trimmedExtractedSeq = extractedSeq;
					readVecTrimmer::GlobalAlnTrimPars trimPars{};
					trimPars.startInclusive_ = 0;
					trimPars.endInclusive_ = len(primaryRefInfo) -1;
					if(pars_.numThreads_ == 1){
						if(len(trimmedExtractedSeq) > orgAlignerObj.parts_.maxSize_){
							orgAlignerObj.parts_.setMaxSize(len(trimmedExtractedSeq));
						}
						readVecTrimmer::trimSeqToRefByGlobalAln(trimmedExtractedSeq, primaryRefInfo, trimPars, orgAlignerObj);
					}else{
						if(len(trimmedExtractedSeq) > aligners[threadNumber].parts_.maxSize_){
							aligners[threadNumber].parts_.setMaxSize(len(trimmedExtractedSeq));
						}
						readVecTrimmer::trimSeqToRefByGlobalAln(trimmedExtractedSeq, primaryRefInfo, trimPars, aligners[threadNumber]);
					}
//					std::cout << trimmedExtractedSeq.name_ << std::endl;
//					std::cout << reg.genBedRecordCore().toDelimStr() << std::endl;
//					std::cout << extenedRegion.genBedRecordCore().toDelimStr() << std::endl;
//					std::cout <<"trimmedExtractedSeq.on_: " <<  bib::colorBool(trimmedExtractedSeq.on_) << std::endl;
//					aligners[threadNumber].alignObjectA_.seqBase_.outPutSeq(std::cout);
//					aligners[threadNumber].alignObjectB_.seqBase_.outPutSeq(std::cout);
					if(trimmedExtractedSeq.on_){
						uint32_t startPos = extractedSeq.seq_.find(trimmedExtractedSeq.seq_);
						uint32_t stopPos = startPos + len(trimmedExtractedSeq);
						uint32_t trimmedOffBack = len(extractedSeq.seq_) - stopPos;
						if(reg.reverseSrand_){
							reg.start_ = extenedRegion.start_ + trimmedOffBack;
							reg.end_ = extenedRegion.end_ - startPos;
						}else{
							reg.start_ = extenedRegion.start_ + startPos;
							reg.end_ = extenedRegion.end_ - trimmedOffBack;
						}
					}
				}
			}

			OutputStream bedOut(OutOptions(bib::files::make_path(refAlignsDir, genome + "_regions.bed")));
			for(const auto & reg : regions){
				++genomeExtractionsResults.at(genome).extractCounts_;
				if(reg.reverseSrand_){
					++genomeExtractionsResults.at(genome).reverseHits_;
				}else{
					++genomeExtractionsResults.at(genome).forwardHits_;
				}
				bedOut << reg.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			}
			TwoBit::TwoBitFile refReader(genomes_.at(genome)->fnpTwoBit_);
			std::string refSeq = "";

			if(regions.size() == 1){
				refReader[regions.front().chrom_]->getSequence(refSeq, regions.front().start_,
						regions.front().end_, regions.front().reverseSrand_);
				MetaDataInName refMeta;
				refMeta.addMeta("genome", genome);
				refMeta.addMeta("chrom", regions.front().chrom_);
				refMeta.addMeta("start", regions.front().start_);
				refMeta.addMeta("end", regions.front().end_);
				refMeta.addMeta("strand", (regions.front().reverseSrand_ ? '-' : '+'));
				{
					std::lock_guard<std::mutex> lock(refSeqsMut);
					refSeqs.emplace_back(seqInfo(genome + " " + refMeta.createMetaName(), refSeq));
				}
			}else{
				uint64_t maxlen = 0;
				for(const auto & reg : regions){
					if(reg.getLen() > maxlen){
						maxlen = reg.getLen();
					}
				}
				readVec::getMaxLength(primaryRefInfo, maxlen);
				aligner alignerObj(maxlen, gapScoringParameters(5,1,5,1,5,1), substituteMatrix(2,-2), true);
				std::vector<std::pair<int32_t,uint32_t>> scores;
				std::vector<seqInfo> genomeSeqs;
									for(const auto  regPos : iter::range(regions.size())){
					const auto & reg  = regions[regPos];
					refReader[reg.chrom_]->getSequence(refSeq, reg.start_,
							reg.end_, reg.reverseSrand_);

					seqInfo genomeSeq(genome, refSeq);
					genomeSeqs.emplace_back(genomeSeq);
					alignerObj.alignCacheGlobal(primaryRefInfo, genomeSeq);
					scores.emplace_back(std::make_pair(alignerObj.parts_.score_, regPos));
										}
				bib::sort(scores, [](const auto & s1, const auto & s2 ){
					return s1.first > s2.first;
				});
				{
					MetaDataInName refMeta;
					refMeta.addMeta("genome", genome);
					refMeta.addMeta("chrom",   regions[scores.front().second].chrom_);
					refMeta.addMeta("start",   regions[scores.front().second].start_);
					refMeta.addMeta("end",     regions[scores.front().second].end_);
					refMeta.addMeta("strand", (regions[scores.front().second].reverseSrand_ ? '-' : '+'));
					genomeSeqs[scores.front().second].name_ += " " + refMeta.createMetaName();
					{
						std::lock_guard<std::mutex> lock(refSeqsMut);
						refSeqs.emplace_back(genomeSeqs[scores.front().second]);
					}
				}
				if(!keepBestOnly){
					for(const auto & scorePos : iter::range<uint32_t>(1, scores.size())){
						MetaDataInName refMeta;
						refMeta.addMeta("genome", genome);
						refMeta.addMeta("chrom",  regions[scores[scorePos].second].chrom_);
						refMeta.addMeta("start",  regions[scores[scorePos].second].start_);
						refMeta.addMeta("end",    regions[scores[scorePos].second].end_);
						refMeta.addMeta("strand", (regions[scores[scorePos].second].reverseSrand_ ? '-' : '+'));
						genomeSeqs[scores[scorePos].second].name_ += "." + estd::to_string(scorePos) +  " " + refMeta.createMetaName();
						{
							std::lock_guard<std::mutex> lock(refSeqsMut);
							refSeqs.emplace_back(genomeSeqs[scores[scorePos].second]);
						}
					}
				}
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars_.numThreads_; ++t){
		threads.emplace_back(std::thread(extractBestGenomeSeq, t));
	}
	bib::concurrent::joinAllThreads(threads);
	table performanceTab(VecStr{"genome", "forwardStrandHits", "reverseStrandHits", "extractionCounts"});
	auto genomeKeys = getVectorOfMapKeys(genomeExtractionsResults);
	bib::sort(genomeKeys);
	for(const auto & genomeKey : genomeKeys){
		performanceTab.addRow(genomeKey,
				genomeExtractionsResults[genomeKey].forwardHits_,
				genomeExtractionsResults[genomeKey].reverseHits_,
				genomeExtractionsResults[genomeKey].extractCounts_);
	}
	auto perTabOpts = TableIOOpts::genTabFileOut(bib::files::make_path(refAlignsDir, "extractionCounts"),true);
	performanceTab.outPutContents(perTabOpts);

	std::vector<seqInfo> ret;
	for (const auto & ref : refSeqs) {
		bool found = false;
		for (auto & otherRef : ret) {
			if (otherRef.seq_ == ref.seq_) {
				otherRef.name_ += "-" + ref.name_;
				found = true;
				break;
			}
		}
		if (!found) {
			ret.emplace_back(ref);
		}
	}

	return ret;
}


std::unordered_map<std::string, bfs::path> MultiGenomeMapper::getBamFnps(
		const std::unordered_map<std::string, MultiGenomeMapper::AlignCmdOutput> & alignOutputs) {
	std::unordered_map<std::string, bfs::path> ret;
	for (const auto & output : alignOutputs) {
		ret[output.first] = output.second.alignedFnp_;
	}
	return ret;
}

}  // namespace bibseq

