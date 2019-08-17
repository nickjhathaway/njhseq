/*
 * TranslatorByAlignment.cpp
 *
 *  Created on: Aug 14, 2019
 *      Author: nicholashathaway
 */


#include "TranslatorByAlignment.hpp"


namespace njhseq {



TranslatorByAlignment::TranslatorByAlignmentPars::TranslatorByAlignmentPars(){
	lzPars_.coverage = 100;
	lzPars_.identity = 70;
}


TranslatorByAlignment::VariantsInfo::VariantsInfo(const std::string & id) : id_(id){

}


void TranslatorByAlignment::VariantsInfo::setFinals(const RunPars & rPars, uint32_t totalPopCount){
	//filter saps and indels by occurrence cut off
	for(auto & snp : snps){
		for(const auto & b : snp.second){
			if(b.second < rPars.occurrenceCutOff){
				continue;
			}
			snpsFinal[snp.first][b.first] = b.second;
			if(b.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
				variablePositons_.emplace(snp.first);
			}
		}
	}
	for(const auto & del : deletions){
		for(const auto & d : del.second){
			if(d.second < rPars.occurrenceCutOff){
				continue;
			}
			deletionsFinal[del.first][d.first] = d.second;
			if(d.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
				variablePositons_.emplace(del.first);
			}
		}
	}
	for(const auto & ins : insertions){
		for(const auto & i : ins.second){
			if(i.second < rPars.occurrenceCutOff){
				continue;
			}
			insertionsFinal[ins.first][i.first] = i.second;
			if(i.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
				variablePositons_.emplace(ins.first);
			}
		}
	}
}


Bed3RecordCore TranslatorByAlignment::VariantsInfo::getVariableRegion() {
	if (!variablePositons_.empty()) {
		return Bed3RecordCore(id_,
				*std::min_element(variablePositons_.begin(), variablePositons_.end()),
				*std::max_element(variablePositons_.begin(), variablePositons_.end()));
	}
	return Bed3RecordCore(id_, std::numeric_limits<uint32_t>::max(),
			std::numeric_limits<uint32_t>::max());
}


void TranslatorByAlignment::VariantsInfo::addVariantInfo(
		const std::string & alignedRefSeq,
		const std::string & alignedQuerySeq,
		uint32_t querySeqCount,
		const comparison & comp,
		uint32_t offSetStart
		){
	uint32_t queryAlnStart = alignedQuerySeq.find_first_not_of("-");
	uint32_t queryAlnEnd = alignedQuerySeq.find_last_not_of("-");
	for(const auto & seqPos : iter::range(queryAlnStart, queryAlnEnd + 1)){
		if('-' != alignedRefSeq[seqPos]){
			uint32_t seqChromPosition = getRealPosForAlnPos(alignedRefSeq, seqPos) + offSetStart;
			allBases[seqChromPosition][alignedQuerySeq[seqPos]] += querySeqCount;
		}
	}
	for(const auto & m : comp.distances_.mismatches_){
		snps[m.second.refBasePos + offSetStart][m.second.seqBase]+= querySeqCount;
	}
	for(const auto & gap : comp.distances_.alignmentGaps_){
		if(gap.second.ref_){
			//insertion
			insertions[gap.second.refPos_ + offSetStart][gap.second.gapedSequence_]+=querySeqCount;
		}else{
			//deletion
			deletions[gap.second.refPos_ + offSetStart][gap.second.gapedSequence_]+=querySeqCount;
		}
	}
}


std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes> TranslatorByAlignment::translateBasedOnAlignment(
			const ReAlignedSeq & realigned,
			const GeneFromGffs & currentGene,
			const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
			aligner & alignerObj){
	std::unordered_map<std::string, TranslateSeqRes> ret;
	for(const auto & transcript : currentGene.mRNAs_){
		auto currentTranscriptInfo = transcriptInfosForGene.at(transcript->getIDAttr());
		auto genePosInfoByGDna = currentTranscriptInfo->getInfosByGDNAPos();
		bool endsAtStopCodon = false;
		uint32_t transStart = 0;
		seqInfo balnSeq(realigned.querySeq_.name_);
		std::vector<uint32_t> codons;
		std::vector<GFFCore> cDNAIntersectedWith;
		for (const auto & cDna : currentGene.CDS_.at(transcript->getIDAttr())) {
			if (realigned.gRegion_.overlaps(*cDna, 3)) {
				cDNAIntersectedWith.emplace_back(*cDna);
			}
		}
		if(cDNAIntersectedWith.size() == 0){

		} else {
			if (cDNAIntersectedWith.size() == 1
					&& realigned.gRegion_.start_ >= cDNAIntersectedWith.front().start_ - 1
					&& realigned.gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
				balnSeq = realigned.querySeq_;
				if (currentGene.gene_->isReverseStrand()) {
					if (genePosInfoByGDna.at(realigned.gRegion_.start_).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
					uint32_t gPos = realigned.gRegion_.end_ - 1;
					auto codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				} else {
					if (genePosInfoByGDna.at(realigned.gRegion_.end_ - 1).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
					uint32_t gPos = realigned.gRegion_.start_;
					uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						++gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				}
			} else {
				njh::sort(cDNAIntersectedWith,
						[](const GenomicRegion & reg1, const GenomicRegion & reg2) {
							if(reg1.start_ < reg2.start_) {
								return true;
							}
							return false;
						});
				if (currentGene.gene_->isReverseStrand()) {
					auto cDnaStop = cDNAIntersectedWith.back().end_;
					uint32_t gPos = std::min(cDnaStop, realigned.gRegion_.end_) - 1;
					auto codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				} else {
					auto cDnaStart = cDNAIntersectedWith.front().start_ - 1;
					uint32_t gPos = std::max(cDnaStart, realigned.gRegion_.start_);
					uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						++gPos;
						codon = njh::mapAt(genePosInfoByGDna, gPos).codonPos_;
						++transStart;
					}
				}
				std::vector<uint32_t> starts;
				std::vector<uint32_t> ends;
				for (const auto & cDna : cDNAIntersectedWith) {
					auto cDnaStart = cDna.start_ - 1;
					auto detStart = std::max(cDnaStart, realigned.gRegion_.start_);
					auto detStop = std::min(cDna.end_, realigned.gRegion_.end_);
					ends.emplace_back(detStop);
					starts.emplace_back(detStart);
					detStart -= realigned.gRegion_.start_;
					detStop -= realigned.gRegion_.start_;
					auto alnStart = getAlnPosForRealPos(realigned.alnRefSeq_.seq_,detStart);
					auto alnStop = getAlnPosForRealPos(realigned.alnRefSeq_.seq_, detStop - 1);
					balnSeq.append(
							realigned.alnQuerySeq_.getSubRead(alnStart,
									alnStop - alnStart + 1));
				}
				uint32_t cDnaStart = *std::min_element(starts.begin(), starts.end());
				uint32_t cDnaStop = *std::max_element(ends.begin(), ends.end());
				if (currentGene.gene_->isReverseStrand()) {
					if (genePosInfoByGDna.at(cDnaStart).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				} else {
					if (genePosInfoByGDna.at(cDnaStop - 1).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				}
				balnSeq.removeGaps();
			}
			if (currentGene.gene_->isReverseStrand()) {
				balnSeq.reverseComplementRead(false, true);
			}
			auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);
			MetaDataInName transMeta;
			transMeta.addMeta("transcript", transcript->getIDAttr());
			balnSeqTrans.name_ += transMeta.createMetaName();

			alignerObj.alignCacheGlobal(currentTranscriptInfo->protein_, balnSeqTrans);
			alignerObj.profilePrimerAlignment(currentTranscriptInfo->protein_, balnSeqTrans);
			TranslateSeqRes tRes;

			tRes.translation_ = balnSeqTrans;
			tRes.refAlnTranslation_ = alignerObj.alignObjectA_.seqBase_;
			tRes.queryAlnTranslation_ = alignerObj.alignObjectB_.seqBase_;
			tRes.comp_ = alignerObj.comp_;

			ret[transcript->getIDAttr()] = tRes;
		}
	}
	return ret;
}


std::unordered_map<std::string, TranslatorByAlignment::TranslateSeqRes> TranslatorByAlignment::translateBasedOnAlignment(
		const BamTools::BamAlignment & bAln,
		const GeneFromGffs & currentGene,
		const std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> & transcriptInfosForGene,
		TwoBit::TwoBitFile & tReader,
		aligner & alignerObj,
		const BamTools::RefVector & refData){

	std::unordered_map<std::string, TranslateSeqRes> ret;

	auto results = std::make_shared<AlignmentResults>(bAln, refData, true);
	results->setRefSeq(tReader);
	results->setComparison(true);


	for(const auto & transcript : currentGene.mRNAs_){
		auto currentTranscriptInfo = transcriptInfosForGene.at(transcript->getIDAttr());
		auto genePosInfoByGDna = currentTranscriptInfo->getInfosByGDNAPos();
		bool endsAtStopCodon = false;
		uint32_t transStart = 0;
		seqInfo balnSeq(bAln.Name);
		std::vector<uint32_t> codons;
		std::vector<GFFCore> cDNAIntersectedWith;
		for (const auto & cDna : currentGene.CDS_.at(transcript->getIDAttr())) {
			if (results->gRegion_.overlaps(*cDna)) {
				cDNAIntersectedWith.emplace_back(*cDna);
			}
		}
		if(cDNAIntersectedWith.size() == 0){

		} else {
			if (cDNAIntersectedWith.size() == 1
					&& results->gRegion_.start_ >= cDNAIntersectedWith.front().start_ - 1
					&& results->gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
				balnSeq = *(results->alnSeq_);
				if (currentGene.gene_->isReverseStrand()) {
					if (genePosInfoByGDna.at(results->gRegion_.start_).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
					uint32_t gPos = results->gRegion_.end_ - 1;
					auto codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				} else {
					if (genePosInfoByGDna.at(results->gRegion_.end_ - 1).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
					uint32_t gPos = results->gRegion_.start_;
					uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						++gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				}
			} else {
				njh::sort(cDNAIntersectedWith,
						[](const GenomicRegion & reg1, const GenomicRegion & reg2) {
							if(reg1.start_ < reg2.start_) {
								return true;
							}
							return false;
						});

				if (currentGene.gene_->isReverseStrand()) {
					auto cDnaStop = cDNAIntersectedWith.back().end_;
					uint32_t gPos = std::min(cDnaStop, results->gRegion_.end_) - 1;
					auto codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						--gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				} else {
					auto cDnaStart = cDNAIntersectedWith.front().start_ - 1;
					uint32_t gPos = std::max(cDnaStart, results->gRegion_.start_);
					uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
					while (0 != codon) {
						++gPos;
						codon = genePosInfoByGDna.at(gPos).codonPos_;
						++transStart;
					}
				}
				std::vector<uint32_t> starts;
				std::vector<uint32_t> ends;
				for (const auto & cDna : cDNAIntersectedWith) {
					auto cDnaStart = cDna.start_ - 1;
					auto detStart = std::max(cDnaStart, results->gRegion_.start_);
					auto detStop = std::min(cDna.end_, results->gRegion_.end_);
					ends.emplace_back(detStop);
					starts.emplace_back(detStart);
					detStart -= results->gRegion_.start_;
					detStop -= results->gRegion_.start_;
					auto alnStart = getAlnPosForRealPos(results->refSeqAligned_->seq_,
							detStart);
					auto alnStop = getAlnPosForRealPos(results->refSeqAligned_->seq_,
							detStop - 1);
					balnSeq.append(
							results->alnSeqAligned_->getSubRead(alnStart,
									alnStop - alnStart + 1));
				}
				uint32_t cDnaStart = *std::min_element(starts.begin(), starts.end());
				uint32_t cDnaStop = *std::max_element(ends.begin(), ends.end());
				if (currentGene.gene_->isReverseStrand()) {
					if (genePosInfoByGDna.at(cDnaStart).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				} else {
					if (genePosInfoByGDna.at(cDnaStop - 1).cDNAPos_
							== currentTranscriptInfo->cDna_.seq_.size() - 1) {
						endsAtStopCodon = true;
					}
				}
				balnSeq.removeGaps();
			}
			if (currentGene.gene_->isReverseStrand()) {
				balnSeq.reverseComplementRead(false, true);
			}

			auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);
			MetaDataInName transMeta;
			transMeta.addMeta("transcript", transcript->getIDAttr());
			balnSeqTrans.name_ += transMeta.createMetaName();

			alignerObj.alignCacheGlobal(currentTranscriptInfo->protein_, balnSeqTrans);
			alignerObj.profilePrimerAlignment(currentTranscriptInfo->protein_, balnSeqTrans);
			TranslateSeqRes tRes;

			tRes.translation_ = balnSeqTrans;
			tRes.refAlnTranslation_ = alignerObj.alignObjectA_.seqBase_;
			tRes.queryAlnTranslation_ = alignerObj.alignObjectB_.seqBase_;
			tRes.comp_ = alignerObj.comp_;

			ret[transcript->getIDAttr()] = tRes;
		}
	}
	return ret;
}


TranslatorByAlignment::TranslatorByAlignment(const TranslatorByAlignmentPars & pars): pars_(pars){
	njh::sys::requireExternalProgramThrow("samtools");
	if(!pars_.useLastz_){
		njh::sys::requireExternalProgramThrow("bowtie2");
	}else{
		njh::sys::requireExternalProgramThrow("lastz");
	}
}


TranslatorByAlignment::TranslatorByAlignmentResult TranslatorByAlignment::run(const SeqIOOptions & seqOpts,
		const std::unordered_map<std::string, uint32_t> & sampCountsForHaps,
		const RunPars & rPars){

	uint32_t totalPopCount = 0;
	for(const auto & sampCount : sampCountsForHaps){
		totalPopCount += sampCount.second;
	}
	TranslatorByAlignmentResult ret;
	std::vector<bfs::path> fnpsToRemove;
	auto seqInputFnp = njh::files::make_path(pars_.workingDirtory_, "inputSeqs.fasta");
	fnpsToRemove.emplace_back(seqInputFnp);
	uint64_t seqMaxLen = 0;

	VecStr names;
	{
		//write out fasta file of input
		seqInfo seq;
		SeqInput reader(seqOpts);
		reader.openIn();
		SeqOutput writer(SeqIOOptions::genFastaOut(seqInputFnp));
		writer.openOut();
		uint32_t pos = 0;
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, seqMaxLen);
			names.emplace_back(seq.name_);
			seq.name_ = estd::to_string(pos);
			++pos;
			writer.write(seq);
		}
	}

//	for(const auto & pos : iter::range(names.size())){
//		std::cout << pos << ":" << names[pos] << std::endl;
//	}

	BioCmdsUtils bRunner(false);
	bRunner.RunFaToTwoBit(pars_.lzPars_.genomeFnp);
	if(!pars_.useLastz_){
		bRunner.RunBowtie2Index(pars_.lzPars_.genomeFnp);
	}

	auto gprefix = bfs::path(pars_.lzPars_.genomeFnp).replace_extension("");
	auto twoBitFnp = gprefix.string() + ".2bit";

	TwoBit::TwoBitFile tReader(twoBitFnp);

	auto uniqueSeqInOpts = SeqIOOptions::genFastaIn(seqInputFnp);
	uniqueSeqInOpts.out_.outFilename_ = njh::files::make_path(pars_.workingDirtory_, "aligned_inputSeqs.sorted.bam");
	uniqueSeqInOpts.out_.outExtention_ = ".sorted.bam";
	fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_);
	fnpsToRemove.emplace_back(uniqueSeqInOpts.out_.outFilename_.string() + ".bai");

	uniqueSeqInOpts.out_.transferOverwriteOpts(seqOpts.out_);
	if(!pars_.useLastz_){
		auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, pars_.lzPars_.genomeFnp, pars_.additionalBowtieArguments_);

		//auto bowtieRunOut = bRunner.bowtie2Align(uniqueSeqInOpts, pars_.lzPars_.genomeFnp, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 --end-to-end");
		BioCmdsUtils::checkRunOutThrow(bowtieRunOut, __PRETTY_FUNCTION__);
	}else{
		auto lastzRunOut = bRunner.lastzAlign(uniqueSeqInOpts, pars_.lzPars_);
		BioCmdsUtils::checkRunOutThrow(lastzRunOut, __PRETTY_FUNCTION__);
	}

	auto regionsCounter = GenomicRegionCounter::countRegionsInBam(uniqueSeqInOpts.out_.outName());
	auto ids = regionsCounter.getIntersectingGffIds(pars_.gffFnp_);
	ret.geneIds_ = ids;

	// get gene information
	auto geneInfoDir = njh::files::make_path(pars_.workingDirtory_, "geneInfos");
	if(pars_.keepTemporaryFiles_){
		njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});
	}
	OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));

	std::unordered_map<std::string, VecStr> idToTranscriptName;
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes = GeneFromGffs::getGenesFromGffForIds(pars_.gffFnp_, ids);
	//std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>>> geneTranscriptInfos;
	;
	uint64_t proteinMaxLen = 0;
	std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;

	for(const auto & gene : genes){
		genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
		for(const auto & transcript : gene.second->mRNAs_){
			idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
		}
		ret.transcriptInfosForGene_[gene.first] = gene.second->generateGeneSeqInfo(tReader, false);
		for(const auto & transcriptInfo : ret.transcriptInfosForGene_[gene.first]){
			readVec::getMaxLength(transcriptInfo.second->protein_, proteinMaxLen);
		}
		if(pars_.keepTemporaryFiles_){
			gene.second->writeOutGeneInfo(tReader, outOpts);
		}
	}

	aligner alignObj(proteinMaxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(1,-1));
	aligner alignObjSeq(seqMaxLen + rPars.realnPars.extendAmount * 2, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(1,-1));

	std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> regionsToGeneIds;
	//targetName, GeneID, AA Position
	std::unordered_map<std::string, std::vector<std::string>> alnRegionToGeneIds;
	for(const auto & gCount : regionsCounter.counts_){
		for (const auto & g : genesByChrom[gCount.second.region_.chrom_]) {
			for(const auto & t : g.second->mRNAs_){
				regionsToGeneIds[gCount.second.region_.createUidFromCoords()][g.first].emplace(t->getIDAttr());
				alnRegionToGeneIds[gCount.second.region_.createUidFromCoords()].emplace_back(g.first);
			}
		}
	}
	BamTools::BamReader bReader;
	bReader.Open(uniqueSeqInOpts.out_.outName().string());
	checkBamOpenThrow(bReader, uniqueSeqInOpts.out_.outName());
	auto refData = bReader.GetReferenceData();
	BamTools::BamAlignment bAln;
	auto chromLengths = tReader.getSeqLens();

	while (bReader.GetNextAlignment(bAln)) {
		if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
			bAln.Name = names[njh::StrToNumConverter::stoToNum<uint32_t>(bAln.Name)];
			auto balnGenomicRegion = GenomicRegion(bAln, refData);
			auto results = ReAlignedSeq::genRealignment(bAln, refData, alignObjSeq, chromLengths, tReader, rPars.realnPars);
			ret.seqAlns_[bAln.Name].emplace_back(results);
			if (!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)) {
				continue;
			}
			for (const auto & g : alnRegionToGeneIds.at(balnGenomicRegion.createUidFromCoords())) {
				const auto & currentGene = genes.at(g);
				const auto & currentGeneInfo = ret.transcriptInfosForGene_.at(g);
				auto translations = translateBasedOnAlignment(results, *currentGene, currentGeneInfo, alignObj);
				for(const auto & trans : translations){
					ret.translations_[bAln.Name].emplace(trans);
				}
			}
		}
	}
	if(!pars_.keepTemporaryFiles_){
		for(const auto & fnp : fnpsToRemove){
			if(bfs::exists(fnp)){
				bfs::remove(fnp);
			}
		}
	}
	//index snps
	for(const auto & seqName : ret.seqAlns_){
		for(const auto & aln : seqName.second){
			if(!njh::in(aln.gRegion_.chrom_, ret.seqVariants_)){
				ret.seqVariants_.emplace(aln.gRegion_.chrom_, VariantsInfo(aln.gRegion_.chrom_));
			}
			for(uint32_t seqPos = 0; seqPos < aln.refSeq_.seq_.size(); ++ seqPos){
				ret.baseForPosition_[aln.gRegion_.chrom_][aln.gRegion_.start_ + seqPos] = aln.refSeq_.seq_[seqPos];
			}
			uint32_t popCount = njh::mapAt(sampCountsForHaps, aln.querySeq_.name_);
			ret.seqVariants_.at(aln.gRegion_.chrom_).addVariantInfo(
					aln.alnRefSeq_.seq_, aln.alnQuerySeq_.seq_,
					popCount, aln.comp_,
					aln.gRegion_.start_);
		}
	}

	for(auto & varPerChrom : ret.seqVariants_){
		varPerChrom.second.setFinals(rPars, totalPopCount);
	}
	//index amino acid changes per transcript
	for(const auto & seqName : ret.translations_){
		for(const auto & transcript : seqName.second){
			if("" == ret.proteinForTranscript_[transcript.first]){
				ret.proteinForTranscript_[transcript.first] = njh::replaceString(transcript.second.refAlnTranslation_.seq_, "-", "");
			}
			auto popCount = sampCountsForHaps.at(seqName.first);
			if(!njh::in(transcript.first, ret.proteinVariants_)){
				ret.proteinVariants_.emplace(transcript.first, VariantsInfo{transcript.first});
			}
			ret.proteinVariants_.at(transcript.first).addVariantInfo(transcript.second.refAlnTranslation_.seq_,
					transcript.second.queryAlnTranslation_.seq_, popCount, transcript.second.comp_);
		}
	}
	for(auto & varPerTrans : ret.proteinVariants_){
		varPerTrans.second.setFinals(rPars, totalPopCount);
	}
	return ret;
}



}  // namespace njhseq
