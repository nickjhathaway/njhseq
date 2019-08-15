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
	{
		//write out fasta file of input
		seqInfo seq;
		SeqInput reader(seqOpts);
		reader.openIn();
		SeqOutput writer(SeqIOOptions::genFastaOut(seqInputFnp));
		writer.openOut();
		while(reader.readNextRead(seq)){
			writer.write(seq);
		}
	}

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
	while (bReader.GetNextAlignment(bAln)) {
		if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
			auto balnGenomicRegion = GenomicRegion(bAln, refData);
			auto results = std::make_shared<AlignmentResults>(bAln, refData, true);
			results->setRefSeq(tReader);
			results->setComparison(true);
			ret.seqAlns_[bAln.Name].emplace_back(results);
			if (!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)) {
				continue;
			}
			for (const auto & g : alnRegionToGeneIds.at(balnGenomicRegion.createUidFromCoords())) {
				const auto & currentGene = genes.at(g);
				const auto & currentGeneInfo = ret.transcriptInfosForGene_.at(g);
				auto translations = translateBasedOnAlignment(bAln, *currentGene, currentGeneInfo, tReader, alignObj, refData);
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
			uint32_t queryAlnStart = aln->alnSeqAligned_->seq_.find_first_not_of("-");
			uint32_t queryAlnEnd = aln->alnSeqAligned_->seq_.find_last_not_of("-");
			auto popCount = sampCountsForHaps.at(seqName.first);
			for(const auto & seqPos : iter::range(queryAlnStart, queryAlnEnd + 1)){
				if('-' != aln->refSeqAligned_->seq_[seqPos]){
					uint32_t seqChromPosition = getRealPosForAlnPos(aln->refSeqAligned_->seq_, seqPos) + aln->gRegion_.start_;
					ret.baseForPosition_[aln->gRegion_.chrom_][seqChromPosition] = aln->refSeqAligned_->seq_[seqPos];
					ret.seqVariants_[aln->gRegion_.chrom_].allBases[seqChromPosition][aln->alnSeqAligned_->seq_[seqPos]] += popCount;
				}
			}
			for(const auto & m : aln->comp_.distances_.mismatches_){
				ret.seqVariants_[aln->gRegion_.chrom_].snps[m.second.refBasePos][m.second.seqBase]+= popCount;
			}
			for(const auto & gap : aln->comp_.distances_.alignmentGaps_){
				if(gap.second.ref_){
					//insertion
					ret.seqVariants_[aln->gRegion_.chrom_].insertions[gap.second.refPos_ - 1][gap.second.gapedSequence_]+=popCount;
				}else{
					//deletion
					ret.seqVariants_[aln->gRegion_.chrom_].deletions[gap.second.refPos_ - 1][gap.second.gapedSequence_]+=popCount;
				}
			}
		}
	}

	for(auto & varPerChrom : ret.seqVariants_){
		//filter saps and indels by occurrence cut off
		for(auto & snp : varPerChrom.second.snps){
			for(const auto & b : snp.second){
				if(b.second < rPars.occurrenceCutOff){
					continue;
				}
				varPerChrom.second.snpsFinal[snp.first][b.first] = b.second;
				if(b.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
					varPerChrom.second.variablePositons_.emplace(snp.first);
				}
			}
		}
		for(const auto & del : varPerChrom.second.deletions){
			for(const auto & d : del.second){
				if(d.second < rPars.occurrenceCutOff){
					continue;
				}
				varPerChrom.second.deletionsFinal[del.first][d.first] = d.second;
				if(d.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
					varPerChrom.second.variablePositons_.emplace(del.first);
				}
			}
		}
		for(const auto & ins : varPerChrom.second.insertions){
			for(const auto & i : ins.second){
				if(i.second < rPars.occurrenceCutOff){
					continue;
				}
				varPerChrom.second.insertionsFinal[ins.first][i.first] = i.second;
				if(i.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
					varPerChrom.second.variablePositons_.emplace(ins.first);
				}
			}
		}
	}


	//index amino acid changes per transcript
	for(const auto & seqName : ret.translations_){
		for(const auto & transcript : seqName.second){
			if("" == ret.proteinForTranscript_[transcript.first]){
				ret.proteinForTranscript_[transcript.first] = njh::replaceString(transcript.second.refAlnTranslation_.seq_, "-", "");
			}
			uint32_t queryAlnStart = transcript.second.queryAlnTranslation_.seq_.find_first_not_of("-");
			uint32_t queryAlnEnd = transcript.second.queryAlnTranslation_.seq_.find_last_not_of("-");
			auto popCount = sampCountsForHaps.at(seqName.first);
			for(const auto & refProteinPos : iter::range(queryAlnStart, queryAlnEnd + 1)){
				if('-' != transcript.second.refAlnTranslation_.seq_[refProteinPos]){
					ret.proteinVariants_[transcript.first].allBases[getRealPosForAlnPos(transcript.second.refAlnTranslation_.seq_, refProteinPos)][transcript.second.queryAlnTranslation_.seq_[refProteinPos]] += popCount;
				}
			}
			for(const auto & m : transcript.second.comp_.distances_.mismatches_){
				ret.proteinVariants_[transcript.first].snps[m.second.refBasePos][m.second.seqBase]+= popCount;
			}
			for(const auto & gap : transcript.second.comp_.distances_.alignmentGaps_){
				if(gap.second.ref_){
					//insertion
					ret.proteinVariants_[transcript.first].insertions[gap.second.refPos_ - 1][gap.second.gapedSequence_]+=popCount;
				}else{
					//deletion
					ret.proteinVariants_[transcript.first].deletions[gap.second.refPos_ - 1][gap.second.gapedSequence_]+=popCount;
				}
			}
		}
	}

	for(auto & varPerTrans : ret.proteinVariants_){
		//filter saps and indels by occurrence cut off
		for(auto & snp : varPerTrans.second.snps){
			for(const auto & b : snp.second){
				if(b.second < rPars.occurrenceCutOff){
					continue;
				}
				varPerTrans.second.snpsFinal[snp.first][b.first] = b.second;
				if(b.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
					varPerTrans.second.variablePositons_.emplace(snp.first);
				}
			}
		}
		for(const auto & del : varPerTrans.second.deletions){
			for(const auto & d : del.second){
				if(d.second < rPars.occurrenceCutOff){
					continue;
				}
				varPerTrans.second.deletionsFinal[del.first][d.first] = d.second;
				if(d.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
					varPerTrans.second.variablePositons_.emplace(del.first);
				}
			}
		}
		for(const auto & ins : varPerTrans.second.insertions){
			for(const auto & i : ins.second){
				if(i.second < rPars.occurrenceCutOff){
					continue;
				}
				varPerTrans.second.insertionsFinal[ins.first][i.first] = i.second;
				if(i.second/static_cast<double>(totalPopCount) > rPars.lowVariantCutOff){
					varPerTrans.second.variablePositons_.emplace(ins.first);
				}
			}
		}
	}


	return ret;
}



}  // namespace njhseq
