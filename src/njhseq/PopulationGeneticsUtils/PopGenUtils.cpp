/*
 * PopGenUtils.cpp
 *
 *  Created on: Oct 16, 2021
 *      Author: nick
 */


#include "PopGenUtils.hpp"

namespace njhseq {



TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars, const std::vector<seqInfo> & input){
  //read in meta if available
  std::unique_ptr<MultipleGroupMetaData> meta;
  if("" != pars.metaFnp){
    meta = std::make_unique<MultipleGroupMetaData>(pars.metaFnp);
  }

  std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(pars.ignoreSubFields);
  auto inputSeqs = CollapsedHaps::collapseReads(input, meta, metaValuesToAvoid);
  return collapseAndCallVariants(pars, inputSeqs);
}

TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars){
  //read in meta if available
  std::unique_ptr<MultipleGroupMetaData> meta;
  if("" != pars.metaFnp){
    meta = std::make_unique<MultipleGroupMetaData>(pars.metaFnp);
  }
  std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(pars.ignoreSubFields);
  auto inputSeqs = CollapsedHaps::readInReads(pars.inOpts, meta, metaValuesToAvoid);
  return collapseAndCallVariants(pars, inputSeqs);
}


TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars, CollapsedHaps & inputSeqs){

	njh::files::checkExistenceThrow(pars.transPars.lzPars_.genomeFnp, __PRETTY_FUNCTION__);
	njh::files::makeDirP(njh::files::MkdirPar(pars.outputDirectory.parent_path()));
	//kinda silly but this will make so any parent directories that need to exist will be made and then
	//overwrite directory will take affect below
	njh::files::makeDir(njh::files::MkdirPar(pars.outputDirectory, pars.overWriteDirectory));
	//samples names
	auto sampNamesPerSeq = inputSeqs.getSampleNamesPerSeqs();


	auto allSamples = inputSeqs.getAllSampleNames();
	//rename based on freq
	inputSeqs.renameBaseOnFreq(pars.identifier);
	//write out seqs
	auto uniqueSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(pars.outputDirectory, "uniqueSeqs.fasta.gz"));
	inputSeqs.writeOutAll(pars.outputDirectory, "uniqueSeqs");
	//key1 = sample, key2 = hap, value = readCount
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> samplesToHapsWithReadCnts;
	VecStr possibleCounts{"readCount", "barcodeCount"};
	for(const auto & e: iter::enumerate(inputSeqs.names_)) {
		for(const auto & name : e.second) {
			auto samp = CollapsedHaps::getSampleNameFromSeqName(name);
			uint32_t readCount = 1;
			if(MetaDataInName::nameHasMetaData(name)) {
				MetaDataInName nameMeta(name);
				for(const auto & posCntField : possibleCounts ) {
					if(nameMeta.containsMeta(posCntField)) {
						readCount = nameMeta.getMeta<uint32_t>(posCntField);
					}
				}
			}
			samplesToHapsWithReadCnts[samp][inputSeqs.seqs_[e.index]->name_] = readCount;
		}
		// std::cout << e.index << std::endl;
		// std::cout << "\t" << njh::conToStr(e.element, ",") << std::endl;
	}
	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);

	std::shared_ptr<aligner> alignerObj = std::make_shared<aligner>(maxLen, gapScoringParameters(7,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj->weighHomopolymers_ = false;
	alignerObj->processAlnInfoInput(pars.alnCacheDir.string(), false);

	std::unordered_map<std::string, uint32_t> seqNameKey = inputSeqs.genSeqNameKey();


	auto variantInfoDir =  njh::files::make_path(pars.outputDirectory, "variantCalls");
	njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
	std::unique_ptr<TranslatorByAlignment> translator = std::make_unique<TranslatorByAlignment>(pars.transPars);
	//translator->pars_.keepTemporaryFiles_ = true;
	translator->pars_.workingDirtory_ = variantInfoDir;
	std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;
	for(const auto pos : iter::range(inputSeqs.size())){
		sampNamesForPopHaps[inputSeqs.seqs_[pos]->name_] = sampNamesPerSeq[pos];
	}
	//pars.transPars.
	//pars.transPars.additionalBowtieArguments_

	translator->pars_.additionalBowtieArguments_ = njh::pasteAsStr(translator->pars_.additionalBowtieArguments_, " -p ", pars.numThreads);
	auto translatedRes = translator->run(SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName()), sampNamesForPopHaps, pars.variantCallerRunPars);
	// std::cout << njh::bashCT::green;
	// for(const auto & seqTransPerScript : translatedRes.translations_) {
	// 	for(const auto & seqTrans : seqTransPerScript.second) {
	// 		std::cout << seqTrans.first << std::endl;
	// 		for(const auto & mis : seqTrans.second.comp_.distances_.mismatches_) {
	// 			std::cout << "\t" << mis.first  << ":" << mis.second.seqBase << std::endl;
	// 		}
	// 		for(const auto & g : seqTrans.second.comp_.distances_.alignmentGaps_) {
	// 			std::cout << "\t" << g.first  << ":" << g.second.gapedSequence_ << std::endl;
	// 		}
	// 	}
	// }
	// std::cout << njh::bashCT::red;
	// for(const auto & seqAlns : translatedRes.seqAlns_) {
	// 	for(const auto & seqAln : seqAlns.second) {
	// 		std::cout << seqAln.refSeq_.name_ << ": " << seqAln.querySeq_.name_ << std::endl;
	// 		seqAln.alnRefSeq_.outPutSeqAnsi(std::cout);
	// 		seqAln.alnQuerySeq_.outPutSeqAnsi(std::cout);
	// 		for(const auto & mis : seqAln.comp_.distances_.mismatches_) {
	// 			std::cout << "\t" << seqAln.gRegion_.start_ + mis.first  << ":" << mis.second.seqBase << std::endl;
	// 		}
	// 		for(const auto & g : seqAln.comp_.distances_.alignmentGaps_) {
	// 			std::cout << "\t" << seqAln.gRegion_.start_ + g.first  << ":" << g.second.gapedSequence_ << std::endl;
	// 		}
	// 	}
	// }
	// std::cout << njh::bashCT::reset;


  OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "inputSeqs.bed"));
	translatedRes.writeOutSeqAlnIndvVars(njh::files::make_path(variantInfoDir, "variantsPerSeqAln.tab.txt.gz"));
	translatedRes.writeSeqLocations(popBedLocs);

	std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;

	OutputStream seqsUnableToBeMappedOut(njh::files::make_path(variantInfoDir, "seqsUnableToBeMapped.txt"));
	seqsUnableToBeMappedOut << njh::conToStr(translatedRes.seqsUnableToBeMapped_, "\n") << std::endl;
	OutputStream seqsTranslationFilteredOut(njh::files::make_path(variantInfoDir, "seqsTranslationFiltered.txt"));
	seqsTranslationFilteredOut << njh::conToStr(translatedRes.seqsTranslationFiltered_, "\n") << std::endl;

	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(!translatedRes.translations_.empty()){
		SeqOutput transwriter(SeqIOOptions::genFastaOutGz(njh::files::make_path(variantInfoDir, "translatedInput.fasta.gz")));
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::unordered_map<std::string, std::vector<seqInfo>> translatedSeqsByTranscript;
		std::unordered_map<std::string, std::vector<std::unordered_set<std::string>>> translatedSeqInputNames;
		auto seqNames = njh::getVecOfMapKeys(translatedRes.translations_);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		njh::sort(seqNames);
		for(const auto & seqName : seqNames){
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			for(const auto & transcript : translatedRes.translations_.at(seqName)){
				transwriter.openWrite(transcript.second.translation_);
				translatedSeqsByTranscript[transcript.first].emplace_back(transcript.second.translation_);
				translatedSeqsByTranscript[transcript.first].back().cnt_ = inputSeqs.names_[seqNameKey[seqName]].size();
				translatedSeqInputNames[transcript.first].emplace_back(inputSeqs.names_[seqNameKey[seqName]]);
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		OutputStream divMeasuresOut(njh::files::make_path(variantInfoDir, "translatedDivMeasures.tab.txt"));
		divMeasuresOut << njh::conToStr(pars.calcPopMeasuresPars.genHeader(), "\t") << std::endl;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		///
		auto fullTypedAAForTranslated = translatedRes.translated_genAATypedStr();
		auto knownsTypedAAForTranslated = translatedRes.translated_genAATypedStrOnlyKnowns();
		auto variableTypedAAForTranslated = translatedRes.translated_genAATypedStrOnlyPopVariant();
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//
		for(const auto & translatedSeqs : translatedSeqsByTranscript){
			auto inputTranslatedSeq = CollapsedHaps::collapseReads(translatedSeqs.second, translatedSeqInputNames[translatedSeqs.first]);
			std::string identifierTranslated = njh::pasteAsStr(pars.identifier, "-translated");
			if(translatedSeqsByTranscript.size() > 1){
				identifierTranslated = njh::pasteAsStr(pars.identifier, "-", translatedSeqs.first, "-translated");
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			for(const auto  seqPos : iter::range(inputTranslatedSeq.seqs_.size())){
//				const auto & seq = inputTranslatedSeq.seqs_[seqPos];
//				std::cout << seq->name_ << std::endl;
//				std::cout << "\tfullTypedAAForTranslated: " << fullTypedAAForTranslated[seq->name_] << std::endl;
//			}
			auto renameRes = inputTranslatedSeq.renameBaseOnFreq(identifierTranslated);
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;

			//write out seqs
			inputTranslatedSeq.writeOutAll(variantInfoDir, njh::pasteAsStr(translatedSeqs.first, "-", "uniqueTranslatedSeqs"));
			//get div measures
			auto calcPopMeasuresPars =  pars.calcPopMeasuresPars;
			calcPopMeasuresPars.numSegSites_ = njh::mapAt(translatedRes.proteinVariants_, translatedSeqs.first).getFinalNumberOfSegratingSites();
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			auto divMeasures = inputTranslatedSeq.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			divMeasuresOut << njh::conToStr(divMeasures.getOut(inputTranslatedSeq, identifierTranslated, calcPopMeasuresPars), "\t")  << std::endl;
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			OutputStream outAATyped(njh::files::make_path(variantInfoDir, njh::pasteAsStr(translatedSeqs.first, "-", "translatedSeqsAATyped.tab.txt.gz") ) );
			outAATyped << "name\tfullTyped\tknownTyped\tvariantTyped" << std::endl;
			for(const auto & seq : inputTranslatedSeq.seqs_){
				outAATyped << seq->name_
						<< "\t" << fullTypedAAForTranslated[renameRes.newNameToOldNameKey_[seq->name_]]
						<< "\t" << knownsTypedAAForTranslated[renameRes.newNameToOldNameKey_[seq->name_]]
						<< "\t" << variableTypedAAForTranslated[renameRes.newNameToOldNameKey_[seq->name_]] << std::endl;
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		OutputStream transBedLocs(njh::files::make_path(variantInfoDir, "translatedInput.bed"));
		translatedRes.writeSeqLocationsTranslation(transBedLocs);
		translatedRes.writeOutTranslatedIndvVars(njh::files::make_path(variantInfoDir, "variantsPerTranslatedSeq.tab.txt.gz"), translator->knownAminoAcidPositions_);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			//protein
			for(auto & varPerTrans : translatedRes.proteinVariants_){
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;

				{
					//writing vcfs
					auto vcfOutputForTrans = varPerTrans.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf")));
					//vcfOut << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read Depth for the ref and alt alleles in the order listed, a count of 0 means not detected\">" << std::endl;
					vcfOutputForTrans.formatEntries_.emplace_back(
						"AT", "1", "Integer",
						"Total Read Depth for this sample, a count of 0 means no coverage in this sample"
					);
					vcfOutputForTrans.formatEntries_.emplace_back(
						"AD", "R", "Integer",
						"Read Depth for the ref and alt alleles in the order listed, a count of 0 means not detected"
					);
					vcfOutputForTrans.formatEntries_.emplace_back(
						"AF", "R", "Float",
						"Read Frequncy for the ref and alt alleles in the order listed, a freq of 0 means not detected"
					);
					std::unordered_set<std::string> chromPositions;
					for(const auto & rec : vcfOutputForTrans.records) {
						chromPositions.emplace(njh::pasteAsStr(rec.chrom_, "-", rec.pos_));
					}
					//key1 = haplotypeName, key2 = chrom, key3 = vcf-based positioning, value = ref,alt
					std::unordered_map<std::string, std::map<std::string, std::map<uint32_t, std::pair<std::string, std::string>>>> vcfAlts;

					for(auto & translatedSeqRes : translatedRes.translations_) {
						if(!njh::in(varPerTrans.first, translatedSeqRes.second)) {
							continue;
						}
						// std::cout << translatedSeqRes.first << std::endl;
						//adding SNPs

						for(const auto & mis : translatedSeqRes.second[varPerTrans.first].comp_.distances_.mismatches_) {
							//adjust for genomic location and for the 1 based positioning of vcf
							auto realTranslatedPos = mis.second.refBasePos;
							auto vcfPosition = realTranslatedPos + 1;
							auto currentVariantChromPos = njh::pasteAsStr(varPerTrans.first, "-", vcfPosition);
							if(njh::in(currentVariantChromPos, chromPositions)) {
								auto ref = varPerTrans.second.getBaseForGenomicRegion(realTranslatedPos);
								auto alt = mis.second.seqBase;
								vcfAlts[translatedSeqRes.first][varPerTrans.first][vcfPosition] = std::make_pair(std::string(1, ref), std::string(1,alt));
							}
						}
						for(const auto & g : translatedSeqRes.second[varPerTrans.first].comp_.distances_.alignmentGaps_) {
							if(g.second.ref_) {
								//insertion
								//substract 1 because vcf does insertions/deletions from the base directly proceding the actual INDEL
								auto realGenomicPos = g.second.refPos_ - 1;
								auto vcfPosition = realGenomicPos + 1;
								auto currentVariantChromPos = njh::pasteAsStr(varPerTrans.first, "-", vcfPosition);
								if(njh::in(currentVariantChromPos, chromPositions)) {
									auto ref = varPerTrans.second.getBaseForGenomicRegion(realGenomicPos);
									auto alt = njh::pasteAsStr(varPerTrans.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
									vcfAlts[translatedSeqRes.first][varPerTrans.first][vcfPosition] = std::make_pair(std::string(1, ref), alt);
								}
							} else {
								//deletion
								//substract 1 because vcf does insertions/deletions from the base directly proceding the actual INDEL
								auto realGenomicPos = g.second.refPos_ - 1;
								auto vcfPosition = realGenomicPos + 1;
								auto currentVariantChromPos = njh::pasteAsStr(varPerTrans.first, "-", vcfPosition);
								if(njh::in(currentVariantChromPos, chromPositions)) {
									auto ref = njh::pasteAsStr(varPerTrans.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
									auto alt = varPerTrans.second.getBaseForGenomicRegion(realGenomicPos);
									vcfAlts[translatedSeqRes.first][varPerTrans.first][vcfPosition] = std::make_pair(ref, std::string(1, alt));
								}
							}
						}
						// for(const auto & chrom : vcfAlts[translatedSeqRes.first]) {
						// 	for(const auto & pos : chrom.second) {
						// 		std::cout << chrom.first << "\t" << pos.first << "\t" << pos.second.first << "\t" << pos.second.second << std::endl;
						// 	}
						// }
						// std::cout << std::endl;
					}


					for (auto&rec: vcfOutputForTrans.records) {
						for (const auto&sample: allSamples) {
							std::vector<uint32_t> dps(1 + rec.alts_.size(), 0);
							// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
							for(const auto & haps : njh::mapAt(samplesToHapsWithReadCnts, sample)) {
								// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
								bool foundAlt = false;
								for(const auto & alt : iter::enumerate(rec.alts_)) {
									if(njh::in(rec.pos_, vcfAlts[haps.first][rec.chrom_]) &&
										rec.ref_ == vcfAlts[haps.first][rec.chrom_][rec.pos_].first &&
										alt.second == vcfAlts[haps.first][rec.chrom_][rec.pos_].second){
										dps[1 + alt.index] += haps.second;
										foundAlt = true;
									}
								}
								bool coveredByHap = false;
								// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
								if(njh::in(haps.first, translatedRes.translations_)  &&
								  njh::in(varPerTrans.first, translatedRes.translations_.at(haps.first)) &&
									translatedRes.translations_[haps.first].at(varPerTrans.first).genBedRec().chrom_ == rec.chrom_ &&
									rec.pos_ + 1 >= translatedRes.translations_[haps.first].at(varPerTrans.first).genBedRec().chromStart_ &&
									rec.pos_ + 1 < translatedRes.translations_[haps.first].at(varPerTrans.first).genBedRec().chromEnd_) {
									coveredByHap = true;
								}
								if(!foundAlt && coveredByHap) {
									//no alternative found, and referene for this position is covered by hap, increase depth for reference
									dps[0] += haps.second;
								}
							}
							auto dpsSum = vectorSum(dps);
							std::vector<double> dpsFreq;
							dpsFreq.reserve(dps.size());
							for (const auto dp: dps) {
								dpsFreq.emplace_back(dp / dpsSum);
							}
							rec.sampleFormatInfos_[sample].addMeta("AT", dpsSum);
							rec.sampleFormatInfos_[sample].addMeta("AD", njh::conToStr(dps, ","));
							rec.sampleFormatInfos_[sample].addMeta("AF", njh::conToStr(dpsFreq, ","));

						}
					}
					{
						OutputStream genomeVcfWithSamples(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-proteinWithSampleInfo.vcf")));
						vcfOutputForTrans.writeOutFixedAndSampleMeta(genomeVcfWithSamples);
					}
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				varPerTrans.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt.gz")), varPerTrans.first, true);
				std::set<uint32_t> knownMutationsLocationsZeroBased;
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				for(const auto & snpPos : varPerTrans.second.allBases){
					if(njh::in(snpPos.first + 1, translator->knownAminoAcidPositions_[varPerTrans.first])){
						knownMutationsLocationsZeroBased.emplace(snpPos.first);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						auto genomicLocationForAAPos = translatedRes.translationInfoForTranscirpt_.at(varPerTrans.first)->genBedFromAAPositions(snpPos.first, snpPos.first + 1);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						for(const auto gPos : iter::range(genomicLocationForAAPos.chromStart_, genomicLocationForAAPos.chromEnd_)){
							knownAAMutsChromPositions[genomicLocationForAAPos.chrom_].emplace(gPos);
						}
					}
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				varPerTrans.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt.gz")), varPerTrans.first, true);
				if(!varPerTrans.second.variablePositons_.empty()){
					GenomicRegion variableRegion = varPerTrans.second.getVariableRegion();
					variableRegion.start_ += 1; //do one based positioning
					OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
					bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				if(!knownMutationsLocationsZeroBased.empty()){
					varPerTrans.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt.gz")), varPerTrans.first, knownMutationsLocationsZeroBased, true);
				}
			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			translatedRes.writeOutAATypedInfo(njh::files::make_path(variantInfoDir, "seqsAATyped.tab.txt.gz"));
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}


	//snps
	uint32_t maxSeqCount = 0;
	auto calcPopMeasuresPars =  pars.calcPopMeasuresPars;
	for(auto & varPerChrom : translatedRes.seqVariants_){
		for(const auto & count : varPerChrom.second.depthPerPosition){
			if(count.second > maxSeqCount){
				//cheap way of doing this for now
				maxSeqCount = count.second;
				calcPopMeasuresPars.numSegSites_ = varPerChrom.second.getFinalNumberOfSegratingSites();
			}
		}

		{
			//writing vcfs
			auto vcfOutputForChrom = varPerChrom.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf")) );
			//vcfOut << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read Depth for the ref and alt alleles in the order listed, a count of 0 means not detected\">" << std::endl;

			vcfOutputForChrom.formatEntries_.emplace_back(
				"AT", "1", "Integer",
				"Total Read Depth for this sample, a count of 0 means no coverage in this sample"
			);
			vcfOutputForChrom.formatEntries_.emplace_back(
				"AD", "R", "Integer",
				"Read Depth for the ref and alt alleles in the order listed, a count of 0 means not detected"
			);
			vcfOutputForChrom.formatEntries_.emplace_back(
				"AF", "R", "Float",
				"Read Frequncy for the ref and alt alleles in the order listed, a freq of 0 means not detected"
			);
			std::unordered_set<std::string> chromPositions;
			for(const auto & rec : vcfOutputForChrom.records) {
				chromPositions.emplace(njh::pasteAsStr(rec.chrom_, "-", rec.pos_));
			}
			//key1 = haplotypeName, key2 = chrom, key3 = vcf-based positioning, value = ref,alt
			std::unordered_map<std::string, std::map<std::string, std::map<uint32_t, std::pair<std::string, std::string>>>> vcfAlts;

			for(const auto & seqAlns : translatedRes.seqAlns_) {
				for(const auto & seqAln : seqAlns.second) {
					// std::cout << seqAln.querySeq_.name_ << std::endl;
					// seqAln.alnRefSeq_.outPutSeqAnsi(std::cout);
					// seqAln.alnQuerySeq_.outPutSeqAnsi(std::cout);
					//adding SNPs
					for(const auto & mis : seqAln.comp_.distances_.mismatches_) {
						//adjust for genomic location and for the 1 based positioning of vcf
						auto realGenomicPos = seqAln.gRegion_.start_ + mis.second.refBasePos;
						auto vcfPosition = realGenomicPos + 1;
						auto currentVariantChromPos = njh::pasteAsStr(seqAln.gRegion_.chrom_, "-", vcfPosition);
						if(njh::in(currentVariantChromPos, chromPositions)) {
							auto ref = varPerChrom.second.getBaseForGenomicRegion(realGenomicPos);
							auto alt = mis.second.seqBase;
							vcfAlts[seqAlns.first][seqAln.gRegion_.chrom_][vcfPosition] = std::make_pair(std::string(1, ref), std::string(1,alt));
						}
					}
					for(const auto & g : seqAln.comp_.distances_.alignmentGaps_) {
						if(g.second.ref_) {
							//insertion
							//substract 1 because vcf does insertions/deletions from the base directly proceding the actual INDEL
							auto realGenomicPos = seqAln.gRegion_.start_ + g.second.refPos_ - 1;
							auto vcfPosition = realGenomicPos + 1;
							auto currentVariantChromPos = njh::pasteAsStr(seqAln.gRegion_.chrom_, "-", vcfPosition);
							if(njh::in(currentVariantChromPos, chromPositions)) {
								auto ref = varPerChrom.second.getBaseForGenomicRegion(realGenomicPos);
								auto alt = njh::pasteAsStr(varPerChrom.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
								vcfAlts[seqAlns.first][seqAln.gRegion_.chrom_][vcfPosition] = std::make_pair(std::string(1, ref), alt);
							}
						} else {
							//deletion
							//substract 1 because vcf does insertions/deletions from the base directly proceding the actual INDEL
							auto realGenomicPos = seqAln.gRegion_.start_ + g.second.refPos_ - 1;
							auto vcfPosition = realGenomicPos + 1;
							auto currentVariantChromPos = njh::pasteAsStr(seqAln.gRegion_.chrom_, "-", vcfPosition);
							if(njh::in(currentVariantChromPos, chromPositions)) {
								auto ref = njh::pasteAsStr(varPerChrom.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
								auto alt = varPerChrom.second.getBaseForGenomicRegion(realGenomicPos);
								vcfAlts[seqAlns.first][seqAln.gRegion_.chrom_][vcfPosition] = std::make_pair(ref, std::string(1, alt));
							}
						}
					}
					// for(const auto & chrom : vcfAlts[seqAlns.first]) {
					// 	for(const auto & pos : chrom.second) {
					// 		std::cout << chrom.first << "\t" << pos.first << "\t" << pos.second.first << "\t" << pos.second.second << std::endl;
					// 	}
					// }
					// std::cout << std::endl;
				}
			}

			for (auto&rec: vcfOutputForChrom.records) {
				for (const auto&sample: allSamples) {
					std::vector<uint32_t> dps(1 + rec.alts_.size(), 0);
					// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
					for(const auto & haps : njh::mapAt(samplesToHapsWithReadCnts, sample)) {
						// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
						bool foundAlt = false;
						for(const auto & alt : iter::enumerate(rec.alts_)) {
							if(njh::in(rec.pos_, vcfAlts[haps.first][rec.chrom_]) &&
								rec.ref_ == vcfAlts[haps.first][rec.chrom_][rec.pos_].first &&
								alt.second == vcfAlts[haps.first][rec.chrom_][rec.pos_].second){
								dps[1 + alt.index] += haps.second;
								foundAlt = true;
							}
						}
						bool coveredByHap = false;
						// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
						if(njh::in(haps.first, translatedRes.seqAlns_)) {
							for(const auto & seqAln : njh::mapAt(translatedRes.seqAlns_, haps.first)) {
								// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
								if(seqAln.gRegion_.chrom_ == rec.chrom_ && rec.pos_ +1 >= seqAln.gRegion_.start_ && rec.pos_ + 1 < seqAln.gRegion_.end_) {
									coveredByHap = true;
									break;
								}
							}
						}
						if(!foundAlt && coveredByHap) {
							//no alternative found, and referene for this position is covered by hap, increase depth for reference
							dps[0] += haps.second;
						}
					}

					auto dpsSum = vectorSum(dps);
					std::vector<double> dpsFreq;
					dpsFreq.reserve(dps.size());
					for (const auto dp: dps) {
						dpsFreq.emplace_back(dp / dpsSum);
					}
					rec.sampleFormatInfos_[sample].addMeta("AT", dpsSum);
					rec.sampleFormatInfos_[sample].addMeta("AD", njh::conToStr(dps, ","));
					rec.sampleFormatInfos_[sample].addMeta("AF", njh::conToStr(dpsFreq, ","));

				}
			}
			{
				OutputStream genomeVcfWithSamples(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomicWithSampleInfo.vcf")));
				vcfOutputForChrom.writeOutFixedAndSampleMeta(genomeVcfWithSamples);
			}
		}
		varPerChrom.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt.gz")), varPerChrom.first);
		if(!knownAAMutsChromPositions[varPerChrom.first].empty()){
			varPerChrom.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt.gz")), varPerChrom.first, knownAAMutsChromPositions[varPerChrom.first]);
		}
		varPerChrom.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt.gz")), varPerChrom.first);
		if(!varPerChrom.second.variablePositons_.empty()){
			GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
			OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
			bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	{
		auto snpTyped = translatedRes.genSeqSNPTypedStr();
		OutputStream snpTypedOut(njh::files::make_path(variantInfoDir, "seqSNPTyped.tab.txt.gz"));
		snpTypedOut << "name\tsnpTyped" << std::endl;
		for(const auto & seqPos : inputSeqs.getOrderByTopCnt()){
			snpTypedOut << inputSeqs.seqs_[seqPos]->name_
					<< "\t" << snpTyped[inputSeqs.seqs_[seqPos]->name_];
			snpTypedOut << std::endl;
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	{
		auto divMeasures = inputSeqs.getGeneralMeasuresOfDiversity(
				calcPopMeasuresPars, alignerObj);

		divMeasures.writeDivMeasures(
				njh::files::make_path(pars.outputDirectory, "divMeasures.tab.txt"),
				inputSeqs, pars.identifier, calcPopMeasuresPars);
	}

	if(!pars.metaFieldsToCalcPopDiffs.empty()){
		auto outputDirPerMeta = njh::files::make_path(pars.outputDirectory, "perMetaFields");
		njh::files::makeDir(njh::files::MkdirPar(outputDirPerMeta));
		std::unordered_map<std::string, std::unordered_map<std::string, CollapsedHaps::GenPopMeasuresRes>> measuresPer;

		for(const auto & metaField : pars.metaFieldsToCalcPopDiffs){
			OutputStream divMeasuresOut(njh::files::make_path(outputDirPerMeta,
								metaField + "_divMeasures.tab.txt"));
			divMeasuresOut << njh::conToStr(calcPopMeasuresPars.genHeader(), "\t") << std::endl;
			auto splitSeqs = inputSeqs.splitOutSeqsByMeta(metaField);
			for(const auto & subField : splitSeqs){
//				std::cout << subField.first << std::endl;
//				for(const auto & seqIter : iter::enumerate(subField.second.seqs_)){
//					auto & seq = seqIter.element;
//					std::cout << "\t" << seq->name_ << std::endl;
//					std::cout << "\t" << seq->cnt_ << std::endl;
//					std::cout << "\t" << seq->frac_ << std::endl;
//					std::cout << "\t" << subField.second.names_[seqIter.index].size() << std::endl;
//				}
				auto divMeasures = subField.second.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
				divMeasuresOut << njh::conToStr(divMeasures.getOut(subField.second, njh::pasteAsStr(pars.identifier, "::", subField.first), calcPopMeasuresPars), "\t") << std::endl;
			}
		}
	}

	alignerObj->processAlnInfoOutput(pars.alnCacheDir.string(), false);

	return translatedRes;
}




}  // namespace njhseq
