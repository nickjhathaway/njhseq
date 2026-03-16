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
  if(!pars.metaFnp.empty()){
    meta = std::make_unique<MultipleGroupMetaData>(pars.metaFnp);
  }

  std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(pars.ignoreSubFields);
  auto inputSeqs = CollapsedHaps::collapseReads(input, meta, metaValuesToAvoid);
  return collapseAndCallVariants(pars, inputSeqs);
}

TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars){
  //read in meta if available
  std::unique_ptr<MultipleGroupMetaData> meta;
  if(!pars.metaFnp.empty()){
    meta = std::make_unique<MultipleGroupMetaData>(pars.metaFnp);
  }
  std::unordered_map<std::string, std::set<std::string>> metaValuesToAvoid = njh::progutils::CmdArgs::sepSubArgsMulti<std::string, std::string>(pars.ignoreSubFields);
  auto inputSeqs = CollapsedHaps::readInReads(pars.inOpts, meta, metaValuesToAvoid);
  return collapseAndCallVariants(pars, inputSeqs);
}


TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars, CollapsedHaps & inputSeqs){

	njh::stopWatch watch;
	watch.setLapName("start");
	njh::files::checkExistenceThrow(pars.transPars.lzPars_.genomeFnp, __PRETTY_FUNCTION__);
	njh::files::makeDirP(njh::files::MkdirPar(pars.outputDirectory.parent_path()));
	//kinda silly but this will make so any parent directories that need to exist will be made and then
	//overwrite directory will take affect below
	njh::files::makeDir(njh::files::MkdirPar(pars.outputDirectory, pars.overWriteDirectory));
	//samples names
	auto sampNamesPerSeq = inputSeqs.getSampleNamesPerSeqs();

	// auto allSamples = inputSeqs.getAllSampleNames();
	//rename based on freq
	inputSeqs.renameBaseOnFreq(pars.identifier);
	//write out seqs
	auto uniqueSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(pars.outputDirectory, "uniqueSeqs.fasta.gz"));
	inputSeqs.writeOutAll(pars.outputDirectory, "uniqueSeqs");
	if(pars.exportLabIsolateSeqs) {
		auto labIsolateSeqsOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(pars.outputDirectory, "refSeqs.fasta.gz"));
		inputSeqs.writeOutLabIsolateSeqs(labIsolateSeqsOpts);
	}
	//key1 = sample, key2 = hap, value = readCount
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> samplesToHapsWithReadCnts;
	//key1 = hap , key2 = sample, value = readCount
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapsToSamplesWithReadCnts;
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
			hapsToSamplesWithReadCnts[inputSeqs.seqs_[e.index]->name_][samp] = readCount;
		}
		// std::cout << e.index << std::endl;
		// std::cout << "\t" << njh::conToStr(e.element, ",") << std::endl;
	}
	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);
	watch.startNewLap("set up for variant calling");
	std::shared_ptr<aligner> alignerObj = std::make_shared<aligner>(maxLen, gapScoringParameters(7,1,0,0,0,0), substituteMatrix(2,-2), false);
	alignerObj->weighHomopolymers_ = false;
	alignerObj->processAlnInfoInput(pars.alnCacheDir.string(), false);

	std::unordered_map<std::string, uint32_t> seqNameKey = inputSeqs.genSeqNameKey();

	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	auto variantInfoDir =  njh::files::make_path(pars.outputDirectory, "variantCalls");
	njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
	std::unique_ptr<TranslatorByAlignment> translator = std::make_unique<TranslatorByAlignment>(pars.transPars);
	//translator->pars_.keepTemporaryFiles_ = true;
	translator->pars_.workingDirtory_ = variantInfoDir;
	// std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;
	// for(const auto pos : iter::range(inputSeqs.size())){
	// 	sampNamesForPopHaps[inputSeqs.seqs_[pos]->name_] = sampNamesPerSeq[pos];
	// }
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//pars.transPars.
	//pars.transPars.additionalBowtieArguments_
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	translator->pars_.additionalBowtieArguments_ = njh::pasteAsStr(translator->pars_.additionalBowtieArguments_, " -p ", pars.numThreads);
	TranslatorByAlignment::TranslatorByAlignmentResult translatedRes;
	watch.startNewLap("run variant calling and translation");
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//if a specific region is supplied, force alignment to that region
	if(!pars.refSeqRegion.chrom_.empty()) {
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		translatedRes = translator->run(SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName()), hapsToSamplesWithReadCnts, pars.refSeqRegion, pars.variantCallerRunPars);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	} else {
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		translatedRes = translator->run(SeqIOOptions::genFastaIn(uniqueSeqsOpts.out_.outName()), hapsToSamplesWithReadCnts, pars.variantCallerRunPars);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	}

	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	watch.startNewLap("set up for writing output");
	auto gprefix = bfs::path(translator->pars_.lzPars_.genomeFnp).replace_extension("");
	auto twoBitFnp = gprefix.string() + ".2bit";
	TwoBit::TwoBitFile tReader(twoBitFnp);
	auto contigLengths = tReader.getSeqLens();

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
	std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;

	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	{
		OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "inputSeqs.bed"));
		translatedRes.writeOutSeqAlnIndvVars(njh::files::make_path(variantInfoDir, "variantsPerSeqAln.tab.txt.gz"));
		translatedRes.writeSeqLocations(popBedLocs);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

		OutputStream seqsUnableToBeMappedOut(njh::files::make_path(variantInfoDir, "seqsUnableToBeMapped.txt"));
		seqsUnableToBeMappedOut << njh::conToStr(translatedRes.seqsUnableToBeMapped_, "\n") << std::endl;
		OutputStream seqsTranslationFilteredOut(njh::files::make_path(variantInfoDir, "seqsTranslationFiltered.txt"));
		seqsTranslationFilteredOut << njh::conToStr(translatedRes.seqsTranslationFiltered_, "\n") << std::endl;
	}

	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(!translatedRes.translations_.empty()){
		watch.startNewLap("writing translation output - initial");
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
		OutputStream translated_prev_freq_out(njh::files::make_path(variantInfoDir, "translated_prev_freq.tsv.gz"));
		translated_prev_freq_out << "target_name\tseq\tfreq\tprev" << std::endl;
		for(const auto & translatedSeqs : translatedSeqsByTranscript){
			watch.startNewLap(njh::pasteAsStr("writing translation output - ", translatedSeqs.first, " - collapse seq"));
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
			watch.startNewLap(njh::pasteAsStr("writing translation output - ", translatedSeqs.first, " - rename"));
			auto renameRes = inputTranslatedSeq.renameBaseOnFreq(identifierTranslated);
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;

			//write out seqs
			watch.startNewLap(njh::pasteAsStr("writing translation output - ", translatedSeqs.first, " - write out seqs"));
			inputTranslatedSeq.writeOutAll(variantInfoDir, njh::pasteAsStr(translatedSeqs.first, "-", "uniqueTranslatedSeqs"));
			//get div measures
			watch.startNewLap(njh::pasteAsStr("writing translation output - ", translatedSeqs.first, " - get div measures"));

			auto calcPopMeasuresPars =  pars.calcPopMeasuresPars;

			calcPopMeasuresPars.numSegSites_ = njh::mapAt(translatedRes.proteinVariants_, translatedSeqs.first).getFinalNumberOfSegregatingSites();
			// std::cout << __FILE__ << " " << __LINE__ << std::endl;
			auto divMeasures = inputTranslatedSeq.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
			// std::cout << __FILE__ << " " << __LINE__ << std::endl;
			divMeasuresOut << njh::conToStr(divMeasures.getOut(inputTranslatedSeq, identifierTranslated, calcPopMeasuresPars), "\t")  << std::endl;
			// std::cout << __FILE__ << " " << __LINE__ << std::endl;
			// std::cout << "variableTypedAAForTranslated names: " << njh::conToStr(njh::getVecOfMapKeys(variableTypedAAForTranslated), ",") << std::endl;

			auto inputTranslatedSeq_prevs = inputTranslatedSeq.getPrevalences();
			auto inputTranslatedSeq_freqs = inputTranslatedSeq.getWeightedAlleleFreqs();


			for (const auto & seq : inputTranslatedSeq.seqs_) {
				translated_prev_freq_out << identifierTranslated << "\t" << seq->seq_
					<< "\t" << inputTranslatedSeq_freqs[seq->seq_]
					<< "\t" << inputTranslatedSeq_prevs[seq->seq_] << std::endl;
			}

			watch.startNewLap(njh::pasteAsStr("writing translation output - ", translatedSeqs.first, " - translatedSeqsAATyped"));
			OutputStream outAATyped(njh::files::make_path(variantInfoDir, njh::pasteAsStr(translatedSeqs.first, "-", "translatedSeqsAATyped.tab.txt.gz") ) );
			outAATyped << "name\tfullTyped\tknownTyped\tvariantTyped" << std::endl;
			for(const auto & seq : inputTranslatedSeq.seqs_){
				// std::cout << '\t' << seq->name_ << std::endl;
				// std::cout << '\t' << renameRes.newNameToOldNameKey_[seq->name_] << std::endl;
				outAATyped << seq->name_
						<< "\t" << fullTypedAAForTranslated[renameRes.newNameToOldNameKey_[seq->name_]]
						<< "\t" << knownsTypedAAForTranslated[renameRes.newNameToOldNameKey_[seq->name_]]
						<< "\t" << variableTypedAAForTranslated[renameRes.newNameToOldNameKey_[seq->name_]] << std::endl;
			}
		}

		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr("writing translation output - ", "writeSeqLocationsTranslation"));

		{
			OutputStream transBedLocs(njh::files::make_path(variantInfoDir, "translatedInput.bed"));
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			translatedRes.writeSeqLocationsTranslation(transBedLocs);
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(njh::pasteAsStr("writing translation output - ", "writeOutTranslatedIndvVars"));
		translatedRes.writeOutTranslatedIndvVars(njh::files::make_path(variantInfoDir, "variantsPerTranslatedSeq.tab.txt.gz"), translator->knownAminoAcidPositions_);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			watch.startNewLap(njh::pasteAsStr("writing translation output - prior to - write out vcf fixed info"));
			//protein
			for(auto & varPerTrans : translatedRes.proteinVariants_){

				{
					watch.startNewLap(njh::pasteAsStr("writing translation output - ", varPerTrans.first, " - write out vcf fixed info"));
					//writing vcfs
					auto vcfOutputForTrans = varPerTrans.second.createVCFOutputFixed();
					vcfOutputForTrans.otherHeaderValuePairs_.emplace("reference", pars.transPars.lzPars_.genomeFnp.string());

					vcfOutputForTrans.addDefaultInfoField("TARGET", pars.identifier,VCFOutput::InfoEntry("TARGET", "1", "String", "the target that covers this variant"));
					vcfOutputForTrans.addDefaultInfoField("GeneID", njh::mapAt(translatedRes.translationInfoForTranscirpt_, varPerTrans.first)->geneID_, VCFOutput::InfoEntry("GeneID", "1", "String", "The Standardized Gene ID"));
					vcfOutputForTrans.addDefaultInfoField("GeneName", njh::mapAt(translatedRes.translationInfoForTranscirpt_, varPerTrans.first)->geneName_, VCFOutput::InfoEntry("GeneName", "1", "String", "A name for the Gene"));
					vcfOutputForTrans.addDefaultInfoField("TranscriptID", njh::mapAt(translatedRes.translationInfoForTranscirpt_, varPerTrans.first)->transcriptID_, VCFOutput::InfoEntry("TranscriptID", "1", "String", "The Transcript ID"));

					//add if location is a known drug resistance location
					vcfOutputForTrans.infoEntries_.emplace("AF", VCFOutput::InfoEntry(
						"KNOWN_AA_CHANGE_POS", "0", "Flag",
						"Position is a known amino acid position of interest, e.g. drug resistance mutation, etc"
					));
					for (auto & rec : vcfOutputForTrans.records_) {
						//if the position, which was entered as 1 based, is in 1 based position list of known muts
						if(njh::in(rec.pos_, translator->knownAminoAcidPositions_[varPerTrans.first])) {
							rec.info_.addMeta("KNOWN_AA_CHANGE_POS",true);
						}
					}

					watch.startNewLap(njh::pasteAsStr("writing translation output - ", varPerTrans.first, " - write out vcf sample info gather"));

					vcfOutputForTrans.contigEntries_.emplace(varPerTrans.first, VCFOutput::ContigEntry(varPerTrans.first, translatedRes.translationInfoForTranscirpt_[varPerTrans.first]->protein_.seq_.length()));
					{
						//auto vcfOutputForTrans = varPerTrans.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf")));
						// OutputStream vcfOut(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf")));
						// vcfOutputForTrans.writeOutFixedOnly(vcfOut);
					}
					vcfOutputForTrans.headerNonSampleFields_.emplace_back("FORMAT");
					vcfOutputForTrans.formatEntries_.emplace("GT", VCFOutput::FormatEntry(
						"GT", "1", "String",
						"Genotype"
					));
					vcfOutputForTrans.formatEntries_.emplace("DP", VCFOutput::FormatEntry(
						"DP", "1", "Integer",
						"Total Read Depth for this sample, a count of 0 means no coverage in this sample"
					));
					vcfOutputForTrans.formatEntries_.emplace("AD", VCFOutput::FormatEntry(
						"AD", "R", "Integer",
						"Read Depth for the ref and alt alleles in the order listed, a count of 0 means not detected"
					));
					vcfOutputForTrans.formatEntries_.emplace("AF", VCFOutput::FormatEntry(
						                                         "AF", "R", "Float",
						                                         "Read Frequency for the ref and alt alleles in the order listed, a freq of 0 means not detected"
					                                         ));


					std::unordered_set<std::string> chromPositions;
					for(const auto & rec : vcfOutputForTrans.records_) {
						chromPositions.emplace(njh::pasteAsStr(rec.chrom_, "-", rec.pos_));
					}
					//key1 = haplotypeName, key2 = chrom, key3 = vcf-based positioning, value = ref,alt
					std::unordered_map<std::string, std::map<std::string, std::map<uint32_t, std::vector<std::pair<std::string, std::string>>>>> vcfAlts;

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
								vcfAlts[translatedSeqRes.first][varPerTrans.first][vcfPosition].emplace_back(std::make_pair(std::string(1, ref), std::string(1,alt)));
							}
						}
						for(const auto & g : translatedSeqRes.second[varPerTrans.first].comp_.distances_.alignmentGaps_) {
							if(g.second.ref_) {
								//insertion
								//subtract 1 because vcf does insertions/deletions from the base directly proceeding the actual INDEL
								auto realGenomicPos = g.second.refPos_ - 1;
								auto vcfPosition = realGenomicPos + 1;
								auto currentVariantChromPos = njh::pasteAsStr(varPerTrans.first, "-", vcfPosition);
								if(njh::in(currentVariantChromPos, chromPositions)) {
									auto ref = varPerTrans.second.getBaseForGenomicRegion(realGenomicPos);
									auto alt = njh::pasteAsStr(varPerTrans.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
									vcfAlts[translatedSeqRes.first][varPerTrans.first][vcfPosition].emplace_back(std::make_pair(std::string(1, ref), alt));
								}
							} else {
								//deletion
								//subtract 1 because vcf does insertions/deletions from the base directly proceeding the actual INDEL
								auto realGenomicPos = g.second.refPos_ - 1;
								auto vcfPosition = realGenomicPos + 1;
								auto currentVariantChromPos = njh::pasteAsStr(varPerTrans.first, "-", vcfPosition);
								if(njh::in(currentVariantChromPos, chromPositions)) {
									auto ref = njh::pasteAsStr(varPerTrans.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
									auto alt = varPerTrans.second.getBaseForGenomicRegion(realGenomicPos);
									vcfAlts[translatedSeqRes.first][varPerTrans.first][vcfPosition].emplace_back(std::make_pair(ref, std::string(1, alt)));
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
					watch.startNewLap(njh::pasteAsStr("writing translation output - ", varPerTrans.first, " - write out vcf sample info add"));
					//add in sample names
					auto allSamplesForVariants = varPerTrans.second.getAllSamples();
					vcfOutputForTrans.samples_ = VecStr(allSamplesForVariants.begin(), allSamplesForVariants.end());
					for (auto&rec: vcfOutputForTrans.records_) {
						for (const auto&sample: allSamplesForVariants) {
							std::vector<uint32_t> dps(1 + rec.alts_.size(), 0);
							// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
							for(const auto & haps : njh::mapAt(samplesToHapsWithReadCnts, sample)) {
								// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
								bool foundAlt = false;
								for (const auto& alt: iter::enumerate(rec.alts_)) {
									if (njh::in(rec.pos_, vcfAlts[haps.first][rec.chrom_])) {
										for (const auto& alts: vcfAlts[haps.first][rec.chrom_][rec.pos_]) {
											if (rec.ref_ == alts.first &&
											    alt.second == alts.second) {
												dps[1 + alt.index] += haps.second;
												foundAlt = true;
											}
										}
									}
								}
								bool coveredByHap = false;
								// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
								if(njh::in(haps.first, translatedRes.translations_)  &&
								  njh::in(varPerTrans.first, translatedRes.translations_.at(haps.first)) &&
									translatedRes.translations_[haps.first].at(varPerTrans.first).genBedRec().chrom_ == rec.chrom_ &&
									rec.pos_ - 1 >= translatedRes.translations_[haps.first].at(varPerTrans.first).genBedRec().chromStart_ &&
									rec.pos_ - 1 < translatedRes.translations_[haps.first].at(varPerTrans.first).genBedRec().chromEnd_) {
									coveredByHap = true;
								}
								if(!foundAlt && coveredByHap) {
									//no alternative found, and reference for this position is covered by hap, increase depth for reference
									dps[0] += haps.second;
								}
							}
							auto dpsSum = vectorSum(dps);
							std::vector<double> dpsFreq;
							dpsFreq.reserve(dps.size());
							for (const auto dp: dps) {
								dpsFreq.emplace_back(dpsSum > 0 ? dp / dpsSum : 0.0);
							}
							if(0 == dpsSum) {
								//if dpsSum equals 0 that most likely means that the sample either had a haplotype that didn't cover
								//this location and/or a haplotype that didn't map
								//if the sample had a mapping haplotype that didn't cover then it either ended before this location or
								//it started after it, either way it's ok to mark as no info for this loc
								rec.sampleFormatInfos_[sample].addMeta("DP", ".");
								rec.sampleFormatInfos_[sample].addMeta("AD", njh::conToStr(std::vector<std::string>(dps.size(), "."), ","));
								rec.sampleFormatInfos_[sample].addMeta("AF", njh::conToStr(std::vector<std::string>(dpsFreq.size(), "."), ","));
							} else {
								rec.sampleFormatInfos_[sample].addMeta("DP", dpsSum);
								rec.sampleFormatInfos_[sample].addMeta("AD", njh::conToStr(dps, ","));
								rec.sampleFormatInfos_[sample].addMeta("AF", njh::conToStr(dpsFreq, ","));
							}
						}
					}

					{
						watch.startNewLap(njh::pasteAsStr("writing translation output - ", varPerTrans.first, " - write out vcf actual writing"));

						// OutputStream genomeVcfWithSamples(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-proteinWithSampleInfo.vcf")));
						OutputStream genomeVcfWithSamples(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf.gz")));
						vcfOutputForTrans.allAddGTFields(pars.variantCallerRunPars.ploidy);
						vcfOutputForTrans.allAutoAddDPFields();
						vcfOutputForTrans.allAutoAddTYPEFields();
						vcfOutputForTrans.allAutoAdd_AN_AC_AF_InfoFields();
						vcfOutputForTrans.allAutoAddWeightedAFRealField();
						vcfOutputForTrans.allAddDefaultFormatField("GQ", 40, VCFOutput::FormatEntry("GQ", "1", "Float", "Genotype Quality"), true);
						vcfOutputForTrans.writeOutFixedAndSampleMeta(genomeVcfWithSamples);
					}
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				watch.startNewLap(njh::pasteAsStr("writing translation output - ", varPerTrans.first, " - writeOutSNPsFinalInfo"));
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
			watch.startNewLap(njh::pasteAsStr("writing translation output - writeOutAATypedInfo"));
			translatedRes.writeOutAATypedInfo(njh::files::make_path(variantInfoDir, "seqsAATyped.tab.txt.gz"));
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}

	watch.startNewLap("writing seq output");
	//snps
	uint32_t maxSeqCount = 0;
	auto calcPopMeasuresPars_local =  pars.calcPopMeasuresPars;

	for(auto & varPerChrom : translatedRes.seqVariants_){
		for(const auto & count : varPerChrom.second.depthPerPosition){
			if(count.second > maxSeqCount){
				//cheap way of doing this for now
				maxSeqCount = count.second;
				calcPopMeasuresPars_local.numSegSites_ = varPerChrom.second.getFinalNumberOfSegregatingSites();
			}
		}

		// {
		// 	auto complex_positions = varPerChrom.second.getComplexPositions(pars.variantCallerRunPars.complexVarPars);
		// 	std::cout << varPerChrom.first  << std::endl;
		// 	for(const auto & pos : complex_positions) {
		// 		std::cout << "\t" << pos.start_ << "\t" << pos.size_ << "\t" << pos.count_ << std::endl;
		// 	}
		// 	std::cout << std::endl;
		// }


		{
			//writing vcfs
			auto vcfOutputForChrom = varPerChrom.second.createVCFOutputFixed();
			vcfOutputForChrom.otherHeaderValuePairs_.emplace("reference", pars.transPars.lzPars_.genomeFnp.string());

			vcfOutputForChrom.contigEntries_.emplace(varPerChrom.first, VCFOutput::ContigEntry(varPerChrom.first, contigLengths[varPerChrom.first]));
			vcfOutputForChrom.addDefaultInfoField("TARGET", pars.identifier,VCFOutput::InfoEntry("TARGET", "1", "String", "the target that covers this variant"));
			{
				//auto vcfOutputForChrom = varPerChrom.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf")) );
				// OutputStream vcfOut(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf")));
				// vcfOutputForChrom.writeOutFixedOnly(vcfOut);
			}
			vcfOutputForChrom.headerNonSampleFields_.emplace_back("FORMAT");
			vcfOutputForChrom.formatEntries_.emplace("GT", VCFOutput::FormatEntry(
				"GT", "1", "String",
				"Genotype"
			));
			vcfOutputForChrom.formatEntries_.emplace("DP", VCFOutput::FormatEntry(
				                                         "DP", "1", "Integer",
				                                         "Total Read Depth for this sample, a count of 0 means no coverage in this sample"
			                                         ));
			vcfOutputForChrom.formatEntries_.emplace("AD", VCFOutput::FormatEntry(
				                                         "AD", "R", "Integer",
				                                         "Read Depth for the ref and alt alleles in the order listed, a count of 0 means not detected"
			                                         ));
			vcfOutputForChrom.formatEntries_.emplace("AF", VCFOutput::FormatEntry(
				                                         "AF", "R", "Float",
				                                         "Read Frequency for the ref and alt alleles in the order listed, a freq of 0 means not detected"
			                                         ));

			//add if location is a known drug resistance location
			vcfOutputForChrom.infoEntries_.emplace("AF", VCFOutput::InfoEntry(
				"KNOWN_AA_CHANGE_POS", "0", "Flag",
				"Position is a known chromosome position that is within the codon for a known amino acid position of interest, e.g. drug resistance mutation, etc"
			));

			for (auto & rec : vcfOutputForChrom.records_) {
				//the position in the vcf output is 1 based but the positions saved in knownAAMutsChromPositions are 0 based
				if(!knownAAMutsChromPositions[varPerChrom.first].empty() && njh::in(rec.pos_ - 1, knownAAMutsChromPositions[varPerChrom.first]) ) {
					rec.info_.addMeta("KNOWN_AA_CHANGE_POS",true);
				}
			}

			std::unordered_set<std::string> chromPositions;
			for(const auto & rec : vcfOutputForChrom.records_) {
				chromPositions.emplace(njh::pasteAsStr(rec.chrom_, "-", rec.pos_));
			}
			//key1 = haplotypeName, key2 = chrom, key3 = vcf-based positioning, value = ref,alt
			std::unordered_map<std::string, std::map<std::string, std::map<uint32_t, std::vector<std::pair<std::string, std::string>>>>> vcfAlts;

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
							vcfAlts[seqAlns.first][seqAln.gRegion_.chrom_][vcfPosition].emplace_back(std::make_pair(std::string(1, ref), std::string(1,alt)));
						}
					}
					for(const auto & g : seqAln.comp_.distances_.alignmentGaps_) {
						if(g.second.ref_) {
							//insertion
							//subtract 1 because vcf does insertions/deletions from the base directly proceeding the actual INDEL
							auto realGenomicPos = seqAln.gRegion_.start_ + g.second.refPos_ - 1;
							auto vcfPosition = realGenomicPos + 1;
							auto currentVariantChromPos = njh::pasteAsStr(seqAln.gRegion_.chrom_, "-", vcfPosition);
							if(njh::in(currentVariantChromPos, chromPositions)) {
								auto ref = varPerChrom.second.getBaseForGenomicRegion(realGenomicPos);
								auto alt = njh::pasteAsStr(varPerChrom.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
								vcfAlts[seqAlns.first][seqAln.gRegion_.chrom_][vcfPosition].emplace_back(std::make_pair(std::string(1, ref), alt));
							}
						} else {
							//deletion
							//subtract 1 because vcf does insertions/deletions from the base directly proceeding the actual INDEL
							auto realGenomicPos = seqAln.gRegion_.start_ + g.second.refPos_ - 1;
							auto vcfPosition = realGenomicPos + 1;
							auto currentVariantChromPos = njh::pasteAsStr(seqAln.gRegion_.chrom_, "-", vcfPosition);
							if(njh::in(currentVariantChromPos, chromPositions)) {
								auto ref = njh::pasteAsStr(varPerChrom.second.getBaseForGenomicRegion(realGenomicPos), g.second.gapedSequence_);
								auto alt = varPerChrom.second.getBaseForGenomicRegion(realGenomicPos);
								vcfAlts[seqAlns.first][seqAln.gRegion_.chrom_][vcfPosition].emplace_back(std::make_pair(ref, std::string(1, alt)));
							}
						}
					}
					// std::cout << "seqAlns.first: " << seqAlns.first << std::endl;
					// for(const auto & chrom : vcfAlts[seqAlns.first]) {
					// 	for(const auto & pos : chrom.second) {
					// 		for(const auto & alts : pos.second) {
					// 			std::cout << chrom.first << "\t" << pos.first << "\t" << alts.first << "\t" << alts.second << std::endl;
					// 		}
					// 	}
					// }
					// std::cout << std::endl;
				}
			}
			//add in sample names
			auto allSamplesForVariant = varPerChrom.second.getAllSamples();
			vcfOutputForChrom.samples_ = VecStr(allSamplesForVariant.begin(), allSamplesForVariant.end());
			for (auto&rec: vcfOutputForChrom.records_) {

				for (const auto&sample: allSamplesForVariant) {
					std::vector<uint32_t> dps(1 + rec.alts_.size(), 0);
					// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
					for (const auto& haps: njh::mapAt(samplesToHapsWithReadCnts, sample)) {
						// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
						bool foundAlt = false;
						for (const auto& alt: iter::enumerate(rec.alts_)) {
							if (njh::in(rec.pos_, vcfAlts[haps.first][rec.chrom_])) {
								for (const auto& alts: vcfAlts[haps.first][rec.chrom_][rec.pos_]) {
									if (rec.ref_ == alts.first &&
									    alt.second == alts.second) {
										dps[1 + alt.index] += haps.second;
										foundAlt = true;
									}
								}
							}
						}
						bool coveredByHap = false;
						// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
						if(njh::in(haps.first, translatedRes.seqAlns_)) {
							for(const auto & seqAln : njh::mapAt(translatedRes.seqAlns_, haps.first)) {
								// std::cout << __FILE__ << " : " << __LINE__ << std::endl;
								if(seqAln.gRegion_.chrom_ == rec.chrom_ && rec.pos_ - 1 >= seqAln.gRegion_.start_ && rec.pos_ - 1 < seqAln.gRegion_.end_) {
									coveredByHap = true;
									break;
								}
							}
						}
						if(!foundAlt && coveredByHap) {
							//no alternative found, and reference for this position is covered by hap, increase depth for reference
							dps[0] += haps.second;
						}
					}

					auto dpsSum = vectorSum(dps);
					std::vector<double> dpsFreq;
					dpsFreq.reserve(dps.size());
					for (const auto dp: dps) {
						dpsFreq.emplace_back(dpsSum > 0 ? dp / dpsSum : 0.0);
					}
					if(0 == dpsSum) {
						//if dpsSum equals 0 that most likely means that the sample either had a haplotype that didn't cover
						//this location and/or a haplotype that didn't map
						//if the sample had a mapping haplotype that didn't cover then it either ended before this location or
						//it started after it, either way it's ok to mark as no info for this loc
						rec.sampleFormatInfos_[sample].addMeta("DP", ".");
						rec.sampleFormatInfos_[sample].addMeta("AD", njh::conToStr(std::vector<std::string>(dps.size(), "."), ","));
						rec.sampleFormatInfos_[sample].addMeta("AF", njh::conToStr(std::vector<std::string>(dpsFreq.size(), "."), ","));
					} else {
						rec.sampleFormatInfos_[sample].addMeta("DP", dpsSum);
						rec.sampleFormatInfos_[sample].addMeta("AD", njh::conToStr(dps, ","));
						rec.sampleFormatInfos_[sample].addMeta("AF", njh::conToStr(dpsFreq, ","));
					}
				}
			}
			{
				// OutputStream genomeVcfWithSamples(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomicWithSampleInfo.vcf")));
				OutputStream genomeVcfWithSamples(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf.gz")));
				vcfOutputForChrom.allAddGTFields(pars.variantCallerRunPars.ploidy);
				vcfOutputForChrom.allAutoAddDPFields();
				vcfOutputForChrom.allAutoAddTYPEFields();
				vcfOutputForChrom.allAutoAdd_AN_AC_AF_InfoFields();
				vcfOutputForChrom.allAutoAddWeightedAFRealField();
				vcfOutputForChrom.allAddDefaultFormatField("GQ", 40, VCFOutput::FormatEntry("GQ", "1", "Float", "Genotype Quality"), true);
				vcfOutputForChrom.writeOutFixedAndSampleMeta(genomeVcfWithSamples);

				OutputStream genomeComplexVcfWithSamples(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-complex-genomic.vcf.gz")));

				auto vcfComplexOuptutForChrom = varPerChrom.second.createVCFOutputComplexFixedWithSampleInfo(pars.variantCallerRunPars.ploidy, varPerChrom.first, contigLengths[varPerChrom.first], varPerChrom.second.getComplexPositions(pars.variantCallerRunPars.complexVarPars));
				vcfComplexOuptutForChrom.otherHeaderValuePairs_.emplace("reference", pars.transPars.lzPars_.genomeFnp.string());
				vcfComplexOuptutForChrom.addDefaultInfoField("TARGET", pars.identifier,VCFOutput::InfoEntry("TARGET", "1", "String", "the target that covers this variant"));
				vcfComplexOuptutForChrom.writeOutFixedAndSampleMeta(genomeComplexVcfWithSamples);
			}
		}
		varPerChrom.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt.gz")), varPerChrom.first);
		if(!knownAAMutsChromPositions[varPerChrom.first].empty()){
			//just positions covered
			std::set<uint32_t> coveredPositions;
			for(const auto pos  : knownAAMutsChromPositions[varPerChrom.first]) {
				if(njh::in(pos, varPerChrom.second.allBases)) {
					coveredPositions.insert(pos);
				}
			}
			// varPerChrom.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt.gz")), varPerChrom.first, knownAAMutsChromPositions[varPerChrom.first]);
			varPerChrom.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-knownAA_SNPs.tab.txt.gz")), varPerChrom.first, coveredPositions);
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
			calcPopMeasuresPars_local, alignerObj);
		divMeasures.writeDivMeasures(
			njh::files::make_path(pars.outputDirectory, "divMeasures.tab.txt"),
			inputSeqs, pars.identifier, calcPopMeasuresPars_local);
	}
	{
		auto inputSeqs_prevs = inputSeqs.getPrevalences();
		auto inputSeqs_freqs = inputSeqs.getWeightedAlleleFreqs();
		OutputStream inputSeqs_prev_freq_out(njh::files::make_path(pars.outputDirectory, "seqs_prev_freq.tsv.gz"));
		inputSeqs_prev_freq_out << "target_name\tseq\tfreq\tprev" << std::endl;
		for (const auto & seq : inputSeqs.seqs_) {
			inputSeqs_prev_freq_out << pars.identifier
				<< "\t" << seq->seq_
				<< "\t" << inputSeqs_freqs[seq->seq_]
				<< "\t" << inputSeqs_prevs[seq->seq_] << std::endl;
		}
	}

	if(!pars.metaFieldsToCalcPopDiffs.empty()){
		auto outputDirPerMeta = njh::files::make_path(pars.outputDirectory, "perMetaFields");
		njh::files::makeDir(njh::files::MkdirPar(outputDirPerMeta));
		std::unordered_map<std::string, std::unordered_map<std::string, CollapsedHaps::GenPopMeasuresRes>> measuresPer;

		for(const auto & metaField : pars.metaFieldsToCalcPopDiffs){
			OutputStream divMeasuresOut(njh::files::make_path(outputDirPerMeta, metaField + "_divMeasures.tab.txt.gz"));
			divMeasuresOut << njh::conToStr(calcPopMeasuresPars_local.genHeader(VecStr{metaField}), "\t") << std::endl;
			auto splitSeqs = inputSeqs.splitOutSeqsByMeta(metaField);
			std::unordered_map<std::string, CollapsedHaps::GenPopMeasuresRes> divMeausresPerPop;
			std::unordered_map<std::string, uint32_t> totalHapsPerPop;
			for(const auto & subField : splitSeqs){
				auto calcPopMeasuresPars_local_forSubField = calcPopMeasuresPars_local;
				maxSeqCount = 0;
				for(auto & varPerChrom : translatedRes.seqVariants_) {
					for(const auto & count : varPerChrom.second.depthPerPosition){
						if(count.second > maxSeqCount){
							//hacky way of doing this for now
							maxSeqCount = count.second;
							calcPopMeasuresPars_local_forSubField.numSegSites_ = varPerChrom.second.getFinalNumberOfSegregatingSites(subField.second.getAllSampleNames());
						}
					}
				}


				auto divMeasures = subField.second.getGeneralMeasuresOfDiversity(calcPopMeasuresPars_local_forSubField, alignerObj);
				divMeausresPerPop[subField.first] = divMeasures;
				totalHapsPerPop[subField.first] = subField.second.getTotalHapCount();
				divMeasuresOut << njh::conToStr(divMeasures.getOut(subField.second, njh::pasteAsStr(pars.identifier), calcPopMeasuresPars_local_forSubField, VecStr{subField.first}), "\t") << std::endl;
			}
			OutputStream diffMeasuresOut(njh::files::make_path(outputDirPerMeta, njh::pasteAsStr(metaField, "_diffMeasures.tab.txt.gz")));
			OutputStream pairwiseDiffMeasuresOut(njh::files::make_path(outputDirPerMeta, njh::pasteAsStr(metaField, "_pairwiseDiffMeasures.tab.txt.gz")));
			diffMeasuresOut << "meta"
					<<"\t"<< "metaSubGroupCount"
					<<"\t"<< "target"
					<<"\t"<<"totalHaps"
					<<"\t"<<"uniqueHaps"
					<<"\t"<<"nsamples"
					<<"\t"<<"HsSample"
					<<"\t"<<"HsEst"
					<<"\t"<<"HtSample"
					<<"\t"<<"HtEst"
					<<"\t"<<"Gst"
					<<"\t"<<"GstEst"
					<<"\t"<<"JostD"
					<<"\t"<<"JostDEst"
					<<"\t"<<"ChaoA"
					<<"\t"<<"ChaoB"
					<<"\t"<<"JostDChaoEst"
					<<"\t"<<"In"<< std::endl;

			pairwiseDiffMeasuresOut << "target"
					<< "\t" << metaField << "1"
					<< "\t" << "popMeta" << "1_totalHaps"
					<< "\t" << "popMeta" << "1_uniqueHaps"
					<< "\t" << "popMeta" << "1_samples"
					<< "\t" << "hapsOnlyIn_popMeta" << "1"
					<< "\t" << "hapsOnlyIn_popMeta" << "1CumFreq"
					<< "\t" << metaField << "2"
					<< "\t" << "popMeta" << "2_totalHaps"
					<< "\t" << "popMeta" << "2_uniqueHaps"
					<< "\t" << "popMeta" << "2_samples"
					<< "\t" << "hapsOnlyIn_popMeta" << "2"
					<< "\t" << "hapsOnlyIn_popMeta" << "2CumFreq"
					<< "\t" << "uniqHapsCombinedPops"
					<< "\t" << "uniqHapsSharedInPops"
					<< "\t" << "HsSample"
									<< "\t" << "HsEst"
									<< "\t" << "HtSample"
									<< "\t" << "HtEst"
									<< "\t" << "Gst"
									<< "\t" << "GstEst"
									<< "\t" << "JostD"
									<< "\t" << "JostDEst"
									<< "\t" << "ChaoA"
									<< "\t" << "ChaoB"
									<< "\t" << "JostDChaoEst"
									<< "\t" << "In"

									<< "\t" << "brayCurtisDissim"
									<< "\t" << "brayCurtisRelativeDissim"
									<< "\t" << "jaccardIndexDissim"
									<< "\t" << "sorensenDistance"
									<< "\t" << "RMSE"
									<< "\t" << "correlationDissim"
									<< "\t" << "matchingCoefficientDistance"
									<< "\t" << "plainAvalance"
									<< std::endl;

			auto sampleCountsPerPop = inputSeqs.getSamplesPerSubFieldsForMetaField(metaField);
			auto hapsForTargetPerPopulation = inputSeqs.getHapsPerSampleMetaSubPopulations(metaField);

			if (hapsForTargetPerPopulation.size() > 1) {
				auto generalDiff = PopGenCalculator::getOverallPopDiffWeighted(hapsForTargetPerPopulation);
				diffMeasuresOut << metaField
						<<"\t"<< hapsForTargetPerPopulation.size()
						<<"\t"<< pars.identifier
						<<"\t"<< inputSeqs.getTotalHapCount()
						<<"\t"<< inputSeqs.seqs_.size()
						<<"\t"<< inputSeqs.getAllSampleNames().size()
						<<"\t"<< generalDiff.hsSample_
						<<"\t"<< generalDiff.hsEst_
						<<"\t"<< generalDiff.htSample_
						<<"\t"<< generalDiff.htEst_
						<<"\t"<< generalDiff.gst_
						<<"\t"<< generalDiff.gstEst_
						<<"\t"<< generalDiff.jostD_
						<<"\t"<< generalDiff.jostDEst_
						<<"\t"<< generalDiff.chaoA_
						<<"\t"<< generalDiff.chaoB_
						<<"\t"<< generalDiff.jostDChaoEst_
						<<"\t"<< generalDiff.informativenessForAssign_<< std::endl;
			} else {
				diffMeasuresOut << metaField
						<<"\t"<< hapsForTargetPerPopulation.size()
						<<"\t"<< pars.identifier
						<<"\t"<< inputSeqs.getTotalHapCount()
						<<"\t"<< inputSeqs.seqs_.size()
						<<"\t"<< inputSeqs.getAllSampleNames().size()
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"
						<<"\t"<< "NA"<< std::endl;
			}

			if (hapsForTargetPerPopulation.size() > 1) {
				std::unordered_map<std::string, std::unordered_map<std::string,
					PopGenCalculator::PopDifferentiationMeasuresPairWise>> pairwiseDiffs = PopGenCalculator::getPairwisePopDiffWeighted(hapsForTargetPerPopulation);
				auto keys = getVectorOfMapKeys(pairwiseDiffs);
				njh::sort(keys);
				for(const auto & key : keys){
					auto subKeys = getVectorOfMapKeys(pairwiseDiffs.at(key));
					njh::sort(subKeys);
					for(const auto & subKey : subKeys){
						pairwiseDiffMeasuresOut << pars.identifier
								<< "\t" << key
								<< "\t" << totalHapsPerPop[key]
								<< "\t" << divMeausresPerPop[key].divMeasures_.alleleNumber_
								<< "\t" << sampleCountsPerPop[key].size()
								<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop1_
								<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop1CumFreq_
								<< "\t" << subKey
								<< "\t" << totalHapsPerPop[subKey]
								<< "\t" << divMeausresPerPop[subKey].divMeasures_.alleleNumber_
								<< "\t" << sampleCountsPerPop[subKey].size()
								<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop2_
								<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsInPop2CumFreq_

								<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsAll_
								<< "\t" << pairwiseDiffs.at(key).at(subKey).uniqueHapsShared_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.hsSample_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.hsEst_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.htSample_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.htEst_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.gst_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.gstEst_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostD_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostDEst_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.chaoA_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.chaoB_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.jostDChaoEst_
													<<"\t"<< pairwiseDiffs.at(key).at(subKey).genDiffMeasures_.informativenessForAssign_


													<< "\t" << pairwiseDiffs.at(key).at(subKey).brayCurtisDissim_
													<< "\t" << pairwiseDiffs.at(key).at(subKey).brayCurtisRelativeDissim_
													<< "\t" << pairwiseDiffs.at(key).at(subKey).jaccardIndexDissim_
													<< "\t" << pairwiseDiffs.at(key).at(subKey).sorensenDistance_
													<< "\t" << pairwiseDiffs.at(key).at(subKey).RMSE_
													<< "\t" << pairwiseDiffs.at(key).at(subKey).halfR_
													<< "\t" << pairwiseDiffs.at(key).at(subKey).matchingCoefficientDistance_
													<< "\t" << pairwiseDiffs.at(key).at(subKey).plainAvalance_

													<< std::endl;
					}
				}
			}
		}
	}

	alignerObj->processAlnInfoOutput(pars.alnCacheDir.string(), false);

	//create summary table
	watch.startNewLap("creating summary table");
	{
		// std::cout << __FILE__ << " " << __LINE__ << std::endl;
		OutputStream summaryTable(njh::files::make_path(pars.outputDirectory, "summaryTable.tab.txt.gz"));
		auto genomicLocs = getBeds(njh::files::make_path(variantInfoDir, "inputSeqs.bed"));
		std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> genomicLocationByName;
		for(const auto & genomicLoc : genomicLocs) {
			genomicLocationByName[genomicLoc->name_] = genomicLoc;
		}
		std::set<std::string> transcripts;
		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>>> transcriptionLocationByName;
		std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<seqInfo>>> translatedSeqsByName;
		std::unordered_map<std::string, std::unordered_map<std::string, std::string>> allKnownTyped;
		std::unordered_map<std::string, std::unordered_map<std::string, std::string>> allFullTyped;
		auto transcriptsBedFnp = njh::files::make_path(variantInfoDir, "translatedInput.bed");

		if(bfs::exists(transcriptsBedFnp)) {
			//some targets will intersect with no genes
			auto proteinLocs = getBeds(transcriptsBedFnp);
			for(const auto & proteinLoc : proteinLocs) {
				transcriptionLocationByName[proteinLoc->chrom_][proteinLoc->name_] = proteinLoc;
				if("*" != proteinLoc->chrom_) {
					transcripts.emplace(proteinLoc->chrom_);
				}
			}
			SeqInput translationReader(SeqIOOptions::genFastaInGz(njh::files::make_path(variantInfoDir, "translatedInput.fasta.gz")));
			auto translatedSeqs = translationReader.readAllReadsPtrs<seqInfo>();
			for(const auto & translatedSeq : translatedSeqs) {
				MetaDataInName meta(translatedSeq->name_);
				auto transcript = meta.getMeta("transcript");
				meta.removeMeta("transcript");
				meta.resetMetaInName(translatedSeq->name_);
				translatedSeqsByName[transcript][translatedSeq->name_] = translatedSeq;
				transcripts.emplace(transcript);
			}

			table seqsAATyped(njh::files::make_path(variantInfoDir, "seqsAATyped.tab.txt.gz"), "\t", true);
			for(const auto & row : seqsAATyped) {
				if(njh::notIn(row[seqsAATyped.getColPos("fullTyped")], VecStr{"Untranslatable", "Unmappable"})) {
					// std::cout << __FILE__ << " " << __LINE__ << std::endl;
					// std::cout << njh::conToStr(row, ",") << std::endl;
					if(!row[seqsAATyped.getColPos("knownTyped")].empty()) {
						auto knownTyped = tokenizeString(row[seqsAATyped.getColPos("knownTyped")], "--");
						allKnownTyped[knownTyped[0]][row[seqsAATyped.getColPos("name")]] = knownTyped[1];
						// std::cout << "knownTyped: " << njh::conToStr(knownTyped, ",")  << std::endl;
					}
					if(!row[seqsAATyped.getColPos("fullTyped")].empty()) {
						auto fullTyped = tokenizeString( row[seqsAATyped.getColPos("fullTyped")], "--");
						allFullTyped[fullTyped[0]][row[seqsAATyped.getColPos("name")]] = fullTyped[1];
						// std::cout << "fullTyped: " << njh::conToStr(fullTyped, ",") << std::endl;
						transcripts.emplace(fullTyped[0]);
					}
				}
			}
		} else if (translatedRes.seqsTranslationFiltered_.size() == inputSeqs.size()) {
			transcripts.emplace("untranslatable");
		} else {
			transcripts.emplace("intergenic");
		}

		auto metaTab = inputSeqs.createMetaFieldsTable(true);
		metaTab.addColumn({pars.identifier}, "target");

		summaryTable << njh::conToStr(metaTab.columnNames_, "\t");
		summaryTable << "\t" << "chrom" << "\t" << "0based_start" << "\t" << "0based_end" << "\t" << "length" << "\t" << "strand";
		summaryTable << "\t" << "transcript";
		summaryTable << "\t" << "translatedSeq";
		summaryTable << "\t" << "transcript_1based_start" << "\t" << "transcript_1based_end" << "\t" << "transcript_length";
		summaryTable << "\t" << "transcript_knownAATyped" << "\t" << "transcript_fullAATyped";
		summaryTable << std::endl;
		for(const auto & row : metaTab) {
			for(const auto & transcript : transcripts) {
				auto seqName = row[metaTab.getColPos("CollapsedName")];

				summaryTable << njh::conToStr(row, "\t");
				if("*" == genomicLocationByName[seqName]->chrom_) {
					summaryTable << "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA";
				} else {
					summaryTable << "\t" << genomicLocationByName[seqName]->chrom_
							<< "\t" << genomicLocationByName[seqName]->chromStart_
							<< "\t" << genomicLocationByName[seqName]->chromEnd_
							<< "\t" << genomicLocationByName[seqName]->length()
							<< "\t" << genomicLocationByName[seqName]->strand_;
				}

				if("*" == genomicLocationByName[seqName]->chrom_) {
					summaryTable << "\t" << "unmappable";
				} else if(njh::notIn(seqName, translatedSeqsByName[transcript]) && "intergenic" != transcript) {
					summaryTable << "\t" << "untranslatable";
				} else {
					summaryTable << "\t" << transcript;
				}

				if(njh::in(seqName, translatedSeqsByName[transcript])) {
					summaryTable << "\t" << translatedSeqsByName[transcript][seqName]->seq_
							<< "\t" << transcriptionLocationByName[transcript][seqName]->chromStart_ + 1
							<< "\t" << transcriptionLocationByName[transcript][seqName]->chromEnd_
							<< "\t" << transcriptionLocationByName[transcript][seqName]->length()
							<< "\t" << (allKnownTyped[transcript][seqName].empty() ? "None" : allKnownTyped[transcript][seqName])
							<< "\t" << (allFullTyped[transcript][seqName].empty() ? "None" : allFullTyped[transcript][seqName]);
				} else {
					summaryTable
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA"
							<< "\t" << "NA";
				}
				summaryTable << std::endl;
			}
		}
	}
	// std::cout << __FILE__ << " " << __LINE__ << std::endl;
	OutputStream timeLog(njh::files::make_path(pars.outputDirectory, "timeLog.txt"));
	watch.logLapTimes(timeLog,true, 6, true);
	return translatedRes;
}




}  // namespace njhseq
