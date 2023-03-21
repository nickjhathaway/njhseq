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

	uint64_t maxLen = readVec::getMaxLength(inputSeqs.seqs_);

	std::shared_ptr<aligner> alignerObj = std::make_shared<aligner>(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
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

  OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "inputSeqs.bed"));
	translatedRes.writeOutSeqAlnIndvVars(njh::files::make_path(variantInfoDir, "variantsPerSeqAln.tab.txt.gz"));
	translatedRes.writeSeqLocations(popBedLocs);

	std::unordered_map<std::string, std::set<uint32_t>> knownAAMutsChromPositions;

	OutputStream seqsUnableToBeMappedOut(njh::files::make_path(variantInfoDir, "seqsUnableToBeMapped.txt"));
	seqsUnableToBeMappedOut << njh::conToStr(translatedRes.seqsUnableToBeMapped_, "\n") << std::endl;
	OutputStream seqsTranslationFilteredOut(njh::files::make_path(variantInfoDir, "seqsTranslationFiltered.txt"));
	seqsTranslationFilteredOut << njh::conToStr(translatedRes.seqsTranslationFiltered_, "\n") << std::endl;


	if(!translatedRes.translations_.empty()){
		SeqOutput transwriter(SeqIOOptions::genFastaOutGz(njh::files::make_path(variantInfoDir, "translatedInput.fasta.gz")));

		std::unordered_map<std::string, std::vector<seqInfo>> translatedSeqsByTranscript;
		std::unordered_map<std::string, std::vector<std::unordered_set<std::string>>> translatedSeqInputNames;
		auto seqNames = njh::getVecOfMapKeys(translatedRes.translations_);
		njh::sort(seqNames);
		for(const auto & seqName : seqNames){
			for(const auto & transcript : translatedRes.translations_.at(seqName)){
				transwriter.openWrite(transcript.second.translation_);
				translatedSeqsByTranscript[transcript.first].emplace_back(transcript.second.translation_);
				translatedSeqsByTranscript[transcript.first].back().cnt_ = inputSeqs.names_[seqNameKey[seqName]].size();
				translatedSeqInputNames[transcript.first].emplace_back(inputSeqs.names_[seqNameKey[seqName]]);
			}
		}
		OutputStream divMeasuresOut(njh::files::make_path(variantInfoDir, "translatedDivMeasures.tab.txt"));
		divMeasuresOut << njh::conToStr(pars.calcPopMeasuresPars.genHeader(), "\t") << std::endl;

		///
		auto fullTypedAAForTranslated = translatedRes.translated_genAATypedStr();
		auto knownsTypedAAForTranslated = translatedRes.translated_genAATypedStrOnlyKnowns();
		auto variableTypedAAForTranslated = translatedRes.translated_genAATypedStrOnlyPopVariant();

		///
		for(const auto & translatedSeqs : translatedSeqsByTranscript){
			auto inputTranslatedSeq = CollapsedHaps::collapseReads(translatedSeqs.second, translatedSeqInputNames[translatedSeqs.first]);
			std::string identifierTranslated = njh::pasteAsStr(pars.identifier, "-translated");
			if(translatedSeqsByTranscript.size() > 1){
				identifierTranslated = njh::pasteAsStr(pars.identifier, "-", translatedSeqs.first, "-translated");
			}

//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			auto divMeasures = inputTranslatedSeq.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			divMeasuresOut << njh::conToStr(divMeasures.getOut(inputTranslatedSeq, identifierTranslated, calcPopMeasuresPars), "\t")  << std::endl;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
		OutputStream transBedLocs(njh::files::make_path(variantInfoDir, "translatedInput.bed"));
		translatedRes.writeSeqLocationsTranslation(transBedLocs);
		translatedRes.writeOutTranslatedIndvVars(njh::files::make_path(variantInfoDir, "variantsPerTranslatedSeq.tab.txt.gz"), translator->knownAminoAcidPositions_);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			//protein
			for(auto & varPerTrans : translatedRes.proteinVariants_){
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				varPerTrans.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein.vcf")));
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
		varPerChrom.second.writeVCF(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-genomic.vcf")));
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
				auto divMeasures = subField.second.getGeneralMeasuresOfDiversity(calcPopMeasuresPars, alignerObj);
				divMeasuresOut << njh::conToStr(divMeasures.getOut(subField.second, njh::pasteAsStr(pars.identifier, "::", subField.first), calcPopMeasuresPars), "\t") << std::endl;
			}
		}
	}

	alignerObj->processAlnInfoOutput(pars.alnCacheDir.string(), false);

	return translatedRes;
}




}  // namespace njhseq
