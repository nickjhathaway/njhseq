#pragma once

/*
 * PopGenUtils.hpp
 *
 *  Created on: Oct 16, 2021
 *      Author: nick
 */




#include "njhseq/objects/seqContainers/CollapsedHaps.hpp"
#include "njhseq/PopulationGenetics/PopGenCalcs.hpp"
#include "njhseq/PopulationGenetics/HapsEncodedMatrix.hpp"
#include "njhseq/objects/Gene/TranslatorByAlignment.hpp"
#include "njhseq/concurrency/PairwisePairFactory.hpp"



namespace njhseq {

struct CollapseAndCallVariantsPars{
	CollapseAndCallVariantsPars(){
		variantCallerRunPars.lowVariantCutOff = 0.005;
		variantCallerRunPars.occurrenceCutOff = 2;
	}

	SeqIOOptions inOpts;
	bfs::path outputDirectory;
	bool overWriteDirectory{false};

  TranslatorByAlignment::TranslatorByAlignmentPars transPars;
  TranslatorByAlignment::RunPars variantCallerRunPars;
  CollapsedHaps::GenPopMeasuresPar calcPopMeasuresPars;


	std::string identifier = "";
	GenomicRegion refSeqRegion;

	uint32_t numThreads = 1;
	bool noDiagAlnPairwiseComps = false;


	bfs::path metaFnp;
	std::set<std::string> ignoreSubFields{"site:LabCross", "site:LabControl", "site:LabContaminated"};

	VecStr metaFieldsToCalcPopDiffs;

	bfs::path alnCacheDir;
};


/*
 *
  pars.inOpts = setUp.pars_.ioOptions_;



	setUp.processAlnInfoInput();
	pars.alnCacheDir = setUp.pars_.alnInfoDirName_;

	setUp.setOption(pars.numThreads, "--numThreads", "number of threads");
	pars.calcPopMeasuresPars.numThreads = pars.numThreads;
	setUp.setOption(pars.calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(pars.noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	pars.calcPopMeasuresPars.diagAlnPairwiseComps = !noDiagAlnPairwiseComps;

	setUp.setOption(pars.variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(pars.variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	pars.calcPopMeasuresPars.lowVarFreq = pars.variantCallerRunPars.lowVariantCutOff;
	pars.transPars.setOptions(setUp, true);

	setUp.setOption(pars.metaFnp,    "--metaFnp",    "Meta data to add to sequences");
	setUp.setOption(pars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");

	setUp.setOption(pars.identifier, "--identifier", "Give a identifier name for info", true);
	setUp.setOption(pars.metaFieldsToCalcPopDiffs, "--metaFieldsToCalcPopDiffs", "Meta Fields To Calculate Pop Diffs");
	setUp.setOption(pars.outputDirectory, "--dout", "Output directory", true);
	setUp.setOption(pars.overWriteDirectory, "--overWriteDir", "over write output Directory");

 */

TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars);
TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars, const std::vector<seqInfo> & input);
TranslatorByAlignment::TranslatorByAlignmentResult collapseAndCallVariants(const CollapseAndCallVariantsPars & pars, CollapsedHaps & inputSeqs);





}  // namespace njhseq

