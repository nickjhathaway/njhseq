//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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
#include <bibcpp/bashUtils.h>
#include <bibcpp/files/fileUtilities.hpp>

#include "seqSetUp.hpp"
#include "bibseq/IO/SeqIO/SeqInput.hpp"
namespace bibseq {



const VecStr seqSetUp::readInFormatsAvailable_ {"sff", "sffBin", "fasta", "fastq",
	"bam", "fastqgz", "fastagz", "fastq1", "fastq2", "fastq1gz", "fastq2gz"};


void seqSetUp::processComparison(comparison & comp, std::string stub) {
	if("" != stub && '-' != stub.back()){
		stub += "-";
	}
	setOption(comp.oneBaseIndel_, "--" + stub + "oneBaseIndel", "Allowable one base indels");
	setOption(comp.twoBaseIndel_, "--" + stub + "twoBaseIndel", "Allowable two base indels");
	setOption(comp.largeBaseIndel_, "--" + stub + "largeBaseIndel",
			"Allowable large base indels (>2 bases)");
	setOption(comp.lqMismatches_, "--" + stub + "lqMismatches",
			"Allowable low quality mismatches");
	setOption(comp.hqMismatches_, "--" + stub + "hqMismatches",
			"Allowable high quality mismatches");
}


void seqSetUp::processQualityFiltering() {

	setOption(pars_.qFilPars_.qualWindow_, "--qualWindow",
				"Sliding Quality Window, format is WindowSize,WindowStep,Threshold");
	seqUtil::processQualityWindowString(pars_.qFilPars_.qualWindow_, pars_.qFilPars_.qualityWindowLength_,
			pars_.qFilPars_.qualityWindowStep_, pars_.qFilPars_.qualityWindowThres_);
	setOption(pars_.qFilPars_.qualCheck_, "--qualCheck", "Qual Check Level");
	setOption(pars_.qFilPars_.qualCheckCutOff_, "--qualCheckCutOff",
			"Cut Off for fraction of bases above qual check of "
					+ estd::to_string(pars_.qFilPars_.qualCheck_));
}

void seqSetUp::processClusteringOptions(){
  processSkipOnNucComp();
  processAdjustHRuns();
  bool firstMatch = false;
  setOption(firstMatch,     "--firstMatch",     "Settle for first Match in Clustering");
  pars_.colOpts_.bestMatchOpts_.findingBestMatch_ = !firstMatch;
  setOption(pars_.colOpts_.bestMatchOpts_.bestMatchCheck_, "--bestMatchCheck", "Best Match Check Number");
}


void seqSetUp::processGap() {
	// check command line for gap settings
	if (setOption(pars_.gap_, "--gapAll",
			"Gap Penalties for All (middle and end gaps)", false, "Alignment")) {
		pars_.gapInfo_ = gapScoringParameters(pars_.gap_);
		pars_.gapLeft_ = pars_.gap_;
		pars_.gapRight_ = pars_.gap_;
	}
	if (setOption(pars_.gap_, "--gap", "Gap Penalties for Middle Gap", false,
			"Alignment")) {
		pars_.gapInfo_.processGapStr(pars_.gap_, pars_.gapInfo_.gapOpen_,
				pars_.gapInfo_.gapExtend_);
	}
	if (setOption(pars_.gapLeft_, "--gapLeft", "Gap Penalties for Left End Gap",
			false, "Alignment")) {
		pars_.gapInfo_.processGapStr(pars_.gapLeft_, pars_.gapInfo_.gapLeftOpen_,
				pars_.gapInfo_.gapLeftExtend_);
	}
	if (setOption(pars_.gapRight_, "--gapRight",
			"Gap Penalties for Right End Gap", false, "Alignment")) {
		pars_.gapInfo_.processGapStr(pars_.gapRight_, pars_.gapInfo_.gapRightOpen_,
				pars_.gapInfo_.gapRightExtend_);
	}

	pars_.gapInfo_.setIdentifer();
}

void seqSetUp::processGapRef() {
	// check command line for gap settings
	if (setOption(pars_.gapRef_, "--refGapAll", "Gap Penalties for Ref All")) {
		pars_.gapInfoRef_ = gapScoringParameters(pars_.gapRef_);
		pars_.gapLeftRef_ = pars_.gapRef_;
		pars_.gapRightRef_ = pars_.gapRef_;
	} else {
		setOption(pars_.gapRef_, "--refGap", "Gap Penalties for Ref Middle Gap");
		setOption(pars_.gapLeftRef_, "--refGapLeft", "Gap Penalties for Ref Left End Gap");
		setOption(pars_.gapRightRef_, "--refGapRight",
				"Gap Penalties for Ref Right End Gap");
		// get the gap penalty
		pars_.gapInfoRef_.processGapStr(pars_.gapRef_, pars_.gapInfoRef_.gapOpen_,
				pars_.gapInfoRef_.gapExtend_);
		pars_.gapInfoRef_.processGapStr(pars_.gapLeftRef_, pars_.gapInfoRef_.gapLeftOpen_,
				pars_.gapInfoRef_.gapLeftExtend_);
		pars_.gapInfoRef_.processGapStr(pars_.gapRightRef_, pars_.gapInfoRef_.gapRightOpen_,
				pars_.gapInfoRef_.gapRightExtend_);
	}
	pars_.gapInfoRef_.setIdentifer();
}

void seqSetUp::processQualThres() {
	// check command line for qualThres setting
	setOption(pars_.qualThres_, "--qualThres",
			"Quality Thresholds, should go PrimaryQual,SecondaryQual", false, "Mismatch Quality Determination");
	// get the qualities
	auto qualToks = tokenizeString(pars_.qualThres_, ",");
	if (qualToks.size() != 2) {
		throw std::runtime_error {
				"QaulThres should be two numbers separated by a comma, eg 20,15, not "
						+ pars_.qualThres_ };
	}
	pars_.qScorePars_.primaryQual_ = std::stoi(qualToks[0]);
	pars_.qScorePars_.secondaryQual_ = std::stoi(qualToks[1]);

	setOption(pars_.qScorePars_.qualThresWindow_, "--qualThresWindow",
			"Quality Threshold Window Length", false, "Mismatch Quality Determination");
}

CollapseIterations seqSetUp::processIteratorMap(const bfs::path & parameters) {

	if (!bfs::exists(parameters)) {
		failed_ = true;
		warnings_.emplace_back(
				bib::bashCT::red + bib::bashCT::bold + "File " + parameters.string()
						+ " doesn't exist" + bib::bashCT::reset);
	}
	//std::cout << __PRETTY_FUNCTION__ << 1 << std::endl;
	//CollapseIterations ret(parameters, false);
	//std::cout << __PRETTY_FUNCTION__ << 2 << std::endl;
	return {parameters.string(), false};
}

CollapseIterations seqSetUp::processIteratorMapOnPerId(const bfs::path & parameters) {

	if (!bfs::exists(parameters)) {
		failed_ = true;
		warnings_.emplace_back(
				bib::bashCT::red + bib::bashCT::bold + "File " + parameters.string()
						+ " doesn't exist" + bib::bashCT::reset);
	}
	//std::cout << __PRETTY_FUNCTION__ << 1 << std::endl;
	//CollapseIterations ret();
	//std::cout << __PRETTY_FUNCTION__ << 2 << std::endl;
	return {parameters.string(), true};
}

void seqSetUp::processKmerLenOptions(){
  setOption(pars_.colOpts_.kmerOpts_.kLength_, "--kLength", "Kmer Length");
}

void seqSetUp::processKmerProfilingOptions() {
	setOption(pars_.colOpts_.kmerOpts_.runCutOffString_, "--runCutOff",
			"Kmer frequency cut off for a mismatch to be considered low frequency",
			false, "Low Frequency Mismatch Determination");
	setOption(pars_.colOpts_.kmerOpts_.kLength_, "--kLength", "Kmer Length",
			false, "Low Frequency Mismatch Determination");
	bool forKmerProfiling = true;
	if (forKmerProfiling && pars_.colOpts_.kmerOpts_.kLength_ % 2 == 0) {
		pars_.colOpts_.kmerOpts_.kLength_--;
		warnings_.emplace_back(
				"--kLength needs to be odd, not even, changing to "
						+ estd::to_string(pars_.colOpts_.kmerOpts_.kLength_));
	}
	bool kAnywhere = false;
	setOption(kAnywhere, "--kAnywhere",
			"Count Kmers without regard for position",
			false, "Low Frequency Mismatch Determination");
	pars_.colOpts_.kmerOpts_.kmersByPosition_ = !kAnywhere;
	setOption(pars_.expandKmerPos_, "--expandKmerPos",
			"Expand Kmer Position Found At",
			false, "Low Frequency Mismatch Determination");
	setOption(pars_.expandKmerSize_, "--expandKmerSize",
			"Expand Kmer Size for extending where kmers where found",
			false, "Low Frequency Mismatch Determination");
}

void seqSetUp::processScoringPars() {
	setOption(pars_.local_, "--local", "Local alignment", false, "Alignment");
	setOption(pars_.colOpts_.alignOpts_.countEndGaps_, "--countEndGaps", "Count End Gaps", false, "Alignment");
	bool noHomopolymerWeighting = false;
	setOption(noHomopolymerWeighting, "--noHomopolymerWeighting",
			"Don't do Homopolymer Weighting", false, "Alignment");
	pars_.colOpts_.iTOpts_.weighHomopolyer_ = !noHomopolymerWeighting;
	std::string scoreMatrixFilename = "";
	bool degenScoring = false;
	bool caseInsensitive = false;
	bool lessN = false;
	if (setOption(scoreMatrixFilename, "--scoreMatrix", "Score Matrix Filename", false, "Alignment")) {
		pars_.scoring_ = substituteMatrix(scoreMatrixFilename);
	} else {
		setOption(pars_.generalMatch_, "--match", "Match score for alignment", false, "Alignment");
		setOption(pars_.generalMismatch_, "--mismatch",
				"Mismatch score for alignment", false, "Alignment");
		setOption(degenScoring, "--degen", "Use Degenerative Base Scoring for alignment", false, "Alignment");
		setOption(lessN, "--lessN", "Use Degenerative Base Scoring but use a lesser score for the degenerative bases", false, "Alignment");
		setOption(caseInsensitive, "--caseInsensitive",
				"Use Case Insensitive Scoring for alignment", false, "Alignment");
		pars_.scoring_ = substituteMatrix::createScoreMatrix(pars_.generalMatch_, pars_.generalMismatch_, degenScoring, lessN, caseInsensitive);
	}
}


void seqSetUp::processSkipOnNucComp(){
  setOption(pars_.colOpts_.skipOpts_.skipOnLetterCounterDifference_, "--fastClustering",
                  "Skip comparisons if their nucleotide composition differs", false, "Clustering");
  setOption(pars_.colOpts_.skipOpts_.fractionDifferenceCutOff_, "--nucCutOff",
                  "Fraction Difference in nucleotide composition cut off for when --fastClustering is used", false, "Clustering");
}

void seqSetUp::processAdjustHRuns(){
  setOption(pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_,
            "--adjustHomopolyerRuns",
            "Adjust Homopolyer Runs To Be Same Quality, useful for IonTorrent");
}


bool seqSetUp::processReadInNames(bool required) {
	return processReadInNames(readInFormatsAvailable_, required);
}
bool seqSetUp::processReadInNames(const VecStr & formats, bool required) {
	setOption(pars_.ioOptions_.processed_, "--processed",
			"Processed, Input Sequence Name has a suffix that contains abundance info", false, "Reading Sequence Input");
	setOption(pars_.ioOptions_.lowerCaseBases_, "--lower",
			"How to handle Lower Case Bases", false, "Reading Sequence Input");
	setOption(pars_.ioOptions_.removeGaps_, "--removeGaps",
			"Remove Gaps from Input Sequences", false, "Reading Sequence Input");
	bool noWhiteSpace = false;
	setOption(noWhiteSpace, "--trimAtWhiteSpace",
			"Remove everything after first whitespace character in input sequence name",
			false, "Reading Sequence Input");
	pars_.ioOptions_.includeWhiteSpaceInName_ = !noWhiteSpace;

	std::stringstream formatWarnings;
	bool foundUnrecFormat = false;

	auto formatChecker = [](const std::string & conVal, const std::string & inVal) {
		return bib::strToLowerRet(bib::lstripRet(conVal, '-')) == bib::strToLowerRet(bib::lstripRet(inVal, '-'));
	};

	for (const auto & format : formats) {
		if (!bib::has(readInFormatsAvailable_, format, formatChecker)) {
			addWarning(
					"Format: " + format + " is not an available sequence input format in "
							+ std::string(__PRETTY_FUNCTION__));
			foundUnrecFormat = true;
		}
	}
	if (foundUnrecFormat) {
		failed_ = true;
		return false;
	}
	VecStr readInFormatsFound;
	if (commands_.gettingFlags()
			|| commands_.printingHelp()
			|| commands_.gettingVersion()) {
		bib::progutils::Flag sffFlagOptions(pars_.ioOptions_.firstName_, "--sff",
				"Input sequence filename, only need 1, sff text file", required, "Reading Sequence Input");
		bib::progutils::Flag sffBinFlagOptions(pars_.ioOptions_.firstName_,
				"--sffBin", "Input sequence filename, only need 1, sff binary file", required, "Reading Sequence Input");
		bib::progutils::Flag fastaFlagOptions(pars_.ioOptions_.firstName_, "--fasta",
				"Input sequence filename, only need 1, fasta text file", required, "Reading Sequence Input");
		bib::progutils::Flag fastqFlagOptions(pars_.ioOptions_.firstName_, "--fastq",
				"Input sequence filename, only need 1, fastq text file", required, "Reading Sequence Input");
		bib::progutils::Flag fastqFirstMateFlagOptions(pars_.ioOptions_.firstName_, "--fastq1",
				"Input sequence filename, only need 1, fastq first mate text file", required, "Reading Sequence Input");
		bib::progutils::Flag fastqSecondMateFlagOptions(pars_.ioOptions_.firstName_, "--fastq2",
				"Input sequence filename, only needed with --fastq1, fastq second mate text file", required, "Reading Sequence Input");
		bib::progutils::Flag fastqFirstMateGzFlagOptions(pars_.ioOptions_.firstName_, "--fastq1gz",
				"Input sequence filename, only need 1, fastq first mate gzipped file", required, "Reading Sequence Input");
		bib::progutils::Flag fastqSecondMateGzFlagOptions(pars_.ioOptions_.firstName_, "--fastq2gz",
				"Input sequence filename, only needed with --fastq1gz, fastq second mate gzipped file", required, "Reading Sequence Input");
		/*bib::progutils::Flag fastqComplSecondMateFlagOptions(pars_.ioOptions_.revComplMate_, "--complementMate",
						"Complement second mate in paired reads", false);
		*/
		bib::progutils::Flag fastqgzFlagOptions(pars_.ioOptions_.firstName_,
				"--fastqgz", "Input sequence filename, only need 1, fastq gzipped file", required, "Reading Sequence Input");
		bib::progutils::Flag fastagzFlagOptions(pars_.ioOptions_.firstName_,
				"--fastagz", "Input sequence filename, only need 1, fasta gzipped file", required, "Reading Sequence Input");
		bib::progutils::Flag bamFlagOptions(pars_.ioOptions_.firstName_, "--bam",
				"Input sequence filename, only need 1, bam file", required, "Reading Sequence Input");
		bib::progutils::flagHolder seqReadInFlags;
		seqReadInFlags.addFlag(sffFlagOptions);
		seqReadInFlags.addFlag(sffBinFlagOptions);
		seqReadInFlags.addFlag(fastaFlagOptions);
		seqReadInFlags.addFlag(fastqFlagOptions);
		seqReadInFlags.addFlag(bamFlagOptions);
		seqReadInFlags.addFlag(fastqgzFlagOptions);
		seqReadInFlags.addFlag(fastagzFlagOptions);
		seqReadInFlags.addFlag(fastqFirstMateFlagOptions);
		seqReadInFlags.addFlag(fastqSecondMateFlagOptions);
		seqReadInFlags.addFlag(fastqFirstMateGzFlagOptions);
		seqReadInFlags.addFlag(fastqSecondMateGzFlagOptions);

		for (const auto & formatFlag : seqReadInFlags.flags_) {
			if (bib::has(formats, formatFlag.second.flags_.front(), formatChecker)) {
				flags_.addFlag(formatFlag.second);
			}
		}
		std::string fastq1Flag = "--fastq1";
		std::string fastq2Flag = "--fastq2";

		std::string fastq1GzFlag = "--fastq1gz";
		std::string fastq2GzFlag = "--fastq2gz";

		if(bib::has(formats, fastq1Flag, formatChecker) ){
			//flags_.addFlag(fastqComplSecondMateFlagOptions);
			if(!bib::has(formats, fastq2Flag, formatChecker)){
				flags_.addFlag(fastqSecondMateFlagOptions);
			}
		}
		if(bib::has(formats, fastq1GzFlag, formatChecker) ){
			//flags_.addFlag(fastqComplSecondMateFlagOptions);
			if(!bib::has(formats, fastq2GzFlag, formatChecker)){
				flags_.addFlag(fastqSecondMateGzFlagOptions);
			}
		}
	}
	//compPerCutOff
	//process format information
	//hasFlagCaseInsen(
	if(commands_.hasFlagCaseInsenNoDash("--fasta")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTA;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTA;
		readInFormatsFound.emplace_back("--fasta");
	}
	if(commands_.hasFlagCaseInsenNoDash("--sff")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::SFFTXT;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
		readInFormatsFound.emplace_back("--sff");
	}
	if(commands_.hasFlagCaseInsenNoDash("--sffBin")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::SFFBIN;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
		readInFormatsFound.emplace_back("--sffBin");
	}
	if(commands_.hasFlagCaseInsenNoDash("--bam")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::BAM;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
		readInFormatsFound.emplace_back("--bam");
	}
	if(commands_.hasFlagCaseInsenNoDash("--fastq")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQ;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
		readInFormatsFound.emplace_back("--fastq");
	}
	if(commands_.hasFlagCaseInsenNoDash("--fastq1")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIRED;
		readInFormatsFound.emplace_back("--fastq1");
	}
	if(commands_.hasFlagCaseInsenNoDash("--fastq1gz")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIREDGZ;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
		readInFormatsFound.emplace_back("--fastq1gz");
	}
	if(commands_.hasFlagCaseInsenNoDash("--fastqgz")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQGZ;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQGZ;
		readInFormatsFound.emplace_back("--fastqgz");
	}
	if(commands_.hasFlagCaseInsenNoDash("--fastagz")){
		pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTAGZ;
		pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTAGZ;
		readInFormatsFound.emplace_back("--fastagz");
	}
	if(readInFormatsFound.size() > 1){
    std::stringstream tempOut;
    tempOut << bib::bashCT::bold
    		<< "Found multiple read in options, should only have one"
            << std::endl;
    tempOut << vectorToString(readInFormatsFound, ",") + bib::bashCT::reset;
    tempOut << std::endl;
    addOtherVec(warnings_, streamToVecStr(tempOut));
    failed_ = true;
    return false;
	} else if (readInFormatsFound.empty()){
		if(required){
	    std::stringstream tempOut;
	    tempOut << bib::bashCT::bold
	    				<< "Did not find a recognizable read in option"
	            << std::endl;
	    tempOut << "Options include: "
	            << bib::bashCT::red << bib::conToStr(formats, ",")
	            << bib::bashCT::reset << std::endl;
	    tempOut << "Command line arguments" << std::endl;
	    writeOutCommandLineArguments(commands_.arguments_, tempOut);
	    addOtherVec(warnings_, streamToVecStr(tempOut));
	    failed_ = true;
		}
    return false;
	} else {
		if (!commands_.gettingFlags()
				&& !commands_.printingHelp()
				&& !commands_.gettingVersion()){
			setOption(pars_.ioOptions_.firstName_, readInFormatsFound.front(), "In Sequence Filename");
		}
		if(readInFormatsFound.front() == "--fasta"){
			if(setOption(pars_.ioOptions_.secondName_, "--qual", "Name of the quality file")){
				pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTAQUAL;
				pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
			}
		}
		if (readInFormatsFound.front() == "--fastq1") {
			if (!setOption(pars_.ioOptions_.secondName_, "--fastq2",
					"Name of the mate file")) {
				addWarning("If supplying -fastq1 need to also have -fastq2");
				failed_ = true;
			}
			/*
			setOption(pars_.ioOptions_.revComplMate_, "--complementMate",
					"Whether to complement the sequence in the mate file");*/
		}
		if (readInFormatsFound.front() == "--fastq1gz") {
			if (!setOption(pars_.ioOptions_.secondName_, "--fastq2gz",
					"Name of the mate file")) {
				addWarning("If supplying -fastq1gz need to also have -fastq2gz");
				failed_ = true;
			}
			/*
			setOption(pars_.ioOptions_.revComplMate_, "--complementMate",
					"Whether to complement the sequence in the mate file");*/
		}
		//setOption(pars_.ioOptions_.out_.outFilename_, "--out", "Name of the out sequence file");
		return true;
	}
}

void seqSetUp::processDirectoryOutputName(const std::string& defaultName,
		bool mustMakeDirectory) {
	//std::cout << defaultName << std::endl;
	setOption(pars_.overWriteDir_, "--overWriteDir",
			"If the directory already exists over write it",false, "Output Directory");
	pars_.directoryName_ = "./";
	if (setOption(pars_.directoryName_, "--dout", "Output Directory Name", false, "Output Directory")) {
		if (!failed_) {
			//std::cout << pars_.directoryName_ << std::endl;
			std::string newDirectoryName = bib::replaceString(pars_.directoryName_,
					"TODAY", getCurrentDate()) + "/";
			if (bib::files::bfs::exists(newDirectoryName) && pars_.overWriteDir_) {
				bib::files::rmDirForce(newDirectoryName);
			} else if (bib::files::bfs::exists(newDirectoryName)
					&& !pars_.overWriteDir_) {
				failed_ = true;
				addWarning(
						"Directory: " + newDirectoryName
								+ " already exists, use --overWriteDir to over write it");
			}
			bib::files::makeDirP(
					bib::files::MkdirPar(newDirectoryName, pars_.overWriteDir_));
			pars_.directoryName_ = newDirectoryName;
		}
	} else {
		if (mustMakeDirectory && !failed_) {
			std::string newDirectoryName = bib::replaceString(defaultName, "TODAY",
					getCurrentDate()) + "/";
			if (bib::files::bfs::exists(newDirectoryName) && pars_.overWriteDir_) {
				bib::files::rmDirForce(newDirectoryName);
			} else if (bib::files::bfs::exists(newDirectoryName)
					&& !pars_.overWriteDir_) {
				failed_ = true;
				addWarning(
						"Directory: " + newDirectoryName
								+ " already exists, use --overWriteDir to over write it");
			}
			pars_.directoryName_ = bib::files::makeDir("./",
					bib::files::MkdirPar(defaultName, pars_.overWriteDir_)).string();
		}
	}
}

void seqSetUp::processDirectoryOutputName(bool mustMakeDirectory) {
  std::string seqName = bfs::basename(pars_.ioOptions_.firstName_) + "_" +
  		bib::replaceString(commands_.getProgramName(), " ", "-") + "_" + getCurrentDate();
  processDirectoryOutputName(seqName, mustMakeDirectory);
}

bool seqSetUp::processDefaultReader(bool readInNamesRequired){
	return processDefaultReader(readInFormatsAvailable_, readInNamesRequired);
}

bool seqSetUp::processDefaultReader(const VecStr & formats, bool readInNamesRequired) {
	bool passed = true;
	if (!processReadInNames(formats, readInNamesRequired)) {
		passed = false;
	}
	std::string outFormat = "";

	if(setOption(outFormat, "--outFormat", "Format of out sequence file", false, "Writing Output")){
		pars_.ioOptions_.outFormat_ = SeqIOOptions::getOutFormat(outFormat);
	}
	if (SeqIOOptions::outFormats::NOFORMAT != pars_.ioOptions_.outFormat_) {
		pars_.ioOptions_.out_.outExtention_ = SeqIOOptions::getOutExtension(
				pars_.ioOptions_.outFormat_);
	}

	setOption(pars_.ioOptions_.out_.outExtention_, "--outExtention", "Extension of out file",  false, "Writing Output");

	processWritingOptions();
	return passed;
}

void seqSetUp::processWritingOptions() {
	setOption(pars_.ioOptions_.out_.overWriteFile_, "--overWrite",
			"Over Write Existing Files",  false, "Writing Output");
	setOption(pars_.ioOptions_.out_.append_, "--appendFile", "Append to file",  false, "Writing Output");
	if (setOption(pars_.ioOptions_.out_.outFilename_, "--out", "Out Filename",  false, "Writing Output")) {
		if("" != pars_.ioOptions_.out_.outFilename_.extension().string()){
			pars_.ioOptions_.out_.outExtention_ =
					pars_.ioOptions_.out_.outFilename_.extension().string();
		}
	}
}

void seqSetUp::processWritingOptions(OutOptions & opts) {
	setOption(opts.overWriteFile_, "--overWrite", "Over Write Existing Files",  false, "Writing Output");
	setOption(opts.append_, "--appendFile", "Append to file",  false, "Writing Output");
	if (setOption(opts.outFilename_, "--out", "Out Filename",  false, "Writing Output")) {
		if("" != opts.outFilename_.extension().string()){
			opts.outExtention_ = opts.outFilename_.extension().string();
		}
	}
}

bool seqSetUp::processRefFilename(bool required) {
	setOption(pars_.refIoOptions_.processed_, "--refProcessed",
			"Reference Name Has Abundance Info");
	if (commands_.hasFlagCaseInsenNoDash("--refFastq")) {
		pars_.refIoOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQ;
		pars_.refIoOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
	} else if (commands_.hasFlagCaseInsenNoDash("--ref")) {
		pars_.refIoOptions_.inFormat_ = SeqIOOptions::inFormats::FASTA;
		pars_.refIoOptions_.outFormat_ = SeqIOOptions::outFormats::FASTA;
	}
	processGapRef();
	return setOption(pars_.refIoOptions_.firstName_,
			"--ref,--refFastq", "Reference Fasta or Fastq File Name", required);
}

bool seqSetUp::processSeq(bool required) {
	return processSeq(pars_.seq_, "--seq", "Sequence", required);
}

bool seqSetUp::processSeq(std::string& inputSeq, const std::string& flag,
		const std::string& parName, bool required) {
	pars_.seqObj_.seqBase_.seq_ = inputSeq;
	bool ret =  processSeq(pars_.seqObj_.seqBase_, flag, parName, required);
	inputSeq = pars_.seqObj_.seqBase_.seq_;
	return ret;
}

bool seqSetUp::processSeq(seqInfo& inputSeq, const std::string& flag,
		const std::string& parName, bool required) {
	bool passed = setOption(inputSeq.seq_, flag, parName, required);
	//std::cout <<"1 "<< inputSeq << std::endl;
	std::string originalSeq = inputSeq.seq_;
	if (bfs::path(inputSeq.seq_).filename().string().length() <= 255 && bfs::exists(inputSeq.seq_)) {
		std::string firstLine = bib::files::getFirstLine(originalSeq);
		inputSeq = seqInfo("seq", firstLine);
		//std::cout << "2 "<< inputSeq << std::endl;
		if (firstLine[0] == '>') {
			//std::cout <<"3 "<< inputSeq << std::endl;
			SeqIOOptions opts = SeqIOOptions::genFastaIn(originalSeq);
			SeqInput reader(opts);
			reader.openIn();
			seqInfo seq;
			reader.readNextRead(seq);
			inputSeq = seq;
		} else if (firstLine[0] == '@') {
			std::ifstream inFile(inputSeq.seq_);
			std::string nextLine = "";
			bib::files::crossPlatGetline(inFile, nextLine);	//name line again
			bib::files::crossPlatGetline(inFile, nextLine);	//seq line
			inputSeq = seqInfo(firstLine.substr(1), nextLine);
		}
	} else {
		inputSeq = seqInfo("seq", inputSeq.seq_);
	}
	//std::cout <<"4 "<< inputSeq << std::endl;
	return passed;
}

bool seqSetUp::processVerbose() {
  return setOption(pars_.verbose_, "-v,--verbose", "Verbose");
}

bool seqSetUp::processDebug() {
  return setOption(pars_.debug_, "--debug", "Debug");
}

bool seqSetUp::processQuiet() {
  return setOption(pars_.quiet_, "--quiet", "quiet");
}

void seqSetUp::processAlignerDefualts() {
  processGap();
  processQualThres();
  processScoringPars();
  processKmerProfilingOptions();
  processAlnInfoInput();
  bool alignScoreBased = false;
  setOption(alignScoreBased, "--scoreBased", "Scored Based Comparison For Alignments(defualt:event based)");
  pars_.colOpts_.alignOpts_.eventBased_ = !alignScoreBased;
}

void seqSetUp::processAlnInfoInput() {
	if (setOption(pars_.alnInfoDirName_, "--alnInfoDir", "Alignment Info Cache Directory Name", false, "Caching Alignments")) {
		if (!setOption(pars_.outAlnInfoDirName_, "--outAlnInfoDir", "Alignment Info Cache Out Directory Name", false, "Caching Alignments")) {
			pars_.outAlnInfoDirName_ = pars_.alnInfoDirName_;
		}
		pars_.writingOutAlnInfo_ = true;
	} else {
		if (setOption(pars_.outAlnInfoDirName_, "--outAlnInfoDir", "Alignment Info Cache Out Directory Name", false, "Caching Alignments")) {
			pars_.writingOutAlnInfo_ = true;
		}
	}
}



//void seqSetUp::printAdditionaInputUsage(std::ostream& out,
//                                        const std::string& lowerRemove) {
//  out << bib::bashCT::bold << "Read in processing options: "<< bib::bashCT::reset  << std::endl;
//  if (lowerRemove == "remove") {
//    out << "1) --lower [option]: can be remove or upper, remove to remove the "
//           "lower case, upper to convert lower to upper case, defaults to "
//           "removing lowercase" << std::endl;
//  } else if (lowerRemove == "upper") {
//    out << "1) --lower [option]: can be remove or upper, remove to remove the "
//           "lower case, upper to convert lower to upper case, defaults to "
//           "converting to upper case" << std::endl;
//  } else {
//    out << "1) --lower [option]: can be remove or upper, remove to remove the "
//           "lower case, upper to convert lower to upper case, defaults to "
//           "doing nothing" << std::endl;
//  }
//  out << "2) --processed : whether the reads being read have frequency info in "
//         "their name, in the form of NAME_t[totalReadNum] or "
//         "NAME_f[fractionAmount], defaults to no info" << std::endl;
//  out << "3) --removeGaps : whether to remove any gap information (-) in the "
//         "input sequence if there is any, defaults to leaving gaps"
//      << std::endl;
//  out << "4) --noSpaceInName,-removeWhiteSpace : whether to remove any any white "
//         "space from the read names or to keep them in, defaults to leaving "
//         "white space" << std::endl;
//}
//void seqSetUp::printFileWritingUsage(std::ostream& out, bool all) {
//  // std::stringstream tempOut;
//  out << bib::bashCT::bold <<"File output options: "<< bib::bashCT::reset << std::endl;
//  out << "1) --out [option]: Name of output file" << std::endl;
//  out << "2) --outFormat [option]: Output file format"
//         ", options are fasta, fastaQual, and fastq" << std::endl;
//  if (all) {
//    out << "3) --overWrite,-overWriteFile : overwrite the file if it "
//           " already exists, defaults to not overwriting" << std::endl;
//    out << "4) --exitOnFailureToWrite : if fail to write the file whether to "
//           "keep going or exit, defaults to keep going " << std::endl;
//  }
//  // out << cleanOut(tempOut.str(), width_, indent_);
//}
//void seqSetUp::printKmerProfilingUsage(std::ostream& out) {
//  // std::stringstream tempOut;
//  out << bib::bashCT::bold <<"kmer options"<< bib::bashCT::reset << std::endl;
//  out << "1) --kLength [option]: The length of the k mer check, defaults to " << pars_.colOpts_.kmerOpts_.kLength_ << std::endl;
//  out << "2) --kAnywhere : check any kmers found anywhere, defaults to "
//         "checking for kmers at the position found" << std::endl;
//  out << "3) --runCutOff [option]: kmer occurrence number cut off "
//         "to count as real sequence, defaults to 1, which means it has "
//         "to occur in more than one read to be considered real" << std::endl;
//  out << "\tif given with a percent sign, will make the cutoff the "
//         "percentage, eg. 0.5% has to occur in more than 0.5% of the "
//         "reads" << std::endl;
//  // out << cleanOut(tempOut.str(), width_, indent_);
//}
//void seqSetUp::printQualThresUsage(std::ostream& out) {
//  // std::stringstream tempOut;
//  out << bib::bashCT::bold <<"Quality Threshold Options"<< bib::bashCT::reset  << std::endl;
//  out << "1) --qualThres [option]: Quality threshold for high qual "
//         "mismatch, given in the format of 20,15, where 20 is the "
//         "primary quality (of mismatch) and 15 is the secondary quality "
//         "(of flanking bases)" << std::endl;
//  out << "2) --qualThresWindow, -qwindow [option]: Number of flanking qualities "
//         "to "
//         "included in quality neighborhood, a value of 5 would mean 5 "
//         "trailing qualities and 5"
//         "leading qualities would be in the window for a total of 11 qualities "
//         "(includig the mismatch"
//         " base, defaults to " << pars_.qScorePars_.qualThresWindow_ << std::endl;
//  // out << cleanOut(tempOut.str(), width_, indent_);
//}
//
//void seqSetUp::printGapUsage(std::ostream & out) const {
//	out << bib::bashCT::bold << "Gap Scoring options" << bib::bashCT::reset
//			<< "\n";
//	out << "--gap [option]: Gap penalty, given in the format 7,0.5 "
//			"where 7 is the gap open penalty and 0.5 is the gap extension"
//			<< "\n";
//	out << "--gapLeft [option]: Gap penalty for putting gaps at the beginning"
//			" of the sequence, same format as -gap" << "\n";
//	out << "--gapRight [option]: Gap penalty for putting gaps at the end"
//			" of the sequence, same format as -gap" << "\n";
//	out << "--gapAll [option]: Gap penalty at all locations, would be like calling"
//			" all three options,-gap, -gapRight, -gapLeft with the same parameters, same format as -gap\n";
//}
//
//
//void seqSetUp::printAlignmentUsage(std::ostream& out) {
//  // std::stringstream tempOut;
//  out << bib::bashCT::bold <<"Alignment options"<< bib::bashCT::reset  << std::endl;
//  printGapUsage(out);
//  out << "4a) --scoreMatrix : A filename for an alignment scoring matrix "
//         "of a custom scoring matrix" << std::endl;
//  out << "4b) --generalMatch [option]: If no score matrix is given the score "
//         " of any match will be this, defaults to 2" << std::endl;
//  out << "4c) --generalMismatch [option]: If no score matrix is given the score "
//         " of any mismatch will be this,"
//         " defaults to -2" << std::endl;
//  out << "5) --local : do local alignment instead of global, defaults to global"
//      << std::endl;
//  out << "6) --countEndGaps: Whether or not to count end gaps in the "
//         "alignment comparison" << std::endl;
//  out << "7) --noHomopolymerWeighting,-noHWeighting: In alignment comparison, "
//         "do not count"
//         " indels in homopolymer differnt from other indels" << std::endl;
//  // out << cleanOut(tempOut.str(), width_, indent_);
//}
//void seqSetUp::printReferenceComparisonUsage(std::ostream& out) {
//  // std::stringstream tempOut;
//  out << bib::bashCT::bold <<"Reference comparison options:"<< bib::bashCT::reset  << std::endl;
//  out << "1a) --ref [option]: Name of a reference file in fasta format "
//         "for"
//         " references sequences to be read from " << std::endl;
//  out << "1b) --refFastq [option]: Same as above but if the file "
//         "is in"
//         " fastq format " << std::endl;
//  out << "2) --refProcessed : If the references have "
//         "frequency info"
//         " in their name, in the form of _t[readNumber] or _f[Fraction]"
//      << std::endl;
//  // out << cleanOut(tempOut.str(), width_, indent_);
//}
//void seqSetUp::printAlnInfoDirUsage(std::ostream& out) {
//  // std::stringstream tempOut;
//  out << bib::bashCT::bold <<"Save Alignments Options: "<< bib::bashCT::reset  << std::endl;
//  out << "1) --alnInfoDir [option] : Name of a directory to save alingments"
//         " if the directory already exists it assumes it has saved alingments "
//         "and these will be read in and used, if the directory does not exist "
//         "a new directory will be created with this name and alingment will be "
//         "save here" << std::endl;
//  out << "2) --outAlnInfoDir,-outAln,-outAlnInfo [option] : Name of a new "
//         "directory to save alignments in, if -alnInfoDir given but nothing "
//         "for this option, the out directory will default to the -alnInfodir "
//         "name" << std::endl;
//  ;
//  // out << cleanOut(tempOut.str(), width_, indent_);
//}
//void seqSetUp::printAdditionalClusteringUsage(std::ostream& out) {
//  // std::stringstream tempOut;
//  out << bib::bashCT::bold <<"Optional Clustering options"<< bib::bashCT::reset << std::endl;
//  out << "1) --bestMatch : Has the program cluster the reads looking "
//         "for the best match, defaults to the first read that meets "
//         "the given parameters" << std::endl;
//  out << "2) --bestMatchCheck [option] : The number of reads to look "
//         "for when -bestMatch has been switched on, defaults to 10"
//      << std::endl;
//  out << "3) --largestFirst : Cluster the reads by comparing the "
//         "largest clusters to each other first, defaults to taking the "
//         "smallest cluster and comparing to the largest and making "
//         "it's way up" << std::endl;
//  // out << cleanOut(tempOut.str(), width_, indent_);
//}

}  // namespace bib
