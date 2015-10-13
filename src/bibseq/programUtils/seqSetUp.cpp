#include "seqSetUp.hpp"
#include "bibseq/programUtils/runningParameters.hpp"
#include <bibcpp/bashUtils.h>
#include <bibcpp/files/fileUtilities.hpp>
namespace bibseq {


seqSetUp::seqSetUp(int argc, char* argv[]) : bib::progutils::programSetUp(argc, argv) {
  initializeDefaults();
}

seqSetUp::seqSetUp(const bib::progutils::commandLineArguments& inputCommands)
    : bib::progutils::programSetUp(inputCommands) {
  initializeDefaults();
}

seqSetUp::seqSetUp(const MapStrStr& inputCommands) : bib::progutils::programSetUp(inputCommands) {
  initializeDefaults();
}

void seqSetUp::initializeDefaults() {
  ioOptions_.firstName_ = "";
  ioOptions_.secondName_ = "";
  ioOptions_.inFormat_ = "";
  ioOptions_.outFormat_ = "";
  ioOptions_.outFilename_ = "out";
  ioOptions_.processed_ = false;
  ioOptions_.append_ = false;
  //
  ioOptions_.removeGaps_ = false;
  ioOptions_.lowerCaseBases_ = "nothing";
  //
  ioOptions_.overWriteFile_ = false;
  ioOptions_.exitOnFailureToWrite_ = true;
  //
  ioOptions_.includeWhiteSpaceInName_ = true;
  //
  seq_ = "";
  seqObj_ = readObject(seqInfo("", seq_));

  //
  directoryName_ = "";
  overWriteDir_ = false;
  //
  refFilename_ = "";
  refSecondName_ = "";
  refFormat_ = "";
  refProcessed_ = false;

  //
  verbose_ = false;
  debug_ = false;
  quiet_ = false;

  //
  gapRef_ = "5,1";
  gapInfoRef_.processGapStr(gapRef_, gapInfoRef_.gapOpen_, gapInfoRef_.gapExtend_);

  gapLeftRef_ = "0,0";
  gapInfoRef_.processGapStr(gapLeftRef_, gapInfoRef_.gapLeftOpen_, gapInfoRef_.gapLeftExtend_);

  gapRightRef_ = "0,0";
  gapInfoRef_.processGapStr(gapRightRef_, gapInfoRef_.gapRightOpen_, gapInfoRef_.gapRightExtend_);
  gapInfoRef_.setIdentifer();
  //
  gap_ = "7,1";
  gapInfo_.processGapStr(gap_, gapInfo_.gapOpen_, gapInfo_.gapExtend_);

  gapLeft_ = "7,1";
  gapInfo_.processGapStr(gapLeft_, gapInfo_.gapLeftOpen_, gapInfo_.gapLeftExtend_);

  gapRight_ = "0,0";
  gapInfo_.processGapStr(gapRight_, gapInfo_.gapRightOpen_, gapInfo_.gapRightExtend_);
  gapInfo_.setIdentifer();
  //
  local_ = false;
  generalMatch_ = 2;
  generalMismatch_ = -2;
  scoring_ = substituteMatrix::createDegenScoreMatrix(1,-1);
  //
  countEndGaps_ = false;
  weightHomopolymers_ = true;

  qualThres_ = "20,15";
  primaryQual_ = 20;
  secondaryQual_ = 15;
  qualThresWindow_ = 2;

  eventBased_ = false;
  //
  alnInfoDirName_ = "";
  outAlnInfoDirName_ = "";
  writingOutAlnInfo_ = false;

  //
  runCutOffString_ = "1";
  runCutoff_ = 1;

  kLength_ = 9;
  kmersByPosition_ = true;
  expandKmerPos_ = false;
  expandKmerSize_ = 5;

  //general clustering
  skipOnLetterCounterDifference_ = false;
  fractionDifferenceCutOff_ = 0.05;
  adjustHomopolyerRuns_ = false;
  largestFirst_ = false;
  firstMatch_ = false;
  bestMatchCheck_ = 10;
}

void seqSetUp::processClusteringOptions(){
  processSkipOnNucComp();
  processAdjustHRuns();
  setOption(largestFirst_,   "-largestFirst",   "Compare largest clusters first");
  setOption(firstMatch_,     "-firstMatch",     "Settle for first Match in Clustering");
  setOption(bestMatchCheck_, "-bestMatchCheck", "Best Match Check Number");
}


void seqSetUp::processGap() {
  // check command line for gap settings
  if (setOption(gap_, "-gapAll", "Gap Penalties for All (middle and end gaps)")) {
    gapInfo_ = gapScoringParameters (gap_);
    gapLeft_ = gap_;
    gapRight_ = gap_;
  }else{
  	setOption(gap_, "-gap", "Gap Penalties for Middle Gap");
  	setOption(gapLeft_, "-gapLeft", "Gap Penalties for Left End Gap");
  	setOption(gapRight_, "-gapRight", "Gap Penalties for Right End Gap");
  	// get the gap penalty
  	gapInfo_.processGapStr(gap_, gapInfo_.gapOpen_, gapInfo_.gapExtend_);
  	gapInfo_.processGapStr(gapLeft_, gapInfo_.gapLeftOpen_, gapInfo_.gapLeftExtend_);
  	gapInfo_.processGapStr(gapRight_, gapInfo_.gapRightOpen_, gapInfo_.gapRightExtend_);
  }
  gapInfo_.setIdentifer();
}

void seqSetUp::processGapRef() {
  // check command line for gap settings
  if (setOption(gapRef_, "-gapRefAll", "Gap Penalties for Ref All")) {
    gapInfoRef_ = gapScoringParameters (gapRef_);
    gapLeftRef_ = gapRef_;
    gapRightRef_ = gapRef_;
  }else{
  	setOption(gapRef_, "-gap", "Gap Penalties for Ref Middle Gap");
  	setOption(gapLeftRef_, "-gapLeft", "Gap Penalties for Ref Left End Gap");
  	setOption(gapRightRef_, "-gapRight", "Gap Penalties for Ref Right End Gap");
  	// get the gap penalty
  	gapInfoRef_.processGapStr(gapRef_, gapInfoRef_.gapOpen_, gapInfoRef_.gapExtend_);
  	gapInfoRef_.processGapStr(gapLeftRef_, gapInfoRef_.gapLeftOpen_, gapInfoRef_.gapLeftExtend_);
  	gapInfoRef_.processGapStr(gapRightRef_, gapInfoRef_.gapRightOpen_, gapInfoRef_.gapRightExtend_);
  }
  gapInfoRef_.setIdentifer();
}

void seqSetUp::processQualThres() {
  // check command line for qualThres setting
  setOption(qualThres_, "-qualThres", "Quality Thresholds, should go PrimaryQual,SecondaryQual");
  // get the qualities
  auto qualToks = tokenizeString(qualThres_, ",");
  if(qualToks.size() != 2){
  	throw std::runtime_error{"QaulThres should be two numbers separated by a comma, eg 20,15, not " + qualThres_};
  }
  primaryQual_ = std::stoi(qualToks[0]);
  secondaryQual_ = std::stoi(qualToks[1]);

  setOption(qualThresWindow_, "-qualThresWindow",
            "Quality Threshold Window Length");
}

void seqSetUp::processIteratorMap(
    std::string& parameters,
    std::map<int, std::vector<double>>& iteratorMap) {

	if(!fexists(parameters)){
		failed_ = true;
		warnings_.emplace_back(bib::bashCT::red + bib::bashCT::bold
				+ "File " + parameters + " doesn't exist"
				+ bib::bashCT::reset);
	}
	iteratorMap = runningParameters::processParameters(parameters);
}

void seqSetUp::processIteratorMapOnPerId(
    std::string& parameters,
    std::map<int, std::vector<double>>& iteratorMap) {

	if(!fexists(parameters)){
		failed_ = true;
		warnings_.emplace_back(bib::bashCT::red + bib::bashCT::bold
				+ "File " + parameters + " doesn't exist"
				+ bib::bashCT::reset);
	}
	iteratorMap = runningParameters::processParametersPerId(parameters);
}

void seqSetUp::processKmerLenOptions(){
  setOption(kLength_, "-kLength", "Kmer Length");
}

void seqSetUp::processKmerProfilingOptions() {
  setOption(runCutOffString_, "-runCutOff", "Kmer_frequencey_cutoff");
  processKmerLenOptions();
  bool forKmerProfiling = true;
  if (forKmerProfiling && kLength_ % 2 == 0) {
    kLength_--;
    warnings_.emplace_back("-kLength needs to be odd, not even, changing to " +
                           estd::to_string(kLength_));
  }
  bool kAnywhere = false;
  setOption(kAnywhere, "-kAnywhere", "Count Kmers without regard for position");
  kmersByPosition_ = !kAnywhere;
  setOption(expandKmerPos_, "-expandKmerPos", "Expand Kmer Position Found At");
  setOption(expandKmerSize_, "-expandKmerSize", "Expand Kmer Size for extending where kmers where found");
}

void seqSetUp::processScoringPars() {
  setOption(local_, "-local", "Local_alignment");
  setOption(countEndGaps_, "-countEndGaps", "CountEndGaps");
  bool noHomopolymerWeighting = false;
  setOption(noHomopolymerWeighting,
                     "-noHomopolymerWeighting",
                     "Don't do Homopolymer Weighting");
  weightHomopolymers_ = !noHomopolymerWeighting;
  std::string scoreMatrixFilename = "";
  bool degenScoring = false;
  bool caseInsensitive = false;
  if(setOption(scoreMatrixFilename, "-scoreMatrix", "Score Matrix Filename")){
  	scoring_ = substituteMatrix(commands_["-scorematrix"]);
  } else {
    setOption(generalMatch_, "-generalMatch,-match", "generalMatch");
    setOption(generalMismatch_, "-generalMismatch,-mismatch", "generalMismatch");
    setOption(degenScoring, "-degen", "Use Degenerative Base Scoring");
    setOption(caseInsensitive, "-caseInsensitive", "Use Case Insensititive Scoring");
    if(degenScoring){
    	if(caseInsensitive){
    		scoring_ =
    		    		substituteMatrix::createDegenScoreMatrixCaseInsensitive(generalMatch_, generalMismatch_);
    	}else{
    		scoring_.setWithCaseInsensitive(generalMatch_, generalMismatch_);
    	}
    }else {
    	if(caseInsensitive){
    		scoring_ =
    		    		substituteMatrix::createDegenScoreMatrixCaseInsensitive(generalMatch_, generalMismatch_);
    	}else{
    		scoring_ =
    		    		substituteMatrix(generalMatch_, generalMismatch_);
    	}
    }
  }
}


void seqSetUp::processSkipOnNucComp(){
  setOption(skipOnLetterCounterDifference_, "-skip",
                  "skipOnLetterCounterDifference");
  setOption(fractionDifferenceCutOff_, "-skipCutOff",
                  "fractionDifferenceCutOff");
}

void seqSetUp::processAdjustHRuns(){
  setOption(adjustHomopolyerRuns_,
            "-adjustHomopolyerRuns",
            "Adjust Homopolyer Runs To Be Same Qual");
}

bool seqSetUp::processReadInNames(bool required) {
  // still need if multiple files given, like -fasta and -fastq at the same time
	VecStr readInFormats {"-sff", "-sffBin", "-fasta", "-fastq", "-fastqgz", "-bam", "-fastqgz"};
	VecStr foundInFormats;
	if (gettingFlags_ || printingHelp_) {
		bib::progutils::flag sffFlagOptions(ioOptions_.firstName_, "-sff",
				"Input sequence filename, sff text file, if required only one format is accepted",
				required);
		bib::progutils::flag sffBinFlagOptions(ioOptions_.firstName_, "-sffBin",
				"Input sequence filename, sff binary file, if required only one format is accepted",
				required);
		bib::progutils::flag fastaFlagOptions(ioOptions_.firstName_, "-fasta",
				"Input sequence filename, fasta text file, if required only one format is accepted",
				required);
		bib::progutils::flag fastqFlagOptions(ioOptions_.firstName_, "-fastq",
				"Input sequence filename, fastq text file, if required only one format is accepted",
				required);
		bib::progutils::flag fastqgzFlagOptions(ioOptions_.firstName_, "-fastqgz",
				"Input sequence filename, fastq gzipped file, if required only one format is accepted",
				required);
		bib::progutils::flag bamFlagOptions(ioOptions_.firstName_, "-bam",
				"Input sequence filename, bam file, if required only one format is accepted",
				required);
		flags_.addFlag(sffFlagOptions);
		flags_.addFlag(sffBinFlagOptions);
		flags_.addFlag(fastaFlagOptions);
		flags_.addFlag(fastqFlagOptions);
		flags_.addFlag(bamFlagOptions);
		flags_.addFlag(fastqgzFlagOptions);
	}
	//process out format
	if(commands_.containsFlagCaseInsensitive("-fasta")){
		ioOptions_.inFormat_ = "fasta";
		ioOptions_.outFormat_ = "fasta";
		ioOptions_.outExtention_ = ".fasta";
		foundInFormats.emplace_back("-fasta");
	}
	if(commands_.containsFlagCaseInsensitive("-sff")){
		ioOptions_.inFormat_ = "sff";
		ioOptions_.outFormat_ = "fastq";
		ioOptions_.outExtention_ = ".fastq";
		foundInFormats.emplace_back("-sff");
	}
	if(commands_.containsFlagCaseInsensitive("-sffBin")){
		ioOptions_.inFormat_ = "sffbin";
		ioOptions_.outFormat_ = "fastq";
		ioOptions_.outExtention_ = ".fastq";
		foundInFormats.emplace_back("-sffBin");
	}
	if(commands_.containsFlagCaseInsensitive("-bam")){
		ioOptions_.inFormat_ = "bam";
		ioOptions_.outFormat_ = "fastq";
		ioOptions_.outExtention_ = ".fastq";
		foundInFormats.emplace_back("-bam");
	}
	if(commands_.containsFlagCaseInsensitive("-fastq")){
		ioOptions_.inFormat_ = "fastq";
		ioOptions_.outFormat_ = "fastq";
		ioOptions_.outExtention_ = ".fastq";
		foundInFormats.emplace_back("-fastq");
	}
	if(commands_.containsFlagCaseInsensitive("-fastqgz")){
		ioOptions_.inFormat_ = "fastqgz";
		ioOptions_.outFormat_ = "fastq";
		ioOptions_.outExtention_ = ".fastq";
		foundInFormats.emplace_back("-fastqgz");
	}
	if(foundInFormats.size() > 1){
    std::stringstream tempOut;
    tempOut << bib::bashCT::bold
    		<< "Found multiple read in options, should only have one"
            << std::endl;
    tempOut << vectorToString(foundInFormats, ",") + bib::bashCT::reset;
    tempOut << std::endl;
    addOtherVec(warnings_, streamToVecStr(tempOut));
    failed_ = true;
    return false;
	} else if (foundInFormats.empty()){
		if(required){
	    std::stringstream tempOut;
	    tempOut << bib::bashCT::bold
	    				<< "Did not find a recognizable read in option"
	            << std::endl;
	    tempOut << "options include: "
	            << bib::bashCT::red << "-fasta, -fasta/-qual, -stub, -fastq, -fastqgz"
	            << bib::bashCT::reset << std::endl;
	    tempOut << "Command line arguments" << std::endl;
	    writeOutCommandLineArguments(commands_.arguments_, tempOut);
	    addOtherVec(warnings_, streamToVecStr(tempOut));
	    failed_ = true;
		}
    return false;
	} else {
		setOption(ioOptions_.firstName_, foundInFormats.front(), "In Sequence Filename");
		if(foundInFormats.front() == "-fasta"){
			if(setOption(ioOptions_.secondName_, "-qual", "Name of the quality file")){
				ioOptions_.inFormat_ = "fastaQual";
				ioOptions_.outFormat_ = "fastq";
			}
		}
		setOption(ioOptions_.outFormat_, "-outFormat", "Format of out sequence file");
		setOption(ioOptions_.outExtention_, "-outExtention", "Extention of out file");
		setOption(ioOptions_.outFilename_, "-out", "Name of the out sequence file");
		return true;
	}
}

void seqSetUp::processDirectoryOutputName(const std::string& defaultName,
                                          bool mustMakeDirectory) {
	setOption(overWriteDir_, "--overWriteDir", "If the directory already exists over write it");
  directoryName_ = "./";
  if (setOption(directoryName_, "-dout", "Output Directory Name")) {
    if (!failed_) {
      std::string newDirectoryName = "./" +
      		replaceString(replaceString(directoryName_, "./", ""), "TODAY", getCurrentDate()) +"/";
    	if(bib::files::bfs::exists(newDirectoryName) && overWriteDir_){
    		bib::files::rmDirForce(newDirectoryName);
    	}
      directoryName_ =
          bib::files::makeDir("./", replaceString(directoryName_, "./", ""));
    }
  } else {
    if (mustMakeDirectory && !failed_) {
      std::string newDirectoryName = "./" +
      		replaceString(defaultName, "TODAY", getCurrentDate()) +"/";
    	if(bib::files::bfs::exists(newDirectoryName) && overWriteDir_){
    		bib::files::rmDirForce(newDirectoryName);
    	}
      directoryName_ = bib::files::makeDir("./", defaultName);
    }
  }
}

void seqSetUp::processDirectoryOutputName(bool mustMakeDirectory) {
  std::string seqName = bib::files::getFileName(ioOptions_.firstName_) + "_" +
                        commands_["-program"] + "_" + getCurrentDate();
  processDirectoryOutputName(seqName, mustMakeDirectory);
}

bool seqSetUp::processDefaultReader(bool readInNamesRequired) {
  bool passed = true;
  if (!processReadInNames(readInNamesRequired)) {
    passed = false;
  }
  //setOption(ioOptions_.forceWrite_, "-overWrite", "Over Write Fill");
  setOption(ioOptions_.processed_, "-processed", "Processed, Input Sequence Name contains abundance info");
  setOption(ioOptions_.outFilename_, "-out", "Out Filename");
  setOption(ioOptions_.lowerCaseBases_, "-lower", "How to handle Lower Case Bases");
  setOption(ioOptions_.removeGaps_, "-removeGaps", "Remove Gaps from Input Sequences");
  setBoolOptionFalse(ioOptions_.includeWhiteSpaceInName_,
                     "-noWSpaceInName",
                     "Remove White Space From Sequence Names");
  processWritingOptions();
  return passed;
}

void seqSetUp::processWritingOptions() {
  setOption(ioOptions_.overWriteFile_, "-overWrite",
            "Over Write Existing Files");
  //setOption(ioOptions_.exitOnFailureToWrite_, "-exitOnFailureToWrite","Exit On Failure To Write");
  setOption(ioOptions_.append_, "-appendFile",
              "Append to file");
}


bool seqSetUp::processRefFilename(bool required) {
  setOption(refProcessed_, "-refProcessed",
            "Reference Name Has Abundance Info");
  if (commands_.containsFlagCaseInsensitive("-refFastq") ||
      commands_.containsFlagCaseInsensitive("-expectFastq")) {
    refFormat_ = "fastq";
  } else if (commands_.containsFlagCaseInsensitive("-ref") ||
             commands_.containsFlagCaseInsensitive("-expect")) {
    refFormat_ = "fasta";
  }
  processGapRef();
  return setOption(refFilename_, "-ref,-refFastq,-expect,-expectFastq",
                   "ReferenceFileName", required);
}

bool seqSetUp::processSeq(bool required) {
  return processSeq(seq_, "-seq", "Sequence", required);
}

bool seqSetUp::processSeq(std::string& inputSeq, const std::string& flag,
                          const std::string& parName, bool required) {
  bool passed = setOption(inputSeq, flag, parName, required);
  //std::cout <<"1 "<< inputSeq << std::endl;
  if (fexists(inputSeq)) {
  	std::ifstream inFile(inputSeq);

    inputSeq = bib::files::getFirstLine(inputSeq);
    seqObj_ = readObject(seqInfo("seq", inputSeq));
    //std::cout << "2 "<< inputSeq << std::endl;
    if(inputSeq[0] == '>'){
    	//std::cout <<"3 "<< inputSeq << std::endl;
    	readObjectIO reader;
    	reader.readFastaStream(inFile, false, false);
    	inputSeq = reader.reads.front().seqBase_.seq_;
    	//NEEDS TO BE FIXED TO BE GIVEN BY REF AND USE THE DEFAULT
    	seqObj_ = reader.reads.front();
    }else if(inputSeq[0] == '@'){
    	//std::cout << "4" << std::endl;
    	std::string nextLine = "";
    	getline(inFile, nextLine);
    	getline(inFile, nextLine);
      seqObj_ = readObject(seqInfo(inputSeq.substr(1), nextLine));
      inputSeq = seqObj_.seqBase_.seq_;
    }
  } else{
  	seqObj_ = readObject(seqInfo("seq", inputSeq));
  }
  //std::cout <<"4 "<< inputSeq << std::endl;
  return passed;
}
bool seqSetUp::processVerbose() {
  return setOption(verbose_, "-v,--verbose", "Verbose");
}

bool seqSetUp::processDebug() {
  return setOption(debug_, "--debug", "Debug");
}

bool seqSetUp::processQuiet() {
  return setOption(quiet_, "--quiet", "quiet");
}

void seqSetUp::processAlignerDefualts() {
  processGap();
  processQualThres();
  processKmerProfilingOptions();
  processScoringPars();
  processAlnInfoInput();
  setOption(eventBased_, "-eventBased", "Event Based Comparison");
}

void seqSetUp::processAlnInfoInput() {
	if (setOption(alnInfoDirName_, "-alnInfoDir", "alnInfoDirName")) {
		if (!setOption(outAlnInfoDirName_, "-outAlnInfoDir", "alnInfoDirName")) {
			outAlnInfoDirName_ = alnInfoDirName_;
		}
		writingOutAlnInfo_ = true;
	} else {
		if (setOption(outAlnInfoDirName_, "-outAlnInfoDir", "alnInfoDirName")) {
			writingOutAlnInfo_ = true;
		}
	}
}

void seqSetUp::printInputUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold << "Input options:" << bib::bashCT::reset << std::endl;
  out << "1a) -stub [option]: Stub name for the .fasta and .fasta.qual files"
      << std::endl;
  out << "1b) -fasta [option]: Full name of the fasta file" << std::endl;
  out << "1c) -qual [option]: Full name of the quality file" << std::endl;
  out << "1d) -fastq [option]: Full name of the fastq file" << std::endl;
  // out << cleanOut(tempOut.str(), width_, indent_);
}
void seqSetUp::printAdditionaInputUsage(std::ostream& out,
                                        const std::string& lowerRemove) {
  out << bib::bashCT::bold << "Read in processing options: "<< bib::bashCT::reset  << std::endl;
  if (lowerRemove == "remove") {
    out << "1) -lower [option]: can be remove or upper, remove to remove the "
           "lower case, upper to convert lower to upper case, defaults to "
           "removing lowercase" << std::endl;
  } else if (lowerRemove == "upper") {
    out << "1) -lower [option]: can be remove or upper, remove to remove the "
           "lower case, upper to convert lower to upper case, defaults to "
           "converting to upper case" << std::endl;
  } else {
    out << "1) -lower [option]: can be remove or upper, remove to remove the "
           "lower case, upper to convert lower to upper case, defaults to "
           "doing nothing" << std::endl;
  }
  out << "2) -processed : whether the reads being read have frequency info in "
         "their name, in the form of NAME_t[totalReadNum] or "
         "NAME_f[fractionAmount], defaults to no info" << std::endl;
  out << "3) -removeGaps : whether to remove any gap infomation (-) in the "
         "input sequence if there is any, defaults to leaving gaps"
      << std::endl;
  out << "4) -noSpaceInName,-removeWhiteSpace : wheter to remove any any white "
         "space from the read names or to keep them in, defaults to leaving "
         "white space" << std::endl;
}
void seqSetUp::printFileWritingUsage(std::ostream& out, bool all) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"File output options: "<< bib::bashCT::reset << std::endl;
  out << "1) -out [option]: Name of output file" << std::endl;
  out << "2) -outFormat [option]: Output file format"
         ", options are fasta, fastaQual, and fastq" << std::endl;
  if (all) {
    out << "3) -overWrite,-overWriteFile : overwrite the file if it "
           " already exists, defaults to not overwriting" << std::endl;
    out << "4) -exitOnFailureToWrite : if fail to write the file whether to "
           "keep going or exit, defaults to keep going " << std::endl;
  }
  // out << cleanOut(tempOut.str(), width_, indent_);
}
void seqSetUp::printKmerProfilingUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"kmer options"<< bib::bashCT::reset << std::endl;
  out << "1) -kLength [option]: The length of the k mer check, defaults to " << kLength_ << std::endl;
  out << "2) -kAnywhere : check any kmers found anywhere, defaults to "
         "checking for kmers at the position found" << std::endl;
  out << "3) -runCutOff [option]: kmer occurrence number cut off "
         "to count as real sequence, defaults to 1, which means it has "
         "to occur in more than one read to be considered real" << std::endl;
  out << "\tif given with a percent sign, will make the cutoff the "
         "percentage, eg. 0.5% has to occur in more than 0.5% of the "
         "reads" << std::endl;
  out << "4) -qualRunCutOff [option]: kmer occurrence number cut off "
      << "to raise the quality threshold for mismatches, same formating as "
         "-runCutOff" << std::endl;
  // out << cleanOut(tempOut.str(), width_, indent_);
}
void seqSetUp::printQualThresUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"Quality Threshold Options"<< bib::bashCT::reset  << std::endl;
  out << "1) -qualThres [option]: Quality threshold for high qual "
         "mismatch, given in the format of 20,15, where 20 is the "
         "primary quality (of mismatch) and 15 is the secondary quality "
         "(of flanking bases)" << std::endl;
  out << "2) -qualThresLowKmer [option]: Same as above but this threshold is "
         "used"
         " instead if low kmer checking is on, defaults to 30,25" << std::endl;
  out << "3) -qualThresWindow, -qwindow [option]: Number of flanking qualities "
         "to "
         "included in quality neighborhood, a value of 5 would mean five "
         "trailing qualities and five"
         "leading qualities would be in the window for a total of 11 qualities "
         "(includig the mismatch"
         " base, defaults to 5" << std::endl;
  // out << cleanOut(tempOut.str(), width_, indent_);
}

void seqSetUp::printGapUsage(std::ostream & out) const {
	out << bib::bashCT::bold << "Gap Scoring options" << bib::bashCT::reset
			<< "\n";
	out << "-gap [option]: Gap penalty, given in the format 7,0.5 "
			"where 7 is the gap open penalty and 0.5 is the gap extension"
			<< "\n";
	out << "-gapLeft [option]: Gap penalty for putting gaps at the beginning"
			" of the sequence, same format as -gap" << "\n";
	out << "-gapRight [option]: Gap penalty for putting gaps at the end"
			" of the sequence, same format as -gap" << "\n";
	out << "-gapAll [option]: Gap penalty at all locations, would be like calling"
			" all three options,-gap, -gapRight, -gapLeft with the same parameters, same format as -gap\n";
}


void seqSetUp::printAlignmentUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"Alignment options"<< bib::bashCT::reset  << std::endl;
  printGapUsage(out);
  out << "4a) -scoreMatrix : A filename for an alignment scoring matrix "
         "of a custom scoring matrix" << std::endl;
  out << "4b) -generalMatch [option]: If no score matrix is given the score "
         " of any match will be this, defaults to 2" << std::endl;
  out << "4c) -generalMismatch [option]: If no score matrix is given the score "
         " of any mismatch will be this,"
         " defaults to -2" << std::endl;
  out << "5) -local : do local alignment instead of global, defaults to global"
      << std::endl;
  out << "6) -countEndGaps: Whether or not to count end gaps in the "
         "alignment comparison" << std::endl;
  out << "7) -noHomopolymerWeighting,-noHWeighting: In alignment comparison, "
         "do not count"
         " indels in homopolymer differnt from other indels" << std::endl;
  // out << cleanOut(tempOut.str(), width_, indent_);
}
void seqSetUp::printReferenceComparisonUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"Reference comparison options:"<< bib::bashCT::reset  << std::endl;
  out << "1a) -ref,-expect [option]: Name of a reference file in fasta format "
         "for"
         " references sequences to be read from " << std::endl;
  out << "1b) -refFastq,-expectFastq [option]: Same as above but if the file "
         "is in"
         " fastq format " << std::endl;
  out << "2) -refprocessed,-expectedProcessed : If the references have "
         "frequency info"
         " in their name, in the form of _t[readNumber] or _f[Fraction]"
      << std::endl;
  // out << cleanOut(tempOut.str(), width_, indent_);
}
void seqSetUp::printAlnInfoDirUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"Save Alignments Options: "<< bib::bashCT::reset  << std::endl;
  out << "1) -alnInfoDir [option] : Name of a directory to save alingments"
         " if the directory already exists it assumes it has saved alingments "
         "and these will be read in and used, if the directory does not exist "
         "a new directory will be created with this name and alingment will be "
         "save here" << std::endl;
  out << "2) -outAlnInfoDir,-outAln,-outAlnInfo [option] : Name of a new "
         "directory to save alignments in, if -alnInfoDir given but nothing "
         "for this option, the out directory will default to the -alnInfodir "
         "name" << std::endl;
  ;
  // out << cleanOut(tempOut.str(), width_, indent_);
}
void seqSetUp::printAdditionalClusteringUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"Optional Clustering options"<< bib::bashCT::reset << std::endl;
  out << "1) -bestMatch : Has the program cluster the reads looking "
         "for the best match, defaults to the first read that meets "
         "the given parameters" << std::endl;
  out << "2) -bestMatchCheck [option] : The number of reads to look "
         "for when -bestMatch has been switched on, defaults to 10"
      << std::endl;
  out << "3) -largestFirst : Cluster the reads by comparing the "
         "largest clusters to each other first, defaults to taking the "
         "smallest cluster and comparing to the largest and making "
         "it's way up" << std::endl;
  // out << cleanOut(tempOut.str(), width_, indent_);
}

}  // namespace bib
