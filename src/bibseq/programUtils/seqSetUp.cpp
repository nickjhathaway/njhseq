//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "seqSetUp.hpp"
#include <bibcpp/bashUtils.h>
#include <bibcpp/files/fileUtils.hpp>
namespace bibseq {

void seqSetUp::processGap() {
  // check command line for gap settings
  if (setOption(gap_, "-gapAll", "Gap_penalitiesAll")) {
    gapInfo_ = gapScoringParameters (gap_);
    gapLeft_ = gap_;
    gapRight_ = gap_;
  }else{
  	setOption(gap_, "-gap", "Gap_penalities");
  	setOption(gapLeft_, "-gapLeft", "Gap_penalities");
  	setOption(gapRight_, "-gapRight", "Gap_penalities");
  	// get the gap penalty
  	gapInfo_.processGapStr(gap_, gapInfo_.gapOpen_, gapInfo_.gapExtend_);
  	gapInfo_.processGapStr(gapLeft_, gapInfo_.gapLeftOpen_, gapInfo_.gapLeftExtend_);
  	gapInfo_.processGapStr(gapRight_, gapInfo_.gapRightOpen_, gapInfo_.gapRightExtend_);
  }
  gapInfo_.setIdentifer();
}

void seqSetUp::processQualThres() {
  // check command line for qualThres setting
  setOption(qualThres_, "-qualThres", "Quality Thresholds, should go PrimaryQual,SecondaryQual");
  setOption(qualThesLowKmer_, "-qualThresLowKmer", "QualThresHold_LowKmer");
  // get the qualities
  auto qualToks = tokenizeString(qualThres_, ",");
  primaryQual_ = std::stoi(qualToks[0]);
  secondaryQual_ = std::stoi(qualToks[1]);

  auto qualToksLowKmer = tokenizeString(qualThesLowKmer_, ",");
  primaryQualLowKmer_ = std::stoi(qualToks[0]);
  secondaryQualLowKmer_ = std::stoi(qualToks[1]);

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
  std::string delim = "";
  delim.push_back(':');
  textFileReader reader(delim);
  reader.readFile(parameters);
  std::map<int, std::vector<double>>::iterator iteratorMapIter;
  int iters = 1;
  for (const auto& row : reader.fileContent.content_) {
    if (row[0][0] == 's') {
      continue;
    }
    std::vector<double> tempVect;
    for (const auto & colPos : iter::range(row.size() )) {
      if (colPos == 0) {
        if (row[colPos].find("%") != std::string::npos) {
          double percent = std::stod(
          		row[colPos].substr(0, row[colPos].size() - 1));
          if (percent >= 100) {
            tempVect.push_back(0.99);
          } else if (percent > 0) {
            tempVect.push_back(percent / 100.00);
          } else {
            tempVect.push_back(0.10);
          }
        } else if (stringToLowerReturn(row[colPos]) == "all") {
          tempVect.push_back(-1.00);
        } else {
          tempVect.push_back(std::stod(row[colPos]));
        }
      } else {
      	tempVect.push_back(std::stod(row[colPos]));
      }
    }
    iteratorMap.insert(std::make_pair(iters, tempVect));
    iters++;
  }
}

void seqSetUp::processKmerOptions() {
  setOption(runCutOffString_, "-runCutOff", "Kmer_frequencey_cutoff");
  setOption(qualRunCutOffString_, "-qualRunCutoff",
            "Kmer_frequencey_cutoff_qual");
  setOption(kLength_, "-kLength", "Kmer_length");
  if (kLength_ % 2 == 0) {
    kLength_--;
    warnings_.emplace_back("-kLength needs to be odd, not even, changing to " +
                           to_string(kLength_));
  }
  setBoolOptionFalse(kmersByPosition_, "-kAnywhere", "Kmers_at_positions");
  setBoolOptionFalse(checkKmers_, "-noKmers", "No_Kmer_checking");
  setOption(expandKmerPos_, "-expandKmerPos", "expandKmerPos");
  setOption(expandKmerSize_, "-expandKmerSize", "expandKmerSize");
}

void seqSetUp::processScoringPars() {
  setOption(local_, "-local", "Local_alignment");
  setOption(countEndGaps_, "-countEndGaps", "CountEndGaps");
  setBoolOptionFalse(weightHomopolymers_,
                     "-noHomopolymerWeighting",
                     "HomopolymerWeighting");
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
void seqSetUp::processRegKmers(){

  setOption(regKmers_, "-regKmers", "regularKmerAnalysis");
}
void seqSetUp::processSkipOnNucComp(){
  setOption(skipOnLetterCounterDifference_, "-skip",
                  "skipOnLetterCounterDifference");
  setOption(fractionDifferenceCutOff_, "-skipCutOff",
                  "fractionDifferenceCutOff");
}
void seqSetUp::processCondensedCollapse(){
  setOption(condensedCollapse_, "-condensedCollapse", "condensedCollapse");
}
void seqSetUp::processAdjustHRuns(){
  setOption(adjustHomopolyerRuns_,
            "-adjustHomopolyerRuns",
            "AdjustHomopolyerRunsToBeSameQual");
}

bool seqSetUp::processReadInNames(bool required) {
  // still need if multiple files given, like -fasta and -fastq at the same time
	VecStr readInFormats {"-sff", "-sffBin", "-fasta", "-fastq", "-bam"};
	VecStr foundInFormats;
	//process out format
	if(commands_.containsFlagCaseInsensitive("-fasta")){
		ioOptions_.inFormat_ = "fasta";
		ioOptions_.outFormat_ = "fasta";

		foundInFormats.emplace_back("-fasta");
	}else if(commands_.containsFlagCaseInsensitive("-sff")){
		ioOptions_.inFormat_ = "sff";
		ioOptions_.outFormat_ = "fastq";
		foundInFormats.emplace_back("-sff");
	}else if(commands_.containsFlagCaseInsensitive("-sffBin")){
		ioOptions_.inFormat_ = "sffbin";
		ioOptions_.outFormat_ = "fastq";
		foundInFormats.emplace_back("-sffBin");
	}else if(commands_.containsFlagCaseInsensitive("-bam")){
		ioOptions_.inFormat_ = "bam";
		ioOptions_.outFormat_ = "fastq";
		foundInFormats.emplace_back("-bam");
	}else if(commands_.containsFlagCaseInsensitive("-fastq")){
		ioOptions_.inFormat_ = "fastq";
		ioOptions_.outFormat_ = "fastq";
		foundInFormats.emplace_back("-fastq");
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
	            << bib::bashCT::red << "-fasta, -fasta/-qual, -stub, -fastq"
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
		setOption(ioOptions_.outFilename_, "-out", "Name of the out sequence file");
		return true;
	}


/*
  if (commands_.containsFlagCaseInsensitive("-fasta")) {
    ioOptions_.firstName_ = commands_["-fasta"];
    if (!commands_.containsFlagCaseInsensitive("-qual") &&
        !commands_.containsFlagCaseInsensitive("-flow")) {
      ioOptions_.inFormat_ = "fasta";
      validOptions_.push_back("-fasta");
    } else {
      if (commands_.containsFlagCaseInsensitive("-qual")) {
        ioOptions_.secondName_ = commands_["-qual"];
        if (commands_.containsFlagCaseInsensitive("-flow")) {
          ioOptions_.thirdName_ = commands_["-flow"];
          ioOptions_.inFormat_ = "fastaQualFlow";
          validOptions_.push_back("-fasta");
          validOptions_.push_back("-qual");
          validOptions_.push_back("-flow");
        } else {
          ioOptions_.inFormat_ = "fastaQual";
          validOptions_.push_back("-fasta");
          validOptions_.push_back("-qual");
        }
      } else if (commands_.containsFlagCaseInsensitive("-flow")) {
        ioOptions_.secondName_ = commands_["-flow"];
        ioOptions_.inFormat_ = "fastaFlow";
        validOptions_.push_back("-fasta");
        validOptions_.push_back("-flow");
      }
    }
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = ioOptions_.inFormat_;
    }

  } else if (commands_.containsFlagCaseInsensitive("-stub")) {
    ioOptions_.firstName_ = commands_["-stub"];
    if (ioOptions_.firstName_.find(".fasta") == std::string::npos) {
      ioOptions_.firstName_ = ioOptions_.firstName_ + ".fasta";
    }
    ioOptions_.secondName_ = ioOptions_.firstName_ + ".qual";

    if (commands_.containsFlagCaseInsensitive("-flow")) {
      ioOptions_.thirdName_ = commands_["-flow"];
      ioOptions_.inFormat_ = "fastaQualFlow";
      validOptions_.push_back("-stub");
      validOptions_.push_back("-flow");
    } else {
      ioOptions_.inFormat_ = "fastaQual";
      validOptions_.push_back("-stub");
    }
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = ioOptions_.inFormat_;
    }
  } else if (commands_.containsFlagCaseInsensitive("-fastq")) {
    ioOptions_.firstName_ = commands_["-fastq"];
    if (commands_.containsFlagCaseInsensitive("-flow")) {
      ioOptions_.secondName_ = commands_["-flow"];
      ioOptions_.inFormat_ = "fastqFlow";

      validOptions_.push_back("-fastq");
      validOptions_.push_back("-flow");
    } else {
      ioOptions_.inFormat_ = "fastq";

      validOptions_.push_back("-fastq");
    }
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = ioOptions_.inFormat_;
    }
  } else if (commands_.containsFlagCaseInsensitive("-bam")) {
    ioOptions_.firstName_ = commands_["-bam"];
    ioOptions_.inFormat_ = "bam";
    validOptions_.push_back("-bam");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fastq";
    }
  } else if (commands_.containsFlagCaseInsensitive("-flow")) {
    ioOptions_.firstName_ = commands_["-flow"];
    if (commands_.containsFlagCaseInsensitive("-qual")) {
      ioOptions_.secondName_ = commands_["-qual"];
      ioOptions_.inFormat_ = "flowQual";
      validOptions_.push_back("-flow");
      validOptions_.push_back("-qual");
    } else {
      ioOptions_.inFormat_ = "flow";
      validOptions_.push_back("-flow");
    }
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      if (ioOptions_.inFormat_ == "flowQual") {
        ioOptions_.outFormat_ = "fastaQualFlow";
      } else {
        ioOptions_.outFormat_ = ioOptions_.inFormat_;
      }
    }
  } else if (commands_.containsFlagCaseInsensitive("-pyrodata")) {
    ioOptions_.firstName_ = commands_["-pyrodata"];
    ioOptions_.inFormat_ = "pyroData";
    ioOptions_.outFormat_ = "fasta";
    validOptions_.push_back("-pyroData");
  } else if (commands_.containsFlagCaseInsensitive("-sff")) {
    ioOptions_.firstName_ = commands_["-sff"];
    ioOptions_.inFormat_ = "sff";
    validOptions_.push_back("-sff");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fastq";
    }
  } else if (commands_.containsFlagCaseInsensitive("-sffbin")) {
    ioOptions_.firstName_ = commands_["-sffbin"];
    ioOptions_.inFormat_ = "sffbin";
    validOptions_.push_back("-sffbin");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fastq";
    }
  } else if (commands_.containsFlagCaseInsensitive("-shorah")) {
    ioOptions_.firstName_ = commands_["-shorah"];
    ioOptions_.inFormat_ = "shorah";
    validOptions_.push_back("-shorah");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fasta";
    }
  } else if (commands_.containsFlagCaseInsensitive("-shorahOld")) {
    ioOptions_.firstName_ = commands_["-shorahOld"];
    ioOptions_.inFormat_ = "shorahOld";
    validOptions_.push_back("-shorahOld");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fasta";
    }
  } else if (commands_.containsFlagCaseInsensitive("-clustal")) {
    ioOptions_.firstName_ = commands_["-clustal"];
    ioOptions_.inFormat_ = "clustal";
    validOptions_.push_back("-clustal");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fasta";
    }
  } else if (commands_.containsFlagCaseInsensitive("-clustal-ng")) {
    ioOptions_.firstName_ = commands_["-clustal-ng"];
    ioOptions_.inFormat_ = "clustal-ng";
    validOptions_.push_back("-clustal-ng");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fasta";
    }
  } else if (commands_.containsFlagCaseInsensitive("-raw")) {
    ioOptions_.firstName_ = commands_["-raw"];
    ioOptions_.inFormat_ = "raw";
    validOptions_.push_back("-raw");
    if (commands_.containsFlagCaseInsensitive("-outformat")) {
      ioOptions_.outFormat_ = commands_["-outformat"];
    } else {
      ioOptions_.outFormat_ = "fasta";
    }
  } else {
    if (required) {
      std::stringstream tempOut;
      tempOut << boldBlackText("Did not find a recognizable read in option")
              << std::endl;
      tempOut << boldBlackText("options include: ") +
                     boldText("-fasta, -fasta/-qual, -stub, -fastq", "31")
              << std::endl;
      tempOut << "Command line arguments" << std::endl;
      writeOutCommandLineArguments(commands_.arguments_, tempOut);
      addOtherVec(warnings_, streamToVecStr(tempOut));
      failed_ = true;
    }
    return false;
  }
  pars_.addParameter("DataFile1", ioOptions_.firstName_, true);
  pars_.addParameter("DataFile2", ioOptions_.secondName_, true);
  pars_.addParameter("DataFile3", ioOptions_.thirdName_, true);
  pars_.addParameter("DataFileInFormat", ioOptions_.inFormat_, true);
  pars_.addParameter("OutFormat", ioOptions_.outFormat_, true);

  return true;*/
}

void seqSetUp::processDirectoryOutputName(const std::string& defaultName,
                                          bool mustMakeDirectory) {
  directoryName_ = "./";
  if (setOption(directoryName_, "-dout", "OutDirectoryName")) {
    if (!failed_) {
      directoryName_ =
          bib::files::makeDir("./", replaceString(directoryName_, "./", ""));
    }
  } else {
    if (mustMakeDirectory && !failed_) {
      directoryName_ = bib::files::makeDir("./", defaultName);
    }
  }
  /*
  bool makeDefault = false;
  if(setOption(makeDefault, "-makeDir,-makeDirectory", "MakeDefaultDirectory")){
    directoryName_ = makeDirectory("./", defaultName);
  }*/
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
  setOption(ioOptions_.processed_, "-processed", "Processed");
  setOption(ioOptions_.outFilename_, "-out", "OutFilename");
  setOption(ioOptions_.lowerCaseBases_, "-lower", "LowerCaseBaseHandling");
  setOption(ioOptions_.removeGaps_, "-removeGaps", "RemoveGaps");
  setBoolOptionFalse(ioOptions_.includeWhiteSpaceInName_,
                     "-noSpaceInName,-removeWhiteSpace,-nowhitespace",
                     "NoWhiteSapceInNmae");
  processWritingOptions();
  return passed;
}

void seqSetUp::processWritingOptions() {
  setOption(ioOptions_.overWriteFile_, "-overWrite",
            "OverWriteExistingFiles");
  setOption(ioOptions_.exitOnFailureToWrite_, "-exitOnFailureToWrite",
            "ExitOnFailureToWrite");
}

bool seqSetUp::processRefFilename(bool required) {
  setOption(refProcessed_, "-refProcessed",
            "ReferenceHasAbundanceInfo");
  if (commands_.containsFlagCaseInsensitive("-refFastq") ||
      commands_.containsFlagCaseInsensitive("-expectFastq")) {
    refFormat_ = "fastq";
  } else if (commands_.containsFlagCaseInsensitive("-ref") ||
             commands_.containsFlagCaseInsensitive("-expect")) {
    refFormat_ = "fasta";
  }
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

    inputSeq = textFileReader::getFirstLine(inputSeq);
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
    	std::string nextLine = "";
    	getline(inFile, nextLine);
    	getline(inFile, nextLine);
      seqObj_ = readObject(seqInfo(inputSeq.substr(1), nextLine));
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
  return setOption(debug_, "-d,--debug", "Debug");
}

bool seqSetUp::processQuiet() {
  return setOption(quiet_, "-q,--quiet", "quiet");
}
void seqSetUp::processAlignerDefualts() {
  processGap();
  processQualThres();
  processKmerOptions();
  processScoringPars();
  processAlnInfoInput();
  setOption(eventBased_, "-eventBased", "Event Based Comparison");
}
void seqSetUp::processAlnInfoInput() {
  if (setOption(alnInfoDirName_, "-alnInfoDir", "alnInfoDirName")) {
    if (!setOption(outAlnInfoDirName_, "-outAlnInfoDir",
                   "alnInfoDirName")) {
      outAlnInfoDirName_ = alnInfoDirName_;
    }
    writingOutAlnInfo_ = true;
  } else {
    if (setOption(outAlnInfoDirName_, "-outAlnInfoDir",
                  "alnInfoDirName")) {
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
void seqSetUp::printKmerUsage(std::ostream& out) {
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
void seqSetUp::printAlignmentUsage(std::ostream& out) {
  // std::stringstream tempOut;
  out << bib::bashCT::bold <<"Alignment options"<< bib::bashCT::reset  << std::endl;
  out << "1) -gap [option]: Gap penalty, given in the format 7,0.5 "
         "where 7 is the gap open penalty and 0.5 is the gap extension"
      << std::endl;
  out << "2) -gapLeft [option]: Gap penalty for putting gaps at the beginning"
         " of the sequence, same format as -gap" << std::endl;
  out << "3) -gapLeft [option]: Gap penalty for putting gaps at the end"
         " of the sequence, same format as -gap" << std::endl;
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
  out << "4) -condensedCollapse : When collapsing only on single base "
         "indels only compare sequences against sequences with the "
         "same homopolymer run" << std::endl;
  // out << cleanOut(tempOut.str(), width_, indent_);
}

}  // namespace bib
