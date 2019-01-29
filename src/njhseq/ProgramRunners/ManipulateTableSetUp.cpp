#include "ManipulateTableSetUp.hpp"
#include "njhseq/IO.h"
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//

namespace njhseq {

ManipulateTableSetUp::ManipulateTableSetUp(int argc, char *argv[]) : njh::progutils::ProgramSetUp(argc, argv) {
  initializeDefaults();
}
ManipulateTableSetUp::ManipulateTableSetUp(const njh::progutils::CmdArgs &inputCommands)
    : njh::progutils::ProgramSetUp(inputCommands) {
  initializeDefaults();
}

void ManipulateTableSetUp::initializeDefaults() {
	ioOptions_.in_.inFilename_ = "STDIN";
	ioOptions_.hasHeader_ = false;
	ioOptions_.inDelim_ = "whitespace";
	ioOptions_.outDelim_ = "whitespace";
	directoryName_ = "";
	sortByColumn_ = "";
	decending_ = false;
	addHeader_ = false;
	ioOptions_.outDelim_ = "";
	ioOptions_.out_.outFilename_ = "";
	ioOptions_.out_.overWriteFile_ = false;
	ioOptions_.out_.exitOnFailureToWrite_ = false;
	ioOptions_.outOrganized_ = false;
}

void ManipulateTableSetUp::setUpCatOrganized() {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "catOrganized" << std::endl;
    tempOut << "Will take in a table and output a column adjusted table to "
               "be more human readable" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -file [option]: The name of the file to be read"
            << std::endl;
    tempOut << "Optional commands" << std::endl;
    tempOut << "1) -delim [option]: The name of the file to be read, use "
               "'-delim tab' for tab delimited files, defaults to whitespace"
            << std::endl;
    tempOut << "2) -out [option]: A name of an out file, defaults to std::cout"
            << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  ioOptions_.outOrganized_ = true;
  processDefaultProgram(true);
  finishSetUp(std::cout);
}
void ManipulateTableSetUp::setUpChangeDelim() {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "changeDelim" << std::endl;
    tempOut
        << "Will take a table with certain delimiter and will change it to "
           "another delimiter given by -outDelim, default behavior is to take "
           "a whitespace delimited file and change it to tab delimited"
        << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -file [option]: The name of the file to be read"
            << std::endl;
    tempOut << "Optional commands" << std::endl;
    tempOut << "1) -delim [option]: The delimiter for the file, use "
               "'-delim tab' for tab delimited files, defaults to whitespace"
            << std::endl;
    tempOut << "2) -out [option]: A name of an out file, output will go to "
               "standard out if no name is given" << std::endl;
    tempOut << "3) -outDelim [option]: A new delimiter, will default to tab "
               "delimited" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  ioOptions_.outDelim_ = "\t";
  processDefaultProgram(true);
  finishSetUp(std::cout);
}

void ManipulateTableSetUp::setUpSortTable() {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "sortTable" << std::endl;
    tempOut
        << "Will take a table will sort it on the column given by -sortByColumn"
        << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -file [option]: The name of the file to be read"
            << std::endl;
    tempOut << "2) -sortByColumn [option]: The name of the column to sort by"
            << std::endl;
    tempOut << "Optinal commands" << std::endl;
    tempOut << "1) -delim [option]: The name of the file to be read, use "
               "'-delim tab' for tab delimited files, defaults to whitespace"
            << std::endl;
    tempOut << "2) -out [option]: A name of an outfile, defaults to tempOut"
            << std::endl;
    tempOut << "3) -header : Whether the file has a header" << std::endl;
    tempOut << "4) -decending : Whether the sort should be decending rather "
               "than ascending" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processFileName(true);
  processNonRquiredDefaults();
  processSorting(true);
  finishSetUp(std::cout);
}

void ManipulateTableSetUp::setUpExtract(
    bool &extractingByMultipleNames, VecStr &columns,
    bool &extractingByColumnElements, std::string &column, std::string &element,
    std::string &elementContains, std::string &elementStartsWith,
    bool &extractingColumnsLoose, std::string &colStartsWith,
    std::string &colContains, bool &getUniqueRows) {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "extract" << std::endl;
    tempOut << "Will take a table and extract specific columns or rows "
               "depending on the extraction parameters" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -file [option]: The name of the file to be read"
            << std::endl;
    tempOut << "Optinal commands" << std::endl;
    tempOut << "1) -delim [option]: The name of the file to be read, use "
               "'-delim tab' for tab delimited files, defaults to whitespace"
            << std::endl;
    tempOut << "2) -out [option]: A name of an outfile, defaults to std::cout"
            << std::endl;
    tempOut << "3) -header : Whether the file has a header" << std::endl;
    tempOut << "4) -sortByColumn [option]: The name of the column to sort by"
            << std::endl;
    tempOut << "5) -decending : Whether the sort should be decending rather "
               "than ascending" << std::endl;
    tempOut << "ExtractionOptions" << std::endl;
    tempOut << "1) -extractColumns [options]: A list of columns to be "
               "extract given separated by a comma ex. Col1,Col2,Col3"
            << std::endl;
    tempOut
        << "2) -extractColumn [option]: The name of the column to be searched"
        << std::endl;
    tempOut << "\t -element [option]: Extract rows that contain this element "
               "in the column name given above" << std::endl;
    tempOut << "\t -elementContains [option]: Extract rows that that element "
               "in the column name given above contains this" << std::endl;
    tempOut << "\t -elementStartsWith [option]: Extract rows that the "
               "elements from the column name given above that starts with "
               "this" << std::endl;
    tempOut << "3) -colStartsWith [option]: Extract the columns that the "
               "name starts with this" << std::endl;
    tempOut << "4) -colContains [option]: Extract the columns that the name "
               "contains with this" << std::endl;
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processFileName();
  processNonRquiredDefaults();
  processSorting();
  std::string extractColumnsStrings = "";
  if (setOption(extractColumnsStrings, "-extractColumns",
                "NamesOfColumnsToExtract")) {
    columns = tokenizeString(extractColumnsStrings, ",");
    extractingByMultipleNames = true;
  }
  if (setOption(column, "-extractColumn", "ColumnNameToBeSearched")) {
    extractingByColumnElements = true;
    if (!setOption(element, "-element", "ElementToFind")) {
      if (!setOption(elementContains, "-elementContains", "ElementContains")) {
        setOption(elementStartsWith, "-elementStartsWtih", "ElementStartsWith");
      }
    }
  }
  if (setOption(colStartsWith, "-colStartsWith,-columnStartsWith",
                "ColumnStartsWtih")) {
    extractingColumnsLoose = true;
  }
  if (setOption(colContains, "-colContains,-columnContains",
                "ColumnContains")) {
    extractingColumnsLoose = true;
  }
  setOption(getUniqueRows, "-getUniqueRows", "GetUniqueRows");
  finishSetUp(std::cout);
}

void ManipulateTableSetUp::setUpTrimContent(std::string &trimAt) {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "trimContent" << std::endl;
    tempOut << "Will take a table and trim each element in it at the first "
               "occurrence of -trimAt if it contains it" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) --file [option]: The name of the file to be read"
            << std::endl;
    tempOut << "1) --trimAt [option]: The pattern to trim all the elements in "
               "the table at" << std::endl;
    printDefaultOptinalOptions(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processDefaultProgram(true);
  setOption(trimAt, "--trimAt", "TrimElementsAt", true);
  finishSetUp(std::cout);
}
void ManipulateTableSetUp::setUpGetStats(std::string &trimAt) {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "getStats" << std::endl;
    tempOut << "Will take a table and get the min,max,mean,sum, and std on "
               "the numeric columns in the table" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -file [option]: The name of the file to be read"
            << std::endl;
    tempOut << "Optional Trim Option" << std::endl;
    tempOut << "1) -trimAt [option]: The pattern to trim all the elements in "
               "the table at before doing the stats, defaults to ( to get "
               "rid of any (%) patterns" << std::endl;
    printDefaultOptinalOptions(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processDefaultProgram();
  setOption(trimAt, "-trimAt", "TrimElementsAt");
  finishSetUp(std::cout);
}

void ManipulateTableSetUp::setUpSplitTable(std::string &column) {
  processDefaultProgram();
  setOption(column, "--column", "Column To Split On, a seperate table will be created for each unique value in this column", true);
  processDirectoryOutputName("./", false);
  finishSetUp(std::cout);
}

void ManipulateTableSetUp::setUpAggregateTable(std::string &functionName,
                                               bool &advance,
                                               std::string &columnName) {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "aggregateTable" << std::endl;
    tempOut << "Will take a table and take all the rows with the same name "
               "in -column and perform the given -function on the other "
               "columns in the table" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -file [option]: The name of the file to be read"
            << std::endl;
    tempOut << "2) -column [option]: The name of the column to aggregate on"
            << std::endl;
    tempOut << "Additional Optional Options" << std::endl;
    tempOut
        << "1) -advance : Default behavior to ignore any non-numeric columns "
           "but with -advance, aggregation will be perfromed by sorting on all "
           "non-numeric columns and not just the -columnName given"
        << std::endl;
    tempOut << "2) -function [option] :A name of the function to perform on "
               "aggregation, defaults to mean, options are mean, sum, min, "
               "max, and std (standard deviation" << std::endl;
    printDefaultOptinalOptions(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processDefaultProgram();
  setOption(columnName, "-colName,-columnName,-column", "ColumnName", true);
  setOption(functionName, "-function,-functionName", "StatFunction");
  setOption(advance, "-advance", "AdvacnedAgrregation");
  finishSetUp(std::cout);
}
void ManipulateTableSetUp::setUpPivotTable(std::string &columnName,
                                           std::string &matchColumn) {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "pivotTable" << std::endl;
    tempOut << "Will split the table on -column and then put the table side "
               "by side and match their rows on -matchColumn" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -file [option]: The name of the file to be read"
            << std::endl;
    tempOut
        << "2) -column [option]: The name of the column to split the table on"
        << std::endl;
    tempOut << "3) -matchColumn [option]: The name of the column match the "
               "split tables on" << std::endl;
    printDefaultOptinalOptions(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processDefaultProgram();
  setOption(columnName, "-colName,-columnName,-column", "ColumnName", true);
  setOption(matchColumn, "-matchColumn,-columnMatch", "Matchcolumn", true);
  finishSetUp(std::cout);
}


///////////////
bool ManipulateTableSetUp::processFileName(bool required) {
  return setOption(ioOptions_.in_.inFilename_, "--file", "Input FileName", required);
}

void ManipulateTableSetUp::processNonRquiredDefaults() {
	setOption(ioOptions_.hasHeader_, "--header", "HeaderPresent");
	if(ioOptions_.hasHeader_){
		addHeader_ = true;
	}
	setOption(ioOptions_.inDelim_, "--delim", "FileDelimiter");
	if(!addHeader_){
		setOption(addHeader_, "--addHeader", "If input doesn't have a header add a header to the output, will default to col.[COL_NUM]");
	}
	if (ioOptions_.inDelim_ == "tab") {
		ioOptions_.inDelim_ = "\t";
	}/* else if (delim_ == "whitespace") {
	 delim_ = " "; writingDistanceCheck
	 }*/
	processWriteOutOptions();
}

bool ManipulateTableSetUp::processSorting(bool required) {
  setOption(decending_, "--Descending", "Sorting In Descending Order");
  return setOption(sortByColumn_, "--sortByColumn", "Sorting By ColumnName",
                   required);
}
void ManipulateTableSetUp::processWriteOutOptions() {
  if(setOption(ioOptions_.out_.outFilename_, "--out", "Out File Name")){
  		if("" != njh::files::getExtension(ioOptions_.out_.outFilename_)){
  			ioOptions_.out_.outExtention_ = "." + njh::files::getExtension(ioOptions_.out_.outFilename_);
  		}else if("" != njh::files::getExtension(ioOptions_.in_.inFilename_)){
  			ioOptions_.out_.outExtention_ = "." + njh::files::getExtension(ioOptions_.in_.inFilename_);
  		}
  }
  setOption(ioOptions_.out_.outExtention_, "--extension", "Out File Extension");
  if (!setOption(ioOptions_.outDelim_, "--outDelim",
                 "Outfile delimiter")) {
    ioOptions_.outDelim_ = ioOptions_.inDelim_;
  }
  if (ioOptions_.outDelim_ == "tab") {
    ioOptions_.outDelim_ = "\t";
  } else if (ioOptions_.outDelim_ == "whitespace") {
    ioOptions_.outDelim_ = " ";
  }
  setOption(ioOptions_.outOrganized_, "--outOrganized",
            "Ouput Human Readable Column Adjusted");
  setOption(ioOptions_.out_.overWriteFile_, "--overWrite",
            "OverWriteFiles");
//  setOption(ioOptions_.out_.exitOnFailureToWrite_, "-exitOnFailureToWrite",
//            "ExitOnFailureToWrite");
  setOption(ioOptions_.out_.append_, "--append", "Append existing file or create new file if non exists, will not repeat header");
}

bool ManipulateTableSetUp::processDefaultProgram(bool fileRequired) {
  bool passed = processFileName(fileRequired);
  processSorting();
  processNonRquiredDefaults();
  //finishSetUp(std::cout);
  return passed;
}//miscTest

void ManipulateTableSetUp::processDirectoryOutputName(
    const std::string &defaultName, bool mustMakeDirectory) {
  if (setOption(directoryName_, "-dout", "Name of Out Directory")) {
    directoryName_ = njh::files::makeDir("./", njh::files::MkdirPar(njh::replaceString(directoryName_, "./", ""))).string();
  } else {
    if (mustMakeDirectory) {
      directoryName_ = njh::files::makeDir("./", njh::files::MkdirPar(defaultName)).string();
    } else {
      directoryName_ = "./";
    }
  }

}

void ManipulateTableSetUp::printDefaultOptinalOptions(std::ostream &out) {
  out << njh::bashCT::bold <<
    		njh::bashCT::black <<
             "Optinal commands for all ManipulateTable programs, some might "
             "not "
             "pertain to certain programs" << njh::bashCT::reset << std::endl;
  out << "1) -delim [option]: The name of the file to be read, use "
         "'-delim tab' for tab delimited files, defaults to whitespace"
      << std::endl;
  out << "2) -header : Whether the file has a header, defaults to false"
      << std::endl;

  out << "3) -sortByColumn [option]: The name of a column to sort the table by"
      << std::endl;
  out << "4) -decending : If sorting, whether the sort should be decending "
         "rather "
         "than ascending" << std::endl;
  out << njh::bashCT::bold <<
    		njh::bashCT::black <<"Optional Write Options" << njh::bashCT::reset << std::endl;
  out << "1) -out [option]: A name of an outfile, defaults to standard out"
      << std::endl;
  out << "2) -outDelim [option] : The delimiter for output, defaults to the "
         "delimiter of the in file" << std::endl;
  out << "3) -outOrganized : Whether if the columns should be aligned to human "
         "readable" << std::endl;
  out << "4) -overWrite : if writing a file whether the file should be "
         "overwritten if it already exits, defaults to false" << std::endl;
  out << "5) -exitOnFailureToWrite : if writing a file and it fails whether "
         "the program should immediately exit without continuing, defaults to "
         "false" << std::endl;
}
void ManipulateTableSetUp::setUpRBind(std::string &contains, bool &recursive) {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "rBind" << std::endl;
    tempOut << "Will take a table and get the min,max,mean,sum, and std on "
               "the numeric columns in the table" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -contains [option]: " << std::endl;

    tempOut << "Optional Option" << std::endl;
    tempOut << "1) -r " << std::endl;
    printDefaultOptinalOptions(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processDefaultProgram(false);
  setOption(contains, "--contains", "contains", true);
  setOption(recursive, "--recursive", "recursive");
  finishSetUp(std::cout);
}
void ManipulateTableSetUp::setUpCBind(std::string &contains, bool &recursive) {
  if (needsHelp()) {
    std::stringstream tempOut;
    tempOut << "cBind" << std::endl;
    tempOut << "Will take a table and get the min,max,mean,sum, and std on "
               "the numeric columns in the table" << std::endl;
    tempOut << "Required commands" << std::endl;
    tempOut << "1) -contains [option]: " << std::endl;

    tempOut << "Optional Option" << std::endl;
    tempOut << "1) -r " << std::endl;
    printDefaultOptinalOptions(tempOut);
    std::cout << cleanOut(tempOut.str(), width_, indent_);
    exit(1);
  }
  processDefaultProgram(false);
  setOption(contains, "-contains", "contains", true);
  setOption(recursive, "-r,-recrusive", "recursive");
  finishSetUp(std::cout);
}


bool ManipulateTableSetUp::processVerbose() {
  return setOption(verbose_, "-v,--verbose", "Verbose");
}

bool ManipulateTableSetUp::processDebug() {
  return setOption(debug_, "--debug", "Debug");
}

}  // namespace njh
