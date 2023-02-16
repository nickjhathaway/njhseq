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


#include "ManipulateTableRunner.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/objects/dataContainers/tables/TableReader.hpp"

#include "njhseq/objects/Meta/MetaDataInName.hpp"
#include "njhseq/programUtils/seqSetUp.hpp"
#include "njhseq/IO.h"

namespace njhseq {

ManipulateTableRunner::ManipulateTableRunner() :
		njh::progutils::ProgramRunner(
				{ addFunc("addColumn", addColumn, false),
				  addFunc("catOrganized",catOrganized, false),
					addFunc("changeDelim", changeDelim, false),
					addFunc("sortTable", sortTable, false),
					addFunc("trimContent", trimContent,false),
					addFunc("getStats", getStats, false),
					addFunc("splitTable", splitTable, false),
					addFunc("aggregateTable", aggregateTable, false),
					addFunc("pivotTable", pivotTable,false),
					addFunc("rBind", rBind, false),
					addFunc("countColumn", countColumn, false),
					addFunc("tableExtractCriteria", tableExtractCriteria, false),
					addFunc("cBind", cBind, false),
					addFunc("countRowLengths",countRowLengths, false),
					addFunc("extractColumnElementLength",extractColumnElementLength, false),
					addFunc("printCol",printCol, false),

					addFunc("tableExtractColumns", tableExtractColumns, false),
					addFunc("tableExtractElementsWithPattern", tableExtractElementsWithPattern, false),
					addFunc("tableExtractElementsStartingWith",tableExtractElementsStartingWith, false),
					addFunc("tableExtractColumnsStartsWith",tableExtractColumnsStartsWith, false),
					addFunc("tableExtractColumnsWithPattern",tableExtractColumnsWithPattern, false),
					addFunc("splitColumnContainingMeta",splitColumnContainingMeta, false),
					addFunc("roughHistogramOfColumn",roughHistogramOfColumn, false),
					addFunc("removeColumns",removeColumns, false),
				}, "ManipulateTable", "1") {
}
//



int ManipulateTableRunner::splitColumnContainingMeta(
		const njh::progutils::CmdArgs & inputCommands) {

	table::splitColWithMetaPars splitingColsPars;
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setOption(splitingColsPars.column_, "--column", "Column to split which contains meta with the following formating [key1=val1;key2=val2;]", true);
	setUp.setOption(splitingColsPars.keepMetaInColumn_, "--keepMetaInColumn", "Keep meta in the original column");
	setUp.setOption(splitingColsPars.removeEmptyColumn_, "--removeEmptyColumn", "remove original column if it becomes an empty column");
	setUp.setOption(splitingColsPars.prefixWithColName_, "--prefixWithColName", "Prefix new columns with column name");
	setUp.processFileName(true);
	setUp.processNonRquiredDefaults();
	splitingColsPars.sorting_ = setUp.processSorting();
	splitingColsPars.descending_ = setUp.decending_;
	splitingColsPars.sortCol_ = setUp.sortByColumn_;
	setUp.finishSetUp(std::cout);
	table inTab(setUp.ioOptions_);

	table outTab = table::splitColWithMeta(inTab, splitingColsPars);
	outTab.hasHeader_ = setUp.addHeader_;
	setUp.ioOptions_.hasHeader_ = setUp.addHeader_;
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}

//
int ManipulateTableRunner::extractColumnElementLength(
		const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	std::string outFilename = "";
	uint64_t length = 6;
	bool longer = false;
	std::string columnName = "";
	setUp.setOption(filename, "-i,--inFile", "In Filename", true);
	setUp.setOption(columnName, "-c,--columnName",
			"Column Name To compare against", true);
	setUp.setOption(outFilename, "-o,--outFile", "Out File Name");
	setUp.processWritingOptions();
	setUp.setOption(length, "-l,--len", "Length to compare column elements to");
	setUp.setOption(longer, "--longer",
			"Extract row if column element is less than --len");
	setUp.finishSetUp(std::cout);
	table inTab(filename, "\t", true);
	table outTab;
	if (longer) {
		outTab = inTab.extractByComp(columnName,
				[&length](const std::string & str) {return str.length() > length;});
	} else {
		outTab = inTab.extractByComp(columnName,
				[&length](const std::string & str) {return str.length() < length;});
	}
	TableIOOpts outOpts(
			OutOptions(outFilename, ".tab.txt", "tab", false,
					setUp.pars_.ioOptions_.out_.overWriteFile_, "false"), "\t",
			outTab.hasHeader_);
	if (outFilename == "") {
		outOpts.outOrganized_ = true;
	}
	outTab.outPutContents(outOpts);
	return 0;
}

int ManipulateTableRunner::countRowLengths(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setOption(setUp.ioOptions_.inDelim_, "--delim", "The delimiter for the rows");
	setUp.setOption(setUp.ioOptions_.hasHeader_, "--header", "File Contains Header");
	setUp.setOption(setUp.ioOptions_.in_.inFilename_, "--file", "Filename", true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	InputStream inFile(setUp.ioOptions_.in_);
	std::string line;
	uint32_t rowNumber = 0;
	table outTab { VecStr { "Row", "length", "elements" } };
	while (njh::files::crossPlatGetline(inFile, line)) {
		if (setUp.verbose_){
			std::cout << "On " << rowNumber << std::endl;
		}
		if (setUp.ioOptions_.hasHeader_ || (rowNumber > 0)) {
			uint32_t index = 0;
			if (setUp.ioOptions_.inDelim_ == "") {
				std::stringstream ss(line);
				std::string out;
				while (!ss.eof()) {
					ss >> out;
					++index;
				}
			} else {
				auto toks = tokenizeString(line, setUp.ioOptions_.inDelim_);
				index = toks.size();
			}
			outTab.content_.emplace_back(toVecStr(rowNumber, line.length(), index));
		}
		++rowNumber;
	}
	outTab.outPutContents(std::cout, "\t");
	return 0;
}

int ManipulateTableRunner::tableExtractCriteria(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.description_ = "Extract table contents by filtering on a column, by default will take values higher than cut off, --lessThan to take values less than cut off";
	std::string columnName = "";
	double cutOff = 1;
	bool lessThan = false;
	setUp.processDefaultProgram(true);
	setUp.setOption(columnName, "--columnName", "Name of the column to extract on using cutOff", true);
	setUp.setOption(cutOff, "--cutOff", "Cut off of column to be extract, not inclusive");
	setUp.setOption(lessThan, "--lessThan", "Take numbers less than value in cutOff flag");
	setUp.finishSetUp(std::cout);

	TableReader inTable(setUp.ioOptions_);
	inTable.header_.containsColumn(columnName);
	OutputStream out(setUp.ioOptions_.out_);
	if(setUp.ioOptions_.hasHeader_){
		inTable.header_.outPutContents(out, setUp.ioOptions_.outDelim_);
	}
	std::function<bool(const std::string & str)> comp;
	if (lessThan) {
		comp = [cutOff](const std::string & str) {
			double numValue = std::stod(str);
			return numValue < cutOff;
		};
	} else {
		comp = [cutOff](const std::string & str) {
			double numValue = std::stod(str);
			return numValue > cutOff;};
	}

	VecStr row;
	auto colPos = inTable.header_.getColPos(columnName);
	while(inTable.getNextRow(row)){
		if(comp(row[colPos])){
			out << njh::conToStr(row, setUp.ioOptions_.outDelim_) << "\n";
		}
	}

	return 0;
}


int ManipulateTableRunner::removeColumns(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	VecStr columns;
	setUp.processDefaultProgram(true);
	setUp.setOption(columns, "--columns","Columns to remove", true);
	setUp.finishSetUp(std::cout);
	table outTab(setUp.ioOptions_);

	outTab.checkForColumnsThrow(columns, __PRETTY_FUNCTION__);
	for(const auto & col : columns){
		outTab.deleteColumn(col);
	}
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}


int ManipulateTableRunner::addColumn(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string columnName = "";
	VecStr element;
	std::string elementStr = "";
	setUp.processDefaultProgram(true);
	setUp.setOption(columnName, "--newColumnName","Name of the new Column to add to table", true);
	setUp.setOption(elementStr, "--element","What to Add to the Table Under Column, can be several comma sep values or just one",true);
	setUp.finishSetUp(std::cout);
	table outTab(setUp.ioOptions_);
	auto toks = tokenizeString(elementStr, ",");
	auto colToks = tokenizeString(columnName, ",");
	if(colToks.size() > 1){
		if(toks.size() != colToks.size()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error, when adding multiple columns the number of elements need to match number of new columns" << "\n"
					<< "New Columns #: " << colToks.size() << ", New Elements #: " << toks.size() <<"\n";
			throw std::runtime_error{ss.str()};
		}
		outTab.addSingleValueColumns(toks, colToks);
	} else {
		if (toks.size() > 1) {
			if (outTab.content_.size() % toks.size() != 0) {
				std::stringstream ss;
				ss << njh::bashCT::red << "Error, table has "
						<< outTab.content_.size() << " and the size of adding elements "
						<< njh::conToStr(elementStr, ",") << " doesn't fit into it"
						<< njh::bashCT::reset << std::endl;
				throw std::runtime_error{ss.str()};
			} else {
				outTab.addColumn(
						repeatVector(toks, { outTab.content_.size() / toks.size() }),
						columnName);
			}
		} else if (toks.size() == 1) {
			outTab.addColumn(toks, columnName);
		}
	}

	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}



int ManipulateTableRunner::roughHistogramOfColumn(
		const njh::progutils::CmdArgs & inputCommands) {
	uint32_t binNumbers = 10;
	uint32_t width = 50;
	ManipulateTableSetUp setUp(inputCommands);
	setUp.processDefaultProgram(true);
	std::string columnName = "";
	setUp.setOption(columnName, "--columnName", "columnName", true);
	setUp.setOption(binNumbers, "--binNumbers", "Bin Numbers");
	if(binNumbers < 2){
		setUp.failed_ = true;
		setUp.addWarning("Number of bins need to be 2 or greater");
	}
	setUp.setOption(width,       "--width", "width");
	if(0 == width){
		setUp.failed_ = true;
		setUp.addWarning("Width can't be zero");
	}
	setUp.finishSetUp(std::cout);

	table inTab(setUp.ioOptions_);
	auto colPos = inTab.getColPos(columnName);

	std::vector<double> columnValues;
	for(const auto & row : inTab.content_){
		if(row[colPos] != "*" && njh::strToLowerRet(row[colPos]) != "na" && njh::strToLowerRet(row[colPos]) != "nan"){
			columnValues.emplace_back(njh::StrToNumConverter::stoToNum<double>(row[colPos]));
		}
	}
	auto min = vectorMinimum(columnValues);
	auto max = vectorMaximum(columnValues);

	auto step = (max - min)/binNumbers;

	struct Bin{
		Bin(double min, double max): min_(min), max_(max){

		}
		double min_;
		double max_;
		uint32_t count_ = 0;

		std::string genId()const{
			return njh::pasteAsStr(min_, "-", max_);
		}
	};

	std::vector<Bin> bins;
	for(uint32_t binNum = 0; binNum < binNumbers; ++binNum){
		if(0 == binNum){
			bins.emplace_back(Bin{min, min + step});
		}else if(binNum == binNumbers - 1){
			bins.emplace_back(Bin{min + step * binNum, max + 0.001});
		}else{
			bins.emplace_back(Bin{min + step * binNum, min + step * (binNum + 1)});
		}
	}

	for(const auto & v : columnValues){
		for( auto & b : bins){
			if(v >= b.min_ && v < b.max_){
				++b.count_;
				break;
			}
		}
	}

	uint32_t maxNameLen = 0;
	double maxCount = 0;
	for(const auto & b : bins){
		if(b.genId().size() > maxNameLen){
			maxNameLen = b.genId().size();
		}
		if(b.count_ > maxCount){
			maxCount = b.count_;
		}
	}

	for(const auto & b : bins){
		uint32_t outWidth = std::round((b.count_/maxCount) * width);
		std::cout << std::string(maxNameLen - b.genId().size(), ' ') << b.genId() << ":" << std::string(outWidth,'*') << std::string(width - outWidth, ' ') << " " << b.count_ << std::endl;
	}

	return 0;
}

int ManipulateTableRunner::countColumn(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.processDefaultProgram(true);
	VecStr columnName;
	setUp.setOption(columnName, "--columnName", "columnName", true);
	setUp.finishSetUp(std::cout);
	TableReader inTabReader(setUp.ioOptions_);

	if(!("STDIN" == setUp.ioOptions_.in_.inFilename_ && !setUp.ioOptions_.hasHeader_)){
		inTabReader.header_.checkForColumnsThrow(columnName, __PRETTY_FUNCTION__);
	}

	setUp.ioOptions_.out_.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);

	table ret;
	VecStr row;
	if (columnName.size() == 1) {
		uint32_t colPos = std::numeric_limits<uint32_t>::max();
		if("STDIN" == setUp.ioOptions_.in_.inFilename_ && !setUp.ioOptions_.hasHeader_){
			std::regex pat("col.([0-9]+)");
			std::smatch match;
			if(std::regex_match(columnName.front(), match, pat)){
				colPos = njh::StrToNumConverter::stoToNum<uint32_t>(match[1]);
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error "
						<< "when reading from standard in (STDIN) with no header, need to supply column name in pattern col.[0-9]+, not "
						<< columnName.front() << "\n";
				throw std::runtime_error{ss.str()};
			}
		}else{
			colPos = inTabReader.header_.getColPos(columnName.front());
		}

		std::unordered_map<std::string, uint64_t> counts;
		while(inTabReader.getNextRow(row)){
			++counts[row[colPos]];
		}
		ret = table(counts, VecStr{columnName.front(), "count"});
	} else {
		std::map<std::string, uint32_t> counts;
		std::vector<uint32_t> colPositions(columnName.size());
		for(const auto pos : iter::range(columnName.size())){
			colPositions[pos] = inTabReader.header_.getColPos(columnName[pos]);
		}
		while(inTabReader.getNextRow(row)){
			std::string codedName = "";
			for(const auto & pos : colPositions){
				if(pos == colPositions.back()){
					codedName+= row[pos];
				}else{
					codedName+= row[pos] + "SPLITONTHIS";
				}
			}
			++counts[codedName];
		}
		VecStr outColumnNames;
		for(const auto & pos : colPositions){
			outColumnNames.emplace_back(inTabReader.header_.columnNames_[pos]);
		}
		outColumnNames.emplace_back("count");
		ret = table(outColumnNames);
		for(const auto & c : counts){
			auto toks = tokenizeString(c.first,"SPLITONTHIS");
			toks.emplace_back(estd::to_string(c.second));
			ret.content_.emplace_back(toks);
		}
	}
	if (setUp.sortByColumn_ != "") {
		ret.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	ret.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::catOrganized(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setUpCatOrganized();
	table inTab(setUp.ioOptions_);
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::changeDelim(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setUpChangeDelim();
	table inTab(setUp.ioOptions_);
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::sortTable(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setUpSortTable();
	table inTab(setUp.ioOptions_);
	inTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}




int ManipulateTableRunner::tableExtractColumns(const njh::progutils::CmdArgs & inputCommands){
	ManipulateTableSetUp setUp(inputCommands);
	VecStr extractColumns;
	bool getUniqueRows = false;
	setUp.processFileName();
	setUp.processNonRquiredDefaults();
	setUp.processSorting();
  std::string extractColumnsStrings = "";
  setUp.setOption(extractColumns, "--columns",
                  "Names Of Columns To Extract, either comma separated input or a file with each line a column name",
									true);
  setUp.setOption(getUniqueRows, "-getUniqueRows", "GetUniqueRows");
  setUp.finishSetUp(std::cout);



  if(!getUniqueRows && "" == setUp.sortByColumn_){
  	// if not sorting and not getting unique columns, just use tab reader
  	TableReader tabReader(setUp.ioOptions_);
  	tabReader.header_.checkForColumnsThrow(extractColumns,__PRETTY_FUNCTION__);

  	OutputStream out(setUp.ioOptions_.out_);
  	out << njh::conToStr(extractColumns, setUp.ioOptions_.outDelim_) << "\n";
  	VecStr row;
  	std::vector<uint32_t> columnPositions;
  	for(const auto & col : extractColumns){
  		columnPositions.emplace_back(tabReader.header_.getColPos(col));
  	}
  	std::function<void(const VecStr &)>extractAndWrite;
  	if(1 == extractColumns.size()){
  		extractAndWrite = [&out,&columnPositions](const VecStr & row){
  			out << row[columnPositions[0]] << "\n";
  		};
  	}else{
  		extractAndWrite = [&out,&columnPositions,&setUp](const VecStr & row){
  			out << row[columnPositions[0]];
  			for(const auto pos : iter::range<uint32_t>(1,columnPositions.size())){
  				out << setUp.ioOptions_.outDelim_ << row[columnPositions[pos]];
  			}
  			out << "\n";
  		};
  	}
  	while(tabReader.getNextRow(row)){
  		extractAndWrite(row);
  	}
  } else {
  	table inTab(setUp.ioOptions_);
  	table outTab = inTab.getColumns(extractColumns);
  	if (setUp.sortByColumn_ != "") {
  		outTab.sortTable(setUp.sortByColumn_, setUp.decending_);
  	}
  	if (getUniqueRows) {
  		outTab = outTab.getUniqueRows();
  	}
  	outTab.outPutContents(setUp.ioOptions_);
  }
	return 0;
}

int ManipulateTableRunner::tableExtractElementsWithPattern(
		const njh::progutils::CmdArgs & inputCommands) {
	std::string column;
	std::string patStr = "";
	bool getUniqueRows = false;
	bool opposite = false;
	ManipulateTableSetUp setUp(inputCommands);

	setUp.processFileName();
	setUp.processNonRquiredDefaults();
	setUp.processSorting();
	setUp.setOption(patStr, "--patStr",
			"Pattern to match elements in column with", true);
	setUp.setOption(column, "--column", "Name of the column to search", true);
	setUp.setOption(getUniqueRows, "--getUniqueRows", "Get Unique Rows");
	setUp.setOption(opposite, "--opposite",
			"Get elements in column that doesn't match this pattern");

	setUp.finishSetUp(std::cout);

	table inTab(setUp.ioOptions_);
	table outTab;
	if (opposite) {
		outTab = inTab.getRowsNotContainingPattern(column, std::regex { patStr });
	} else {
		outTab = inTab.getRowsMatchingPattern(column, std::regex { patStr });
	}
	if (setUp.sortByColumn_ != "") {
		outTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	if (getUniqueRows) {
		outTab = inTab.getUniqueRows();
	}
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::tableExtractElementsStartingWith(const njh::progutils::CmdArgs & inputCommands){
	std::string column;
	std::string patStr = "";
	bool getUniqueRows = false;
	ManipulateTableSetUp setUp(inputCommands);

	setUp.processFileName();
	setUp.processNonRquiredDefaults();
	setUp.processSorting();
	setUp.setOption(patStr, "--patStr", "Pattern to match to the beginning of the elements in column with", true);
	setUp.setOption(column, "--column", "Name of the column to search", true);
	setUp.setOption(getUniqueRows, "--getUniqueRows", "GetUniqueRows");
	setUp.finishSetUp(std::cout);

	table inTab(setUp.ioOptions_);
	table outTab = inTab.getRowsMatchingPattern(column, std::regex{"^" + patStr});

	if (setUp.sortByColumn_ != "") {
		outTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	if (getUniqueRows) {
		outTab = inTab.getUniqueRows();
	}
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::tableExtractColumnsStartsWith(const njh::progutils::CmdArgs & inputCommands){
	std::string patStr = "";
	bool getUniqueRows = false;
	ManipulateTableSetUp setUp(inputCommands);

	setUp.processFileName();
	setUp.processNonRquiredDefaults();
	setUp.processSorting();
	setUp.setOption(patStr, "--patStr", "Pattern to match columns with", true);
	setUp.setOption(getUniqueRows, "--getUniqueRows", "GetUniqueRows");
	setUp.finishSetUp(std::cout);

	table inTab(setUp.ioOptions_);
	table outTab = inTab.getColumnsMatchingPattern(std::regex{patStr});

	if (setUp.sortByColumn_ != "") {
		outTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	if (getUniqueRows) {
		outTab = inTab.getUniqueRows();
	}
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::tableExtractColumnsWithPattern(const njh::progutils::CmdArgs & inputCommands){
	std::string patStr = "";
	bool getUniqueRows = false;
	ManipulateTableSetUp setUp(inputCommands);

	setUp.processFileName();
	setUp.processNonRquiredDefaults();
	setUp.processSorting();
	setUp.setOption(patStr, "--patStr", "Pattern to match to the beginning of columns with", true);
	setUp.setOption(getUniqueRows, "--getUniqueRows", "GetUniqueRows");
	setUp.finishSetUp(std::cout);

	table inTab(setUp.ioOptions_);
	table outTab = inTab.getColumnsMatchingPattern(std::regex{"^" + patStr});

	if (setUp.sortByColumn_ != "") {
		outTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	if (getUniqueRows) {
		outTab = inTab.getUniqueRows();
	}
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}


int ManipulateTableRunner::trimContent(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string trimAt = "";
	setUp.setUpTrimContent(trimAt);
	table inTab(setUp.ioOptions_);
	inTab.trimElementsAtFirstOccurenceOf(trimAt);
	if (setUp.sortByColumn_ != "") {
		inTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}
int ManipulateTableRunner::getStats(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string trimAt = "(";
	setUp.setUpGetStats(trimAt);
	table inTab(setUp.ioOptions_);
	inTab.trimElementsAtFirstOccurenceOf(trimAt);
	table outTable = inTab.getStatsTable();
	if (setUp.sortByColumn_ != "") {
		outTable.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	outTable.outPutContents(setUp.ioOptions_);
	return 0;
}




int ManipulateTableRunner::splitTable(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string column = "";
	setUp.setUpSplitTable(column);
	if("./" != setUp.directoryName_){
		setUp.startARunLog(setUp.directoryName_);
	}
	TableReader tabReader(setUp.ioOptions_);
	tabReader.header_.checkForColumnsThrow({column}, __PRETTY_FUNCTION__);

	MultiOutputStreamCache writer;
	VecStr line;
	while(tabReader.getNextRow(line)){
		std::string colVal = line[tabReader.header_.getColPos(column)];
		if(!writer.containsReader(colVal)){
			writer.addOutputStream(colVal, OutOptions(njh::files::make_path(setUp.directoryName_, colVal + tabReader.tabOpts_.out_.outExtention_)));
			writer.add(colVal, njh::conToStr(tabReader.header_.columnNames_, tabReader.tabOpts_.outDelim_));
		}
		writer.add(colVal, njh::conToStr(line, tabReader.tabOpts_.outDelim_));
	}
	table inTab(setUp.ioOptions_);
	return 0;
}


int ManipulateTableRunner::aggregateTable(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string columnName = "";
	std::string functionName = "mean";
	bool advance = false;
	setUp.setUpAggregateTable(functionName, advance, columnName);

	table inTab(setUp.ioOptions_);

	if (advance) {
		table outTable = inTab.aggregateAdvance(columnName, functionName);
		outTable.outPutContents(setUp.ioOptions_);
	} else {
		table outTable = inTab.aggregateSimple(columnName, functionName);
		outTable.outPutContents(setUp.ioOptions_);
	}

	return 0;
}
int ManipulateTableRunner::pivotTable(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string column = "";
	std::string matchColumn = "";
	setUp.setUpPivotTable(column, matchColumn);
	table inTab(setUp.ioOptions_);

	// split the table on the column
	auto tables = inTab.splitTableOnColumn(column);
	// then bind the split table and force the rows to match
	auto ans = table::cbind(tables, matchColumn);
	ans.outPutContents(setUp.ioOptions_);
	return 0;
}


/**@b Function to list all the files of a directory with the option to search recursively and name filtering
 *
 * @param dirName the name of the directory to search
 * @param recursive Whether the search should be recursive
 * @param contains A vector of strings that the path names must contains to be returned
 * @param levels The maximum number of levels to search
 * @return A map of the directory paths with key being the file path and the value being a bool indicating if it is a directory or not
 */
inline std::map<njh::files::bfs::path, bool> listAllFiles(
		const std::string & dirName, bool recursive,
		const std::vector<std::string>& contains,
		const std::vector<std::string>& doesNotContain,
		uint32_t levels = std::numeric_limits < uint32_t > ::max()) {
	std::map < bfs::path, bfs::path > filesGathering;
	njh::files::listAllFilesHelper(dirName, recursive, filesGathering, 1, levels);
	std::map<bfs::path, bool> files = njh::files::convertMapFnpFnpToFnpIsDir(
			filesGathering);

	if (!contains.empty()) {
		std::map<njh::files::bfs::path, bool> specificFiles;
		for (const auto & f : files) {
			if (njh::checkForSubStrs(f.first.string(), contains)
					&& !njh::checkForSubStrs(f.first.string(), doesNotContain)) {
				specificFiles.emplace(f);
			}
		}
		return specificFiles;
	}
	return files;
}



int ManipulateTableRunner::rBind(
		const njh::progutils::CmdArgs & inputCommands) {
	std::string contains = "";
	std::string doesNotContains = "";
	bool recursive = false;
	bool verbose = false;
	bool fill = false;
	bool skipNonExistFiles = false;
	std::string files = "";
	uint32_t levels = std::numeric_limits < uint32_t > ::max();
	std::string directory = "./";
	ManipulateTableSetUp setUp(inputCommands);
	setUp.description_ = "Program to bind together text files representing tables with similar information";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --files file1.tsv,file2.tsv,file3.tsv --header --delim tab --out combined.tsv --overWrite #combined a list of files");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --contains _info.txt --header --delim tab --out combined.tsv --overWrite #combine together all files in current directory that contain _info.txt in their name");
	setUp.setOption(levels, "--depth", "Depth of search");
	setUp.setOption(directory, "--directory", "Directory to search, defaults to current directory");
	setUp.setOption(verbose, "--verbose,-v", "verbose");
	setUp.setOption(doesNotContains, "--doesNotContains", "Does Not Contains");
	setUp.setOption(files, "--files", "A direct list of files to combine");
	setUp.processDefaultProgram(false);
	setUp.setOption(contains, "--contains", "contains", "" == files);
	setUp.setOption(recursive, "--recursive", "recursive");
	setUp.setOption(fill, "--fill", "fill in missing columns with NAs when combining tables");
	setUp.setOption(skipNonExistFiles, "--skipNonExistFiles", "skip Non Existant Files");


	if(setUp.needsHelp()){
		setUp.printFlags(std::cout);
		exit(1);
	}
	setUp.finishSetUp(std::cout);
	setUp.ioOptions_.out_.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);
	VecStr containsVec = tokenizeString(contains, ",");
	VecStr doesNotContainsVec = tokenizeString(doesNotContains, ",");
	std::map<njh::files::bfs::path, bool> allFiles;
	if("" != files){
		auto vals = getInputValues(files, ",");
		for(const auto & f : vals){
			allFiles.emplace(f, false);
		}
	}else{
		if (doesNotContains != "") {
			allFiles = listAllFiles(directory, recursive, containsVec, doesNotContainsVec,
					levels);
		} else {
			allFiles = njh::files::listAllFiles(directory, recursive, containsVec, levels);
		}
	}





	if(!fill){
		njh::files::bfs::path firstFileFnp;
		for (const auto &file : allFiles) {
			if (njh::files::bfs::is_directory(file.first)) {
				continue;
			}
			if (skipNonExistFiles && !bfs::exists(file.first)){
				continue;
			}
			if (0 == njh::files::bfs::file_size(file.first)) {
				continue;
			}
			firstFileFnp = file.first;
			break;
		}

		TableReader firstTable(TableIOOpts(InOptions(firstFileFnp), setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_));
		OutputStream out(setUp.ioOptions_.out_);
		VecStr row;
		if(setUp.ioOptions_.hasHeader_){
			out << njh::conToStr(firstTable.header_.columnNames_, setUp.ioOptions_.outDelim_) << '\n';
		}
		while(firstTable.getNextRow(row)){
			out << njh::conToStr(row, setUp.ioOptions_.outDelim_) << '\n';
		}
		for (const auto &file : allFiles) {
			if (verbose) {
				std::cout << file.first.string() << std::endl;
			}
			if (njh::files::bfs::is_directory(file.first)) {
				if (verbose) {
					std::cout << "Skipping directory: " << file.first.string() << std::endl;
				}
				continue;
			}
			if (skipNonExistFiles && !bfs::exists(file.first)){
				continue;
			}
			if (0 == njh::files::bfs::file_size(file.first)) {
				if (verbose) {
					std::cout << "Skipping empty file: " << file.first.string() << std::endl;
				}
				continue;
			}
			if(file.first != firstFileFnp){
				TableReader currentTable(TableIOOpts(InOptions(file.first), setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_));
				VecStr row;
				if(!std::equal(firstTable.header_.columnNames_.begin(), firstTable.header_.columnNames_.end(),
						currentTable.header_.columnNames_.begin(), currentTable.header_.columnNames_.end())){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "header for " << file.first << " doesn't match other columns"<< "\n";
					ss << "expected header: " << njh::conToStr(firstTable.header_.columnNames_) << "\n";
					ss << "found    header: " << njh::conToStr(currentTable.header_.columnNames_) << '\n';
					throw std::runtime_error{ss.str()};
				}
				while(currentTable.getNextRow(row)){
					out << njh::conToStr(row, setUp.ioOptions_.outDelim_) << '\n';
				}
			}
		}
	}else{
		table mainTable;
		uint32_t count = 0;
		for (const auto &file : allFiles) {
			if (verbose) {
				std::cout << file.first.string() << std::endl;
			}
			if (njh::files::bfs::is_directory(file.first)) {
				if (verbose) {
					std::cout << "Skipping directory: " << file.first.string() << std::endl;
				}
				continue;
			}
			if (skipNonExistFiles && !bfs::exists(file.first)){
				continue;
			}
			if (0 == njh::files::bfs::file_size(file.first)) {
				if (verbose) {
					std::cout << "Skipping empty file: " << file.first.string() << std::endl;
				}
				continue;
			}
			if (count == 0) {
				table inTab(file.first.string(), setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_);
				mainTable = inTab;
			} else {
				table inTab(file.first.string(), setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_);
				try {
					mainTable.rbind(inTab, fill);
				}catch (std::exception & e) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", failed to add table from " << file.first << "\n";
					ss << e.what();
					throw std::runtime_error{ss.str()};
				}
			}
			++count;
		}
		mainTable.outPutContents(setUp.ioOptions_);
	}

	return 0;
}

int ManipulateTableRunner::cBind(
		const njh::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string contains = "";
	bool recursive = false;
	bool verbose = false;
	bool doingSelectColumns = false;
	std::string selColNamesStr = "";
	std::string selColPosStr = "";
	std::vector < uint32_t > selColPos;
	VecStr selColNames;
	if (setUp.setOption(selColPosStr,
			"-selPos,-selColPos,-selecteColumnsPos,-selColumnPos", "selectedColumns")
			|| setUp.setOption(selColNamesStr,
					"-selName,-selColName,-selecteColumnsName,-selColumnName,-selNames,-selColNames,-selecteColumnsNames,-selColumnNames",
					"selectedColumns")) {
		doingSelectColumns = true;
		if (selColPosStr != "") {
			selColPos = vecStrToVecNum<uint32_t>(tokenizeString(selColPosStr, ","));
		} else if (selColNamesStr != "") {
			selColNames = tokenizeString(selColNamesStr, ",");
		}
	}

	setUp.setOption(verbose, "-verbose,-v", "verbose");
	setUp.setUpRBind(contains, recursive);
	VecStr containsVec = tokenizeString(contains, ",");
	auto allFiles = getFiles(".", containsVec, "files", true, recursive);
	if (verbose) {
		std::cout << std::left;
		uint32_t firstMax = 0;
		uint32_t secondMax = 0;
		std::map<std::string, std::pair<std::string, bool>>::iterator fileIter;
		for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
			if (fileIter->first.length() > firstMax) {
				firstMax = fileIter->first.length();
			}
			if (fileIter->second.first.length() > secondMax) {
				secondMax = fileIter->first.length();
			}
		}
		for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
			std::cout << std::setw((int) firstMax) << std::right << fileIter->first
					<< " " << std::left << std::setw((int) secondMax)
					<< fileIter->second.first << std::endl;
		}
	}
	table inTab(setUp.ioOptions_);
	table mainTable;
	uint32_t count = 0;
	for (const auto &file : allFiles) {
		if (verbose)
			std::cout << count << std::endl;
		if (count == 0) {
			table inTab(file.first, setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_);
			if (doingSelectColumns) {
				if (selColNames.empty()) {
					mainTable = inTab.getColumns(selColPos);
				} else {
					mainTable = inTab.getColumns(selColNames);
				}
			} else {
				mainTable = inTab;
			}
		} else {
			table inTab(file.first, setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_);
			if (doingSelectColumns) {
				if (selColNames.empty()) {
					mainTable.cbind(inTab.getColumns(selColPos), false);
				} else {
					mainTable.cbind(inTab.getColumns(selColNames), false);
				}
			} else {
				mainTable.cbind(inTab, false);
			}
		}
		++count;
	}
	mainTable.outPutContents(setUp.ioOptions_);
	return 0;
}


int ManipulateTableRunner::printCol(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path fnp = "";
	std::string columnName = "";
	std::string delim = "\t";
	bool header = false;
	bool includeHeader = false;
	bool sort = false;
	bool unique = false;
	bool removeEmptyValues = false;
	bool removeNAs = false;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".txt";
	seqSetUp setUp(inputCommands);
	setUp.setOption(fnp, "--file,--fnp", "Filename path", true);
	setUp.setOption(columnName, "--columnName", "columnName", true);
	setUp.setOption(delim, "--delim", "delim");
	setUp.setOption(sort, "--sort", "sort output");
	setUp.setOption(unique, "--unique", "get unique values only");
	setUp.setOption(header, "--header", "Has header");
	setUp.setOption(removeEmptyValues, "--removeEmptyValues", "Remove Empty Values");
	setUp.setOption(removeNAs, "--removeNAs", "remove NAs");
	setUp.setOption(includeHeader, "--includeHeader", "Include Header");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	table tab(fnp.string(), delim, header);
	auto justColumn = tab.getColumns({columnName});
	OutputStream out(outOpts);
	if(sort){
		justColumn.sortTable(columnName, false);
	}
	if(unique){
		justColumn = justColumn.getUniqueRows();
	}
	auto col = justColumn.getColumn(columnName);

	if (includeHeader) {
		out << columnName << std::endl;
	}
	if(removeEmptyValues){
		removeElement(col, std::string(""));
	}
	if(removeNAs){
		//remove several possible NA values
		removeElements(col, VecStr{"NA", "na", "N/A", "Na"});
	}
	printVector(col, "\n", out);
	return 0;
}
}  // namespace njh
