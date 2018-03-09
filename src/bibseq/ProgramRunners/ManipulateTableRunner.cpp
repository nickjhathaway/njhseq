

#include "ManipulateTableRunner.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"
#include "bibseq/objects/Meta/MetaDataInName.hpp"
#include "bibseq/programUtils/seqSetUp.hpp"
#include "bibseq/IO.h"

namespace bibseq {

ManipulateTableRunner::ManipulateTableRunner() :
		bib::progutils::ProgramRunner(
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
				}, "ManipulateTable", "1") {
}
//
int ManipulateTableRunner::splitColumnContainingMeta(
		const bib::progutils::CmdArgs & inputCommands) {

	std::string column = "";
	bool keepMetaInColumn = false;
	bool removeEmptyColumn = false;
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setOption(column, "--column", "Column to split which contains meta with the following formating [key1=val1;key2=val2;]", true);
	setUp.setOption(keepMetaInColumn, "--keepMetaInColumn", "Keep meta in the original column");
	setUp.setOption(removeEmptyColumn, "--removeEmptyColumn", "remove original column if it becomes an empty column");
	setUp.processFileName(true);
	setUp.processNonRquiredDefaults();
	bool sorting = setUp.processSorting();
	setUp.finishSetUp(std::cout);
	table inTab(setUp.ioOptions_);
	inTab.checkForColumnsThrow({column}, __PRETTY_FUNCTION__);
	std::set<std::string> metaFields;
	std::unordered_map<std::string, VecStr> metaValues;

	bool noneContainMeta = false;
	for (const auto & row : inTab.content_) {
		if(MetaDataInName::nameHasMetaData(row[inTab.getColPos(column)])){
			MetaDataInName rowMeta(row[inTab.getColPos(column)]);
			auto metas = getVectorOfMapKeys(rowMeta.meta_);
			std::copy(metas.begin(), metas.end(),
					std::inserter(metaFields, metaFields.end()));
			noneContainMeta = false;
		}
	}
	bool allEmpty = true;
	for (auto & row : inTab.content_) {
		if(MetaDataInName::nameHasMetaData(row[inTab.getColPos(column)])){
			MetaDataInName rowMeta(row[inTab.getColPos(column)]);
			for(const auto & m : rowMeta.meta_){
				metaValues[m.first].emplace_back(m.second);
			}
			for(const auto & metaField : metaFields){
				if(!bib::in(metaField, rowMeta.meta_)){
					metaValues[metaField].emplace_back("NA");
				}
			}
			if(!keepMetaInColumn){
				MetaDataInName::removeMetaDataInName(row[inTab.getColPos(column)]);
				if("" != row[inTab.getColPos(column)]){
					allEmpty = false;
				}
			}
		}else{
			for(const auto & metaField : metaFields){
				metaValues[metaField].emplace_back("NA");
			}
		}
	}
	if(allEmpty && removeEmptyColumn){
		inTab.deleteColumn(column);
	}
	for(const auto & m : metaValues){
		inTab.addColumn(m.second, m.first);
	}
	if(sorting){
		inTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	inTab.hasHeader_ = setUp.addHeader_;
	setUp.ioOptions_.hasHeader_ = setUp.addHeader_;
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}

//
int ManipulateTableRunner::extractColumnElementLength(
		const bib::progutils::CmdArgs & inputCommands) {
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
		const bib::progutils::CmdArgs & inputCommands) {
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
	while (bib::files::crossPlatGetline(inFile, line)) {
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
		const bib::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string columnName = "";
	double cutOff = 1;
	bool lessThan = false;
	setUp.processDefaultProgram(true);
	setUp.setOption(columnName, "--columnName",
			"Name of the column to extract on using cutOff", true);
	setUp.setOption(cutOff, "--cutOff",
			"Cut off of column to be extract, not inclusive, can be negative");
	setUp.setOption(lessThan, "--lessThan",
			"Take numbers less than value in cutOff flag");
	setUp.finishSetUp(std::cout);
	table inTab(setUp.ioOptions_);
	table outTab;
	if (lessThan) {
		auto compLess = [cutOff](const std::string & str) {
			double numValue = std::stod(str);
			return numValue < cutOff;
		};
		outTab = inTab.extractByComp(columnName, compLess);
	} else {
		auto compGreater = [cutOff](const std::string & str) {
			double numValue = std::stod(str);
			return numValue > cutOff;};
		outTab = inTab.extractByComp(columnName, compGreater);
	}
	outTab.hasHeader_ = setUp.ioOptions_.hasHeader_;
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::addColumn(
		const bib::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string columnName = "";
	VecStr element;
	std::string elementStr = "";
	setUp.processDefaultProgram(true);
	setUp.setOption(columnName, "-newColumnName",
			"Name of the new Column to add to table", true);
	setUp.setOption(elementStr, "-element",
			"What to Add to the Table Under Column, can be several comma sep values or just one",
			true);
	setUp.finishSetUp(std::cout);
	table outTab(setUp.ioOptions_);
	auto toks = tokenizeString(elementStr, ",");
	if (toks.size() > 1) {
		if (outTab.content_.size() % toks.size() != 0) {
			std::stringstream ss;
			ss << bib::bashCT::red << "Error, table has "
					<< outTab.content_.size() << " and the size of adding elements "
					<< bib::conToStr(elementStr, ",") << " doesn't fit into it"
					<< bib::bashCT::reset << std::endl;
			throw std::runtime_error{ss.str()};
		} else {
			outTab.addColumn(
					repeatVector(toks, { outTab.content_.size() / toks.size() }),
					columnName);
		}
	} else if (toks.size() == 1) {
		outTab.addColumn(toks, columnName);
	} else {
		//shouldn't be zero...
	}
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}



int ManipulateTableRunner::roughHistogramOfColumn(
		const bib::progutils::CmdArgs & inputCommands) {
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
		if(row[colPos] != "*" && bib::strToLowerRet(row[colPos]) != "na" && bib::strToLowerRet(row[colPos]) != "nan"){
			columnValues.emplace_back(bib::StrToNumConverter::stoToNum<double>(row[colPos]));
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
			return bib::pasteAsStr(min_, "-", max_);
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
		const bib::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.processDefaultProgram(true);
	std::string columnName = "";
	setUp.setOption(columnName, "--columnName", "columnName", true);
	setUp.finishSetUp(std::cout);
	auto toks = tokenizeString(columnName, ",");
	table inTab(setUp.ioOptions_);
	if (toks.size() == 1) {
		auto counts = inTab.countColumn(columnName);
		table ret(counts, { "element", "count" });
		if (setUp.sortByColumn_ != "") {
			ret.sortTable(setUp.sortByColumn_, setUp.decending_);
		}
		ret.outPutContents(setUp.ioOptions_);
	} else {
		auto ret = inTab.countColumn(toks);
		if (setUp.sortByColumn_ != "") {
			ret.sortTable(setUp.sortByColumn_, setUp.decending_);
		}
		ret.outPutContents(setUp.ioOptions_);
	}

	return 0;
}

int ManipulateTableRunner::catOrganized(
		const bib::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setUpCatOrganized();
	table inTab(setUp.ioOptions_);
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::changeDelim(
		const bib::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setUpChangeDelim();
	table inTab(setUp.ioOptions_);
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::sortTable(
		const bib::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	setUp.setUpSortTable();
	table inTab(setUp.ioOptions_);
	inTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	inTab.outPutContents(setUp.ioOptions_);
	return 0;
}




int ManipulateTableRunner::tableExtractColumns(const bib::progutils::CmdArgs & inputCommands){
	ManipulateTableSetUp setUp(inputCommands);
	std::string columns;
	bool getUniqueRows = false;
	setUp.processFileName();
	setUp.processNonRquiredDefaults();
	setUp.processSorting();
  std::string extractColumnsStrings = "";
  setUp.setOption(columns, "--columns",
                  "Names Of Columns To Extract, either comma separated input or a file with each line a column name",
									true);
  setUp.setOption(getUniqueRows, "-getUniqueRows", "GetUniqueRows");
  setUp.finishSetUp(std::cout);

  auto extractColumns = getInputValues(columns, ",");

	table inTab(setUp.ioOptions_);
	table outTab = inTab.getColumns(extractColumns);

	if (setUp.sortByColumn_ != "") {

		outTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}

	if (getUniqueRows) {
		outTab = outTab.getUniqueRows();
	}

	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::tableExtractElementsWithPattern(const bib::progutils::CmdArgs & inputCommands){
	std::string column;
	std::string patStr = "";
	bool getUniqueRows = false;
	ManipulateTableSetUp setUp(inputCommands);

	setUp.processFileName();
	setUp.processNonRquiredDefaults();
	setUp.processSorting();
	setUp.setOption(patStr, "--patStr", "Pattern to match elements in column with", true);
	setUp.setOption(column, "--column", "Name of the column to search", true);
	setUp.setOption(getUniqueRows, "--getUniqueRows", "GetUniqueRows");
	setUp.finishSetUp(std::cout);

	table inTab(setUp.ioOptions_);
	table outTab = inTab.getRowsMatchingPattern(column, std::regex{patStr});

	if (setUp.sortByColumn_ != "") {
		outTab.sortTable(setUp.sortByColumn_, setUp.decending_);
	}
	if (getUniqueRows) {
		outTab = inTab.getUniqueRows();
	}
	outTab.outPutContents(setUp.ioOptions_);
	return 0;
}

int ManipulateTableRunner::tableExtractElementsStartingWith(const bib::progutils::CmdArgs & inputCommands){
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

int ManipulateTableRunner::tableExtractColumnsStartsWith(const bib::progutils::CmdArgs & inputCommands){
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

int ManipulateTableRunner::tableExtractColumnsWithPattern(const bib::progutils::CmdArgs & inputCommands){
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
		const bib::progutils::CmdArgs & inputCommands) {
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
		const bib::progutils::CmdArgs & inputCommands) {
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
		const bib::progutils::CmdArgs & inputCommands) {
	ManipulateTableSetUp setUp(inputCommands);
	std::string column = "";
	bool getStatsInstead = false;
	bool splitingLoose = false;
	std::string splitLooseOccurrence = "";

	setUp.setUpSplitTable(column, getStatsInstead, splitingLoose,
			splitLooseOccurrence);
	table inTab(setUp.ioOptions_);
	inTab.fillWithZeros();
	std::map < std::string, table > tables;
	if (splitingLoose) {
		if (splitLooseOccurrence == "") {
			tables = inTab.splitTableOnColumnLoose(column);
		} else {
			tables = inTab.splitTableOnColumnLoose(column, splitLooseOccurrence);
		}
	} else {
		tables = inTab.splitTableOnColumn(column);
	}
	std::map < std::string, table > outTables;
	if (getStatsInstead) {
		for (const auto &tabIter : tables) {
			outTables.insert( { tabIter.first, tabIter.second.getStatsTable() });
		}
	} else {
		outTables = tables;
	}
	if (setUp.directoryName_ == "./") {
		table::printOutSplitTable(outTables, std::cout,
				setUp.ioOptions_.outDelim_, setUp.ioOptions_.outOrganized_);
	} else {
		for (const auto &tab : outTables) {
			std::ofstream outFile;
			openTextFile(outFile,
					bib::files::make_path(setUp.directoryName_, tab.first).string(),
					".txt", true, false);
			if (setUp.ioOptions_.outOrganized_) {
				tab.second.outPutContentOrganized(outFile);
			} else {
				tab.second.outPutContents(outFile, setUp.ioOptions_.outDelim_);
			}
		}
	}
	return 0;
}

int ManipulateTableRunner::aggregateTable(
		const bib::progutils::CmdArgs & inputCommands) {
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
		const bib::progutils::CmdArgs & inputCommands) {
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
inline std::map<bib::files::bfs::path, bool> listAllFiles(
		const std::string & dirName, bool recursive,
		const std::vector<std::string>& contains,
		const std::vector<std::string>& doesNotContain,
		uint32_t levels = std::numeric_limits < uint32_t > ::max()) {
	std::map < bfs::path, bfs::path > filesGathering;
	bib::files::listAllFilesHelper(dirName, recursive, filesGathering, 1, levels);
	std::map<bfs::path, bool> files = bib::files::convertMapFnpFnpToFnpIsDir(
			filesGathering);

	if (!contains.empty()) {
		std::map<bib::files::bfs::path, bool> specificFiles;
		for (const auto & f : files) {
			if (bib::checkForSubStrs(f.first.string(), contains)
					&& !bib::checkForSubStrs(f.first.string(), doesNotContain)) {
				specificFiles.emplace(f);
			}
		}
		return specificFiles;
	}
	return files;
}

int ManipulateTableRunner::rBind(
		const bib::progutils::CmdArgs & inputCommands) {
	std::string contains = "";
	std::string doesNotContains = "";
	bool recursive = false;
	bool verbose = false;
	std::string files = "";
	uint32_t levels = std::numeric_limits < uint32_t > ::max();

	ManipulateTableSetUp setUp(inputCommands);
	setUp.description_ = "Program to bind together text files representing tables with similar information";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --files file1.tsv,file2.tsv,file3.tsv --header --delim tab --out combined.tsv --overWrite #combined a list of files");
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --contains _info.txt --header --delim tab --out combined.tsv --overWrite #combine together all files in current directory that contain _info.txt in their name");
	setUp.setOption(levels, "--depth", "Depth of search");
	setUp.setOption(verbose, "--verbose,-v", "verbose");
	setUp.setOption(doesNotContains, "--doesNotContains", "Does Not Contains");
	setUp.setOption(files, "--files", "A direct list of files to combine");
	setUp.processDefaultProgram(false);
	setUp.setOption(contains, "--contains", "contains", "" == files);
	setUp.setOption(recursive, "--recursive", "recursive");
	if(setUp.needsHelp()){
		setUp.printFlags(std::cout);
		exit(1);
	}
	setUp.finishSetUp(std::cout);
	setUp.ioOptions_.out_.throwIfOutExistsNoOverWrite(__PRETTY_FUNCTION__);
	VecStr containsVec = tokenizeString(contains, ",");
	VecStr doesNotContainsVec = tokenizeString(doesNotContains, ",");
	std::map<bib::files::bfs::path, bool> allFiles;
	if("" != files){
		auto vals = getInputValues(files, ",");
		for(const auto & f : vals){
			allFiles.emplace(f, false);
		}
	}else{
		if (doesNotContains != "") {
			allFiles = listAllFiles("./", recursive, containsVec, doesNotContainsVec,
					levels);
		} else {
			allFiles = bib::files::listAllFiles("./", recursive, containsVec, levels);
		}
	}


	table mainTable;
	uint32_t count = 0;
	for (const auto &file : allFiles) {
		if (verbose) {
			std::cout << file.first.string() << std::endl;
		}
		if (bib::files::bfs::is_directory(file.first)) {
			if (verbose) {
				std::cout << "Skipping directory: " << file.first.string() << std::endl;
			}
			continue;
		}
		if (count == 0) {
			table inTab(file.first.string(), setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_);
			mainTable = inTab;
		} else {
			table inTab(file.first.string(), setUp.ioOptions_.inDelim_, setUp.ioOptions_.hasHeader_);
			try {
				mainTable.rbind(inTab, false);
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
	return 0;
}

int ManipulateTableRunner::cBind(
		const bib::progutils::CmdArgs & inputCommands) {
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
		const bib::progutils::CmdArgs & inputCommands) {
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
	setUp.setOption(fnp, "--fnp", "Filename path", true);
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
}  // namespace bib
