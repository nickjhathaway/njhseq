#include "table.hpp"
#include "bibseq/IO/fileUtils.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {

void table::setColNamePositions(){
	colNameToPos_.clear();
	for(const auto & pos : iter::range(columnNames_.size())){
		colNameToPos_[columnNames_[pos]] = pos;
	}
}
uint32_t table::getColPos(const std::string & colName) const {
	auto search = colNameToPos_.find(colName);
	if (search == colNameToPos_.end()) {
		throw std::runtime_error { bib::bashCT::boldRed(
				"No column " + colName + " in header," + "available colnames are: "
						+ bib::bashCT::blue + bib::conToStr(columnNames_, ",")) };
	}
	return search->second;
}

VecStr table::getColumnLevels(uint32_t colPos){
	std::set<std::string> retSet;
	for(const auto & row : content_){
		retSet.emplace(row[colPos]);
	}
	return VecStr(retSet.begin(), retSet.end());
}

VecStr table::getColumnLevels(const std::string & colName){
	return getColumnLevels(getColPos(colName));
}

table::table(std::istream & in, const std::string &inDelim,
        bool header){
	VecStr columnNames;
	inDelim_ = inDelim;
	std::vector<VecStr> tempFileContent;
	if (inDelim == " " || inDelim == "whitespace") {
		std::string currentLine;
		std::getline(in, currentLine);
		if (header) {
			std::stringstream tempStream;
			tempStream << currentLine;
			VecStr tempVect;
			while (!tempStream.eof()) {
				std::string tempName;
				tempStream >> tempName;
				tempVect.emplace_back(tempName);
			}
			columnNames = tempVect;
			std::getline(in, currentLine);
		}
		while (!in.eof()) {
			std::stringstream tempStream;
			tempStream << currentLine;
			VecStr tempVect;
			while (!tempStream.eof()) {
				std::string tempName;
				tempStream >> tempName;
				tempVect.emplace_back(tempName);
			}
			std::getline(in, currentLine);
			tempFileContent.emplace_back(tempVect);
		}
	} else {
		if(inDelim == "tab"){
			inDelim_ = "\t";
		}
		std::string currentLine;
		std::getline(in, currentLine);
		if (header) {
			VecStr temp = tokenizeString(currentLine, inDelim_, true);
			columnNames = temp;
			std::getline(in, currentLine);
		}
		while (!in.eof()) {
			tempFileContent.emplace_back(tokenizeString(currentLine, inDelim_, true));
			std::getline(in, currentLine);
		}
	}
	if (header) {
		*this = table(tempFileContent, columnNames);
	} else {
		*this = table(tempFileContent);
	}
	setColNamePositions();
}
table::table(const std::string &filename, const std::string &inDelim,
             bool header) {
  std::ifstream textFile(filename.c_str());
  if (!textFile) {
  	std::stringstream ss;
    ss << bib::bashCT::red << bib::bashCT::bold
		<< "Error in opening " << filename
		<< bib::bashCT::reset << "\n";
    throw std::runtime_error{ss.str()};
  }
  *this = table(textFile, inDelim, header);
  setColNamePositions();
}
void table::addPaddingToEndOfRows() {
  size_t maxRowLength = 0;
  for (auto &iter : content_) {
    if (iter.size() > maxRowLength) {
      maxRowLength = iter.size();
    }
  }
  for (auto &iter : content_) {
    if (iter.size() < maxRowLength) {
      addOtherVec(iter, VecStr((int)maxRowLength - (int)iter.size(), ""));
    }
  }
}
void table::addPaddingZeros() {
  for (const auto &colPos : iter::range(columnNames_.size())) {
    // size_t colPos = getFirstPositionOfTarget(columnNames_, colIter);
    std::vector<size_t> blankPositions;
    VecStr nonBlanks;
    size_t pos = 0;
    for (const auto &fIter : content_) {
      if (fIter[colPos] == "") {
        blankPositions.emplace_back(pos);
      } else {
        nonBlanks.emplace_back(fIter[colPos]);
      }
      ++pos;
    }
    if (isVecOfDoubleStr(nonBlanks)) {
      for (const auto &i : blankPositions) {
        content_[i][colPos] = "0";
      }
    }
  }
  return;
}
void table::addColumn(const VecStr & columnNew, const std::string & name){
	//check size
	if(columnNew.size() != content_.size() && columnNew.size() != 1){
		std::cout << "new column's size doesn't match the size of the table not adding" << "\n";
		std::cout << "tableSize: " << content_.size() << "\n";
		std::cout << "addColumnSize: " << columnNew.size() << "\n";
	}else{
		columnNames_.emplace_back(name);
		if(columnNew.size() == 1){
			for(const auto & rowPos : iter::range(content_.size())){
				content_[rowPos].emplace_back(columnNew[0]);
			}
		}else{
			for(const auto & rowPos : iter::range(content_.size())){
				content_[rowPos].emplace_back(columnNew[rowPos]);
			}
		}
	}
	setColNamePositions();
	return;
}
void table::addZerosToEnds() {
  uint32_t maxLength = 0;
  for (const auto &row : content_) {
    if (row.size() > maxLength) {
      maxLength = row.size();
    }
  }
  for (auto &row : content_) {
    if (row.size() < maxLength) {
    	for (uint32_t i = 0; i < maxLength - row.size(); ++i){
      //for (const auto diff : iter::range(maxLength - row.size())) {
        row.emplace_back("0");
      }
    }
  }
}
table table::getColumns(const VecStr &specificColumnNames) const{
  std::vector<uint32_t> positions =
      getPositionsMultipleTargets(columnNames_, specificColumnNames);
  return getColumns(positions);
}

table table::getColumns(const std::vector<uint32_t> &specificColumnPositions)const {
  std::vector<VecStr> ans;
  for (const auto &fileIter : content_) {
    VecStr currentRow;
    ans.emplace_back(getTargetsAtPositions(fileIter, specificColumnPositions));
  }
  table outTab(ans, getTargetsAtPositions(columnNames_, specificColumnPositions));
  outTab.hasHeader_ = hasHeader_;
  return outTab;
}
table table::getColumnsNotAtPositions(
    const std::vector<uint32_t> &specificColumnPositions) const{
  std::vector<VecStr> ans;
  for (const auto &fileIter : content_) {
    VecStr currentRow;
    ans.emplace_back(getTargetsNotAtPositions(fileIter, specificColumnPositions));
  }
  return table(ans,
               getTargetsNotAtPositions(columnNames_, specificColumnPositions));
}

VecStr table::getColumn(const std::string &specifcColumnName) const {
  uint32_t colPos = getFirstPositionOfTarget(columnNames_, specifcColumnName);
  if(colPos == 4294967295){
  	std::stringstream ss;
  	ss << bib::bashCT::bold
				<< "Can't find col: " << bib::bashCT::red << specifcColumnName
				<< bib::bashCT::reset << "\n";
  	throw std::runtime_error{bib::bashCT::boldRed(ss.str())};
  }
  return getColumn(colPos);
}
VecStr table::getColumn(uint32_t colPos) const {
  VecStr ans;
  if (colPos >= columnNames_.size()) {
  	std::stringstream ss;
  	ss << "positions: " << colPos << " is out of the bounds of "
              << columnNames_.size() << " return nothing" << "\n";
  	throw std::runtime_error{bib::bashCT::boldRed(ss.str())};
  } else {
    for (const auto &fIter : content_) {
      ans.emplace_back(fIter[colPos]);
    }
  }
  return ans;
}

std::vector<std::string *> table::getColumnPointer(
    const std::string &specifcColumnName) {
  size_t colPos = getFirstPositionOfTarget(columnNames_, specifcColumnName);
  std::vector<std::string *> ans;
  for (auto &fIter : content_) {
    ans.emplace_back(&fIter[colPos]);
  }
  return ans;
}

void table::deleteColumn(const std::string &columnName) {
  size_t pos = getFirstPositionOfTarget(columnNames_, columnName);
  deleteColumn(pos);
}
void table::deleteColumn(size_t columnIndex) {
  if (columnIndex < columnNames_.size()) {
    columnNames_.erase(columnNames_.begin() + columnIndex);
    for (auto &rIter : content_) {
      rIter.erase(rIter.begin() + columnIndex);
    }
  } else {
    std::cout << "Column " << columnIndex << " doesn't exist, table left alone"
              << "\n";
  }
}
void table::deleteRow(size_t rowIndex) {
  if (rowIndex < content_.size()) {
    content_.erase(content_.begin() + rowIndex);
  } else {
    std::cout << "Row " << rowIndex << " is out of bounds, table left alone"
              << "\n";
  }
}

table table::getRows(const std::string &forColumn,
                     const std::string &element) const {
  // std::vector<VecStr> ans;
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions = getPositionsOfTarget(col, element);

  return table(getTargetsAtPositions(content_, positions), columnNames_);
}
table table::getRowsLoose(const std::string &forColumn,
                          const std::string &subString) const {
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions = getPositionsOfSubStrTarget(col, subString);

  return table(getTargetsAtPositions(content_, positions), columnNames_);
}
table table::getRowsStartsWith(const std::string &forColumn,
                               const std::string &startsWtih) const{
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions =
      getPositionsOfTargetStartsWith(col, startsWtih);

  return table(getTargetsAtPositions(content_, positions), columnNames_);
}

table table::getRows(const std::vector<uint32_t> &specificRowPositions) const{

  return table(getTargetsAtPositions(content_, specificRowPositions),
               columnNames_);
}
void table::outPutContents(std::ostream &out, std::string delim) const {
  if (delim == "tab") {
    delim = "\t";
  } else if (delim == "whitespace") {
    delim = " ";
  }
  if (hasHeader_) {
    out << vectorToString(columnNames_, delim) << "\n";
  }
  outputVectorOfVectors(content_, delim, out);
}

void table::sortTable(const std::string &byThisColumn, bool decending,
                      bool guessFormat) {
  if (!vectorContains(columnNames_, byThisColumn)) {
    std::cout << "Table does not contain " << byThisColumn
              << " not sorting table" << "\n";
    std::cout << "options are" << "\n";
    printVector(columnNames_, ",");
    return;
  }
  VecStr col = getColumn(byThisColumn);
  if (guessFormat && isVecOfDoubleStr(col)) {
    if (isVecOfIntStr(col)) {
      std::multimap<int32_t, VecStr> mapSorter;
      size_t count = 0;
      for (const auto &iter : content_) {
        mapSorter.insert(std::make_pair(std::stoi(col[count]), iter));
        ++count;
      }
      content_.clear();
      for (const auto &iter : mapSorter) {
        content_.emplace_back(iter.second);
      }
    } else {
      std::multimap<double, VecStr> mapSorter;
      size_t count = 0;
      for (const auto &iter : content_) {
        mapSorter.insert(std::make_pair(std::stod(col[count]), iter));
        ++count;
      }
      content_.clear();
      for (const auto &iter : mapSorter) {
        content_.emplace_back(iter.second);
      }
    }
  } else {
    std::multimap<std::string, VecStr> mapSorter;
    size_t count = 0;
    for (std::vector<VecStr>::iterator iter = content_.begin();
         iter != content_.end(); ++iter) {
      mapSorter.insert(std::make_pair(col[count], *iter));
      ++count;
    }
    content_.clear();
    for (std::multimap<std::string, VecStr>::iterator iter = mapSorter.begin();
         iter != mapSorter.end(); ++iter) {
      content_.emplace_back(iter->second);
    }
  }
  if (decending) {
    reverseRows();
  }
  return;
}

std::map<std::string, table> table::splitTableOnColumn(
    const std::string &colName) const {

  VecStr factorNames = getUniqueStrings(getColumn(colName));
  //std::cout << factorNames << "\n";
  std::map<std::string, table> ans;
  for (const auto &iter : factorNames) {
    ans.insert({iter, getRows(colName, iter)});
  }
  return ans;
}
std::map<std::string, table> table::splitTableOnColumnLoose(
    const std::string &colName) const {
  VecStr factorNames = getUniqueStrings(getColumn(colName));
  return splitTableOnColumnLoose(colName, factorNames);
}
std::map<std::string, table> table::splitTableOnColumnLoose(
    const std::string &colName, const VecStr &useTheseNames) const {
  std::map<std::string, table> ans;
  for (const auto &iter : useTheseNames) {
    ans.insert({iter, getRowsLoose(colName, iter)});
  }
  return ans;
}

std::map<std::string, table> table::splitTableOnColumnLoose(
    const std::string &colName, const std::string &firstOccurnceOf) const {
  auto names = getColumn(colName);
  trimStringsAtFirstOccurence(names, firstOccurnceOf);
  names = getUniqueStrings(names);
  return splitTableOnColumnLoose(colName, names);
}

void table::rbind(const table &otherTable) {
  // should add in column name checking
  addOtherVec(content_, otherTable.content_);
  addPaddingToEndOfRows();
}
void table::cbind(const table &otherTable) {
  size_t count = 0;
  size_t maxRows = content_.size() - 1;
  size_t numCol = columnNames_.size();
  for (const auto &tabIter : otherTable.content_) {
    if (count > maxRows) {
      VecStr empty(numCol, "");
      addOtherVec(empty, tabIter);
      content_.emplace_back(empty);
    } else {
      addOtherVec(content_[count], tabIter);
    }
    ++count;
  }
  addOtherVec(columnNames_, otherTable.columnNames_);
}

table table::getColumnsLoose(const std::string &subStr) const{
  std::vector<uint32_t> positions =
      getPositionsOfSubStrTarget(columnNames_, subStr);
  if (positions.size() == 0) {
    std::cout << "Couldn't find any columns that contained: " << subStr
              << " returning the whole table" << "\n";
    return *this;
  } else {
    return getColumns(positions);
  }
}

table table::getColumnsStartWith(const std::string &startsWith) const{
  std::vector<uint32_t> positions =
      getPositionsOfTargetStartsWith(columnNames_, startsWith);
  if (positions.size() == 0) {
    std::cout << "Couldn't find any columns that started with: " << startsWith
              << " returning the whole table" << "\n";
    return *this;
  } else {
    return getColumns(positions);
  }
}

table table::getUniqueRows() const{
  return table(collapseUniqueVectors(content_), columnNames_);
}

void table::reverseColumns() {
  for (auto &rowIter : content_) {
    std::reverse(rowIter.begin(), rowIter.end());
  }
  std::reverse(columnNames_.begin(), columnNames_.end());
}
void table::reverseRows() { std::reverse(content_.begin(), content_.end()); }

void table::outPutContentOrganized(std::ostream &out) const {
  if (hasHeader_) {
    printTableOrganized(content_, columnNames_, out);
  } else {
    printTableOrganized(content_, out);
  }
  return;
}

void table::printOutSplitTable(const std::map<std::string, table> &tabSplit,
                               std::ostream &out, const std::string &delim,
                               bool outOrganized) {
  for (const auto &iter : tabSplit) {
    out << iter.first << "\n";
    if (outOrganized) {
      iter.second.outPutContentOrganized(out);
    } else {
      iter.second.outPutContents(out, delim);
    }
  }
  return;
}

void table::trimElementsAtFirstOccurenceOf(const std::string &trimAt) {
  for (auto &row : content_) {
    trimStringsAtFirstOccurence(row, trimAt);
  }
  return;
}

table table::getStatsTable() const {
  table outTable;
  outTable.hasHeader_ = true;
  outTable.columnNames_.emplace_back("stat");
  for (auto sIter : {"sum", "mean", "median", "min", "max", "std"}) {
    outTable.content_.emplace_back(VecStr{sIter});
  }
  for (auto &colNameIter : columnNames_) {
    VecStr currentColumn = getColumn(colNameIter);
    if (isVecOfDoubleStr(currentColumn)) {
      outTable.columnNames_.emplace_back(colNameIter);
      std::vector<double> currentNumbers;
      for (const auto &numStr : currentColumn) {
        currentNumbers.emplace_back(std::stod(numStr));
      }
      auto currentStats = getStatsOnVec(currentNumbers);
      for (auto i : iter::range(0, (int)outTable.content_.size())) {
        if (i == 0) {
          outTable.content_[i].emplace_back(estd::to_string(currentStats["sum"]));
        } else if (i == 1) {
          outTable.content_[i].emplace_back(estd::to_string(currentStats["mean"]));
        } else if (i == 2) {
          outTable.content_[i]
              .emplace_back(estd::to_string(currentStats["median"]));
        } else if (i == 3) {
          outTable.content_[i].emplace_back(estd::to_string(currentStats["min"]));
        } else if (i == 4) {
          outTable.content_[i].emplace_back(estd::to_string(currentStats["max"]));
        } else if (i == 5) {
          outTable.content_[i].emplace_back(estd::to_string(currentStats["std"]));
        }
      }
    }
  }
  return outTable;
}

table table::aggregateSimple(const std::string &columnName,
                             const std::string &function, bool addZeros) {
  if (function == "mean" || function == "median" || function == "min" ||
      function == "max" || function == "sum" || function == "std") {
    if (addZeros) {
      addPaddingZeros();
    }
    auto tables = splitTableOnColumn(columnName);
    std::map<std::string, table> stats;
    for (const auto &tabIter : tables) {
      stats.insert({tabIter.first, tabIter.second.getStatsTable()});
    }
    table outTable;
    outTable.hasHeader_ = true;
    outTable.columnNames_.emplace_back(columnName);
    VecStr otherColNames = stats.begin()->second.columnNames_;
    std::string remove = "stat";
    removeElement(otherColNames, remove);
    addOtherVec(outTable.columnNames_, otherColNames);
    for (const auto &statsIter : stats) {
      VecStr currentLine;
      currentLine.emplace_back(statsIter.first);
      table info = statsIter.second.getRows("stat", function);
      VecStr currentInfo = info.content_.front();
      removeElement(currentInfo, function);
      addOtherVec(currentLine, currentInfo);
      outTable.content_.emplace_back(currentLine);
    }
    return outTable;
  } else {
    std::cout << "unrecognized function for aggregate," << function
              << "\n";
    std::cout << "Available optins are sum, mean, median, max, min, std"
              << "\n";
    return *this;
  }
}

void table::outPutContents(outOptions options) const {
  if (options.outDelim_ == "tab") {
    options.outDelim_ = "\t";
  } else if (options.outDelim_ == "whitespace") {
    options.outDelim_ = " ";
  }
  if (options.outFilename_ == "") {
    if (options.outOrganized_) {
      outPutContentOrganized(std::cout);
    } else {
      outPutContents(std::cout, options.outDelim_);
    }
  } else {
    std::ofstream outFile;
    openTextFile(outFile, options.outFilename_, options.extention_,
                 options.overWriteFile_, options.exitOnFailureToWrite_);
    if (options.outOrganized_) {
      outPutContentOrganized(outFile);
    } else {
      outPutContents(outFile, options.outDelim_);
    }
  }
}

std::vector<uint32_t> table::getNumericColumnPositions(bool addZeros) {
  if (addZeros) {
    addPaddingZeros();
  }
  std::vector<uint32_t> ans;
  uint32_t pos = 0;
  for (const auto &col : columnNames_) {
    VecStr currentColumn = getColumn(col);
    if (isVecOfDoubleStr(currentColumn)) {
      ans.emplace_back(pos);
    }
    ++pos;
  }
  return ans;
}
table table::getAllNumericColumns(bool addZeros) {
  std::vector<uint32_t> numericPosisitions =
      getNumericColumnPositions(addZeros);
  return getColumns(numericPosisitions);
}
table table::getAllNonNumericColumns(bool addZeros) {
  std::vector<uint32_t> numericPosisitions =
      getNumericColumnPositions(addZeros);
  return getColumnsNotAtPositions(numericPosisitions);
}
table table::aggregateAdvance(const std::string &columnName,
                              const std::string &function, bool addZeros) {
  std::string sep = "______";
  table numericColumns = getAllNumericColumns();
  table nonNumericColumns = getAllNonNumericColumns();
  std::vector<VecStr> collapsedColumns;
  for (const auto row : nonNumericColumns.content_) {
    collapsedColumns.emplace_back(VecStr{vectorToString(row, sep)});
  }

  table tempTable(collapsedColumns, {"collapsed"});
  tempTable.cbind(numericColumns);
  table aggregated = tempTable.aggregateSimple("collapsed", function);
  VecStr collapsedColumn = aggregated.getColumn("collapsed");
  aggregated.deleteColumn("collapsed");
  std::vector<VecStr> expanded;
  for (const auto &element : collapsedColumn) {
    expanded.emplace_back(tokenizeString(element, sep));
  }
  table expandedTable(expanded, nonNumericColumns.columnNames_);
  expandedTable.cbind(aggregated);
  return expandedTable;
}

table table::cbind(const std::vector<table> &tables,
                   const std::string &columnForceMatch, bool addZeros) {
  table ans;
  std::map<int, std::map<std::string, table>> splitTables;
  int count = 0;
  VecStr names;
  uint64_t numOfCols = 0;
  VecStr spNames;
  for (const auto &tab : tables) {
    splitTables[count] = tab.splitTableOnColumn(columnForceMatch);
    for (auto &secondTab : splitTables[count]) {
      secondTab.second.deleteColumn(columnForceMatch);
      numOfCols = secondTab.second.columnNames_.size();
      spNames = secondTab.second.columnNames_;
    }
    addOtherVec(names, getVectorOfMapKeys(splitTables[count]));
    ++count;
  }
  names = getUniqueStrings(names);
  std::sort(names.begin(), names.end());
  ans.columnNames_.emplace_back(columnForceMatch);
  for (const auto &name : names) {
    ans.content_.emplace_back(VecStr{name});
  }

  addOtherVec(ans.columnNames_, repeatVector(spNames, {splitTables.size()}));
  for (auto &row : ans.content_) {
    for (auto &spTab : splitTables) {
      if (spTab.second.find(row[0]) != spTab.second.end()) {
        addOtherVec(row, spTab.second[row[0]].content_[0]);
      } else {
        addOtherVec(row, VecStr(numOfCols, ""));
      }
    }
  }

  if (addZeros) {
    ans.addPaddingZeros();
  }
  return ans;
}
std::map<std::string, uint32_t> table::countColumn(
    const std::string &columnName) {

  auto colPos = getPositionsOfTarget(columnNames_, columnName);
  if (colPos.empty()) {
  	std::stringstream ss;
    ss << "No column named, "
    		<< bib::bashCT::bold << columnName
    		<< bib::bashCT::reset << "\n";
    ss << "Options are : "; printVector(columnNames_, ", ");
    throw std::runtime_error{ss.str()};
  }
  return countColumn(colPos.front());
}
std::map<std::string, uint32_t> table::countColumn(uint32_t colPos) {
  std::map<std::string, uint32_t> ans;
  if (colPos > columnNames_.size()) {
  	std::stringstream ss;
    ss << "colPos: " << colPos << ", out of bounds of "
              << columnNames_.size();
    throw std::runtime_error{ss.str()};
  }
  for (const auto &row : content_) {
    ++ans[row[colPos]];
  }
  return ans;
}
table table::cbind(const std::map<std::string, table> &tables,
                   const std::string &columnForceMatch, bool addZeros) {
  std::vector<table> tableVector = getVectorOfMapValues(tables);
  return cbind(tableVector, columnForceMatch, addZeros);
}
void table::removeEmpty(bool addPadding) {
	std::vector<uint32_t> emptyPositions;
	for(const auto & rowPos : iter::range(content_.size())){
		if(content_.empty()){
			emptyPositions.emplace_back(rowPos);
		}else {
			bool empty =true;
			for(const auto & colPos : iter::range(content_[rowPos].size())){
				if(!allWhiteSpaceStr(content_[rowPos][colPos])){
					empty = false;
					break;
				}
			}
			if(empty){
				emptyPositions.emplace_back(rowPos);
			}
		}
	}
	if(!emptyPositions.empty()){
		bib::sort(emptyPositions);
		for(const auto & pos : iter::reverse(emptyPositions)){
			content_.erase(content_.begin() + pos);
		}
	}

}


table table::countColumn(const VecStr &columnNames){
	std::vector<uint32_t> colPositions;
	for(const auto & col :columnNames){
		colPositions.emplace_back(getColPos(col));
	}
	return countColumn(colPositions);
}
table table::countColumn(const std::vector<uint32_t> & colPositions){
	bool pass = true;
	for(const auto & pos : colPositions){
		if(pos >= columnNames_.size()){
			pass = false;
			break;
		}
	}
	if(!pass){
		std::stringstream ss;
		ss << bib::bashCT::red << bib::bashCT::bold
				<< "Error in table::countColumn(const std::vector<uint32_t> & colPositions), out of range error\n"
				<< bib::bashCT::blue << vectorToString(colPositions, ",")
				<< bib::bashCT::red << " out of range of " << columnNames_.size()
				<< bib::bashCT::reset << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::map<std::string, uint32_t> counts;
	for(const auto & row : content_){
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
		outColumnNames.emplace_back(columnNames_[pos]);
	}
	outColumnNames.emplace_back("count");
	table ret(outColumnNames);
	for(const auto & c : counts){
		auto toks = tokenizeString(c.first,"SPLITONTHIS");
		toks.emplace_back(estd::to_string(c.second));
		ret.content_.emplace_back(toks);
	}
	return ret;
}

table table::extractNumColGreater(const std::string & colName, double cutOff)const{
  auto comp = [&](const std::string & str){
  	double numValue = std::stod(str);
  	return numValue > cutOff;
  };
  return extractByComp(colName, comp);
}
table table::extractNumColGreater(uint32_t colPos, double cutOff)const{
  auto comp = [&](const std::string & str){
  	double numValue = std::stod(str);
  	return numValue > cutOff;
  };
  return extractByComp(colPos, comp);
}

}  // namespace bib
