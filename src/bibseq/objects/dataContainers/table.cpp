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

#include "table.hpp"
#include "bibseq/IO/fileUtils.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {
table::table(std::istream & in, const std::string &inDelim,
        bool header){
	VecStr columnNames;
	inDelim_ = inDelim;
	std::vector<VecStr> tempFileContent;
	if (inDelim == " " || inDelim == "whitespace") {
		std::string currentLine;
		/*std::fstream textFile(filename.c_str());
		if (!textFile) {
			std::cout << "Error in opening " << filename << std::endl;
			exit(1);
		}*/
		getline(in, currentLine);
		if (header) {
			std::stringstream tempStream;
			tempStream << currentLine;
			VecStr tempVect;
			while (!tempStream.eof()) {
				std::string tempName;
				tempStream >> tempName;
				tempVect.push_back(tempName);
			}
			columnNames = tempVect;
			getline(in, currentLine);
		}
		while (!in.eof()) {
			std::stringstream tempStream;
			tempStream << currentLine;
			VecStr tempVect;
			while (!tempStream.eof()) {
				std::string tempName;
				tempStream >> tempName;
				tempVect.push_back(tempName);
			}
			getline(in, currentLine);
			tempFileContent.push_back(tempVect);
		}
	} else {
		std::string currentLine;
		/*std::ifstream textFile(filename.c_str());
		if (!textFile) {
			std::cout << "Error in opening " << filename << std::endl;
			exit(1);
		}*/
		getline(in, currentLine);
		if (header) {
			VecStr temp = tokenizeString(currentLine, inDelim, true);
			columnNames = temp;
			getline(in, currentLine);
		}
		while (!in.eof()) {
			VecStr temp = tokenizeString(currentLine, inDelim, true);
			tempFileContent.push_back(tokenizeString(currentLine, inDelim));
			getline(in, currentLine);
		}
	}
	if (header) {
		*this = table(tempFileContent, columnNames);
	} else {
		*this = table(tempFileContent);
	}
}
table::table(const std::string &filename, const std::string &inDelim,
             bool header) {
  VecStr columnNames;
  inDelim_ = inDelim;
  std::vector<VecStr> tempFileContent;
  if (inDelim == " " || inDelim == "whitespace") {
    std::string currentLine;
    std::fstream textFile(filename.c_str());
    if (!textFile) {
      std::cout << "Error in opening " << filename << std::endl;
      exit(1);
    }
    getline(textFile, currentLine);
    if (header) {
      std::stringstream tempStream;
      tempStream << currentLine;
      VecStr tempVect;
      while (!tempStream.eof()) {
        std::string tempName;
        tempStream >> tempName;
        tempVect.push_back(tempName);
      }
      columnNames = tempVect;
      getline(textFile, currentLine);
    }
    while (!textFile.eof()) {
      std::stringstream tempStream;
      tempStream << currentLine;
      VecStr tempVect;
      while (!tempStream.eof()) {
        std::string tempName;
        tempStream >> tempName;
        tempVect.push_back(tempName);
      }
      getline(textFile, currentLine);
      tempFileContent.push_back(tempVect);
    }
  } else {
    std::string currentLine;
    std::ifstream textFile(filename.c_str());
    if (!textFile) {
      std::cout << "Error in opening " << filename << std::endl;
      exit(1);
    }
    getline(textFile, currentLine);
    if (header) {
      VecStr temp = tokenizeString(currentLine, inDelim, true);
      columnNames = temp;
      getline(textFile, currentLine);
    }
    while (!textFile.eof()) {
      VecStr temp = tokenizeString(currentLine, inDelim, true);
      tempFileContent.push_back(tokenizeString(currentLine, inDelim));
      getline(textFile, currentLine);
    }
  }
  if (header) {
    *this = table(tempFileContent, columnNames);
  } else {
    *this = table(tempFileContent);
  }
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
        blankPositions.push_back(pos);
      } else {
        nonBlanks.push_back(fIter[colPos]);
      }
      ++pos;
    }
    if (vectorOfNumberStringsDouble(nonBlanks)) {
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
		std::cout << "new column's size doesn't match the size of the table not adding" << std::endl;
		std::cout << "tableSize: " << content_.size() << std::endl;
		std::cout << "addColumnSize: " << columnNew.size() << std::endl;
	}else{
		columnNames_.emplace_back(name);
		if(columnNew.size() == 1){
			for(const auto & rowPos : iter::range(len(content_))){
				content_[rowPos].emplace_back(columnNew[0]);
			}
		}else{
			for(const auto & rowPos : iter::range(len(content_))){
				content_[rowPos].emplace_back(columnNew[rowPos]);
			}
		}
	}
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
table table::getColumns(const VecStr &specificColumnNames) {
  std::vector<uint32_t> positions =
      getPositionsMultipleTargets(columnNames_, specificColumnNames);
  return getColumns(positions);
}

table table::getColumns(const std::vector<uint32_t> &specificColumnPositions) {
  std::vector<VecStr> ans;
  for (const auto &fileIter : content_) {
    VecStr currentRow;
    ans.push_back(getTargetsAtPositions(fileIter, specificColumnPositions));
  }
  table outTab(ans, getTargetsAtPositions(columnNames_, specificColumnPositions));
  outTab.hasHeader_ = hasHeader_;
  return outTab;
}
table table::getColumnsNotAtPositions(
    const std::vector<uint32_t> &specificColumnPositions) {
  std::vector<VecStr> ans;
  for (const auto &fileIter : content_) {
    VecStr currentRow;
    ans.push_back(getTargetsNotAtPositions(fileIter, specificColumnPositions));
  }
  return table(ans,
               getTargetsNotAtPositions(columnNames_, specificColumnPositions));
}

VecStr table::getColumn(const std::string &specifcColumnName) const {
  uint32_t colPos = getFirstPositionOfTarget(columnNames_, specifcColumnName);
  return getColumn(colPos);
}
VecStr table::getColumn(uint32_t colPos) const {
  VecStr ans;
  if (colPos >= columnNames_.size()) {
    std::cout << "positions: " << colPos << " is out of the bounds of "
              << columnNames_.size() << " return nothing" << std::endl;
  } else {
    for (const auto &fIter : content_) {
      ans.push_back(fIter[colPos]);
    }
  }
  return ans;
}

std::vector<std::string *> table::getColumnPointer(
    const std::string &specifcColumnName) {
  size_t colPos = getFirstPositionOfTarget(columnNames_, specifcColumnName);
  std::vector<std::string *> ans;
  for (auto &fIter : content_) {
    ans.push_back(&fIter[colPos]);
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
              << std::endl;
  }
}
void table::deleteRow(size_t rowIndex) {
  if (rowIndex < content_.size()) {
    content_.erase(content_.begin() + rowIndex);
  } else {
    std::cout << "Row " << rowIndex << " is out of bounds, table left alone"
              << std::endl;
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
                               const std::string &startsWtih) {
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions =
      getPositionsOfTargetStartsWith(col, startsWtih);

  return table(getTargetsAtPositions(content_, positions), columnNames_);
}

table table::getRows(const std::vector<uint32_t> &specificRowPositions) {

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
    out << vectorToString(columnNames_, delim) << std::endl;
  }
  outputVectorOfVectors(content_, delim, out);
}

void table::sortTable(const std::string &byThisColumn, bool decending,
                      bool guessFormat) {
  if (!vectorContains(columnNames_, byThisColumn)) {
    std::cout << "Table does not contain " << byThisColumn
              << " not sorting table" << std::endl;
    std::cout << "options are" << std::endl;
    printVector(columnNames_, ",");
    return;
  }
  VecStr col = getColumn(byThisColumn);
  if (guessFormat && vectorOfNumberStringsDouble(col)) {
    if (vectorOfNumberStringsInt(col)) {
      std::multimap<int, VecStr> mapSorter;
      size_t count = 0;
      for (const auto &iter : content_) {
        mapSorter.insert(std::make_pair(std::stoi(col[count]), iter));
        ++count;
      }
      content_.clear();
      for (const auto &iter : mapSorter) {
        content_.push_back(iter.second);
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
        content_.push_back(iter.second);
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
      content_.push_back(iter->second);
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
  std::cout << factorNames << std::endl;
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
      content_.push_back(empty);
    } else {
      addOtherVec(content_[count], tabIter);
    }
    ++count;
  }
  addOtherVec(columnNames_, otherTable.columnNames_);
}

table table::getColumnsLoose(const std::string &subStr) {
  std::vector<uint32_t> positions =
      getPositionsOfSubStrTarget(columnNames_, subStr);
  if (positions.size() == 0) {
    std::cout << "Couldn't find any columns that contained: " << subStr
              << " returning the whole table" << std::endl;
    return *this;
  } else {
    return getColumns(positions);
  }
}

table table::getColumnsStartWith(const std::string &startsWith) {
  std::vector<uint32_t> positions =
      getPositionsOfTargetStartsWith(columnNames_, startsWith);
  if (positions.size() == 0) {
    std::cout << "Couldn't find any columns that started with: " << startsWith
              << " returning the whole table" << std::endl;
    return *this;
  } else {
    return getColumns(positions);
  }
}

table table::getUniqueRows() {
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
    out << iter.first << std::endl;
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
  outTable.columnNames_.push_back("stat");
  for (auto sIter : {"sum", "mean", "median", "min", "max", "std"}) {
    outTable.content_.push_back(VecStr{sIter});
  }
  for (auto &colNameIter : columnNames_) {
    VecStr currentColumn = getColumn(colNameIter);
    if (vectorOfNumberStringsDouble(currentColumn)) {
      outTable.columnNames_.push_back(colNameIter);
      std::vector<double> currentNumbers;
      for (const auto &numStr : currentColumn) {
        currentNumbers.push_back(std::stod(numStr));
      }
      auto currentStats = getStatsOnVec(currentNumbers);
      for (auto i : iter::range(0, (int)outTable.content_.size())) {
        if (i == 0) {
          outTable.content_[i].push_back(std::to_string(currentStats["sum"]));
        } else if (i == 1) {
          outTable.content_[i].push_back(std::to_string(currentStats["mean"]));
        } else if (i == 2) {
          outTable.content_[i]
              .push_back(std::to_string(currentStats["median"]));
        } else if (i == 3) {
          outTable.content_[i].push_back(std::to_string(currentStats["min"]));
        } else if (i == 4) {
          outTable.content_[i].push_back(std::to_string(currentStats["max"]));
        } else if (i == 5) {
          outTable.content_[i].push_back(std::to_string(currentStats["std"]));
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
    outTable.columnNames_.push_back(columnName);
    VecStr otherColNames = stats.begin()->second.columnNames_;
    std::string remove = "stat";
    removeElement(otherColNames, remove);
    addOtherVec(outTable.columnNames_, otherColNames);
    for (const auto &statsIter : stats) {
      VecStr currentLine;
      currentLine.push_back(statsIter.first);
      table info = statsIter.second.getRows("stat", function);
      VecStr currentInfo = info.content_.front();
      removeElement(currentInfo, function);
      addOtherVec(currentLine, currentInfo);
      outTable.content_.push_back(currentLine);
    }
    return outTable;
  } else {
    std::cout << "unrecognized function for aggregate," << function
              << std::endl;
    std::cout << "Available optins are sum, mean, median, max, min, std"
              << std::endl;
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
    if (vectorOfNumberStringsDouble(currentColumn)) {
      ans.push_back(pos);
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
    collapsedColumns.push_back({vectorToString(row, sep)});
  }

  table tempTable(collapsedColumns, {"collapsed"});
  tempTable.cbind(numericColumns);
  table aggregated = tempTable.aggregateSimple("collapsed", function);
  VecStr collapsedColumn = aggregated.getColumn("collapsed");
  aggregated.deleteColumn("collapsed");
  std::vector<VecStr> expanded;
  for (const auto &element : collapsedColumn) {
    expanded.push_back(tokenizeString(element, sep));
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
  ans.columnNames_.push_back(columnForceMatch);
  for (const auto &name : names) {
    ans.content_.push_back({name});
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
    std::cout << "No column named, "
    		<< bib::bashCT::bold << columnName
    		<< bib::bashCT::reset << std::endl;
    std::cout << "Options are : "; printVector(columnNames_, ", ");
    exit(1);
  }
  return countColumn(colPos.front());
}
std::map<std::string, uint32_t> table::countColumn(uint32_t colPos) {
  std::map<std::string, uint32_t> ans;
  if (colPos > len(columnNames_)) {
    std::cout << "colPos: " << colPos << ", out of bounds of "
              << len(columnNames_);
    exit(1);
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
	for(const auto & rowPos : iter::range(len(content_))){
		if(content_.empty()){
			emptyPositions.emplace_back(rowPos);
		}else {
			bool empty =true;
			for(const auto & colPos : iter::range(len(content_[rowPos]))){
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
		sort(emptyPositions);
		for(const auto & pos : iter::reverse(emptyPositions)){
			content_.erase(content_.begin() + pos);
		}
	}

}

table table::extractNumColGreater(const std::string & colName, double cutOff){
  auto comp = [&](const std::string & str){
  	double numValue = std::stod(str);
  	return numValue > cutOff;
  };
  return extractByComp(colName, comp);
}
table table::extractNumColGreater(uint32_t colPos, double cutOff){
  auto comp = [&](const std::string & str){
  	double numValue = std::stod(str);
  	return numValue > cutOff;
  };
  return extractByComp(colPos, comp);
}

}  // namespace bib
