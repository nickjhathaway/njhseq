#include "table.hpp"
#include "njhseq/IO/fileUtils.hpp"
#include "njhseq/IO/InputStream.hpp"
#include "njhseq/IO/OutputStream.hpp"
#include <njhcpp/bashUtils.h>

#include <utility>
#include "njhseq/objects/Meta/MetaDataInName.hpp"


//
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



//adding rows

void table::addRows(const std::vector<VecStr> & rows){
	for(const auto & row : rows){
		addRow(row);
	}
}

void table::addRowFill(VecStr row, const std::string & fill){
  if(len(row) > nCol()){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << "row length: " << len(row) << " is greater than number of columns: " << nCol() << "\n";
    throw std::runtime_error { ss.str() };
  }
  if(len(row) < nCol()){
    addOtherVec(row, std::vector<std::string>(nCol() - len(row), fill));
  }
  content_.emplace_back(row);
}


void table::addRow(const VecStr & row){
	if(len(row) != nCol()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error size of adding row doesn't match column number" << std::endl;
		ss << "Row size: " << len(row) << std::endl;
		ss << "Number of columns " << nCol() << std::endl;
		throw std::runtime_error{ss.str()};
	}
	content_.emplace_back(row);
}

//get lengths
uint32_t table::nCol()const{
	return columnNames_.size();
}

uint32_t table::nRow()const{
	return content_.size();
}


bool table::empty()const{
	return content_.empty();
}


table::table(const TableIOOpts & opts) :
		table(opts.in_.inFilename_.string(), opts.inDelim_, opts.hasHeader_) {

}


void table::addSingleValueColumns(const VecStr & columnValues,
		const VecStr & columnNames) {
	if (columnValues.size() != columnNames.size()) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< ", columnValues size doesn't equal columnNames size" << std::endl;
		ss << "columnValues size: " << columnValues.size() << ", columnNames size: "
				<< columnNames.size() << std::endl;
		throw std::runtime_error { ss.str() };
	} else {
		for (auto pos : iter::range(columnNames.size())) {
			addColumn(VecStr { columnValues[pos] }, columnNames[pos]);
		}
	}
}

bool table::containsColumn(const std::string & colName)const{
	return njh::in(colName, columnNames_);
}



bool table::containsColumn(const std::string & colName,
		std::function<bool(const std::string&, const std::string &)> comp) const {
	return njh::has(columnNames_, colName, comp);
}

bool table::containsColumns(const VecStr & colNames)const{
	for(const auto & col : colNames){
		if(!containsColumn(col)){
			return false;
		}
	}
	return true;
}

bool table::containsColumn(const uint32_t & colPos)const{
	return colPos <= columnNames_.size();
}

void table::setRowSize(uint32_t rowSize){
	for (auto i : iter::range(rowSize)) {
		columnNames_.emplace_back("col." + leftPadNumStr(i, rowSize));
	}
	setColNamePositions();
}

void table::setColNamePositions(){
	colNameToPos_.clear();
	for(const auto pos : iter::range(columnNames_.size())){
		colNameToPos_[columnNames_[pos]] = pos;
	}
}
uint32_t table::getColPos(const std::string & colName) const {
	auto search = colNameToPos_.find(colName);
	if (search == colNameToPos_.end()) {
		throw std::runtime_error { njh::bashCT::boldRed(
				"No column " + colName + " in header," + "available colnames are: "
						+ njh::bashCT::blue + njh::conToStr(columnNames_, ",")) };
	}
	return search->second;
}

VecStr table::getColumnLevels(uint32_t colPos) const {
	std::set<std::string> retSet;
	for (const auto & row : content_) {
		retSet.emplace(row[colPos]);
	}
	return VecStr(retSet.begin(), retSet.end());
}

VecStr table::getColumnLevels(const std::string & colName) const {
	return getColumnLevels(getColPos(colName));
}

void table::changeHeaderToLowerCase(){
	njh::for_each(columnNames_,[](std::string & col){
		njh::strToLower(col);
	});
	setColNamePositions();
}

void table::populateTable(std::istream & in, const std::string &inDelim, bool header){
	inDelim_ = inDelim;
	if(inDelim == "tab"){
		inDelim_ = "\t";
	}
	hasHeader_ = header;
	content_.clear();
	columnNames_.clear();

	std::string currentLine = "";
	uint32_t lineCount = 0;
	while(njh::files::crossPlatGetline(in, currentLine)){
		if(lineCount == 0 && header){
			columnNames_ = tokenizeString(currentLine, inDelim_, true);
		}else{
			content_.emplace_back(tokenizeString(currentLine, inDelim_, true));
		}
		++lineCount;
	}
	if(!hasHeader_){
		for (auto i : iter::range(content_[0].size())) {
			columnNames_.emplace_back("col." + leftPadNumStr(i, content_[0].size()));
		}
	}

	addPaddingToEndOfRows();
	setColNamePositions();

}

table::table(std::istream & in, const std::string &inDelim,
        bool header){
	populateTable(in, inDelim, header);
}

table::table(const bfs::path &filename, const std::string &inDelim,bool header) {

	if(filename != "STDIN"){
		if("" == filename){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "input file is blank" << "\n";
			throw std::runtime_error{ss.str()};
		}
		njh::files::checkExistenceThrow(filename, __PRETTY_FUNCTION__);
	}
	InputStream in(filename);
	populateTable(in, inDelim, header);
}

void table::addPaddingToEndOfRows(const std::string & padding) {
  size_t maxRowLength = columnNames_.size();
  for (auto &iter : content_) {
    if (iter.size() > maxRowLength) {
      maxRowLength = iter.size();
    }
  }
  for (auto &iter : content_) {
    if (iter.size() < maxRowLength) {
      addOtherVec(iter, VecStr(maxRowLength - iter.size(), padding));
    }
  }
}

void table::fillWithZeros() {
  for (const auto colPos : iter::range(columnNames_.size())) {
    std::vector<size_t> blankPositions;
    VecStr nonBlanks;
    size_t pos = 0;
    for (const auto &row : content_) {
      if (row[colPos] == "") {
        blankPositions.emplace_back(pos);
      } else {
        nonBlanks.emplace_back(row[colPos]);
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
		std::stringstream ss;
		ss << "new column's size doesn't match the size of the table" << "\n";
		ss << "tableSize: " << content_.size() << "\n";
		ss << "addColumnSize: " << columnNew.size() << "\n";
		throw std::runtime_error{ss.str()};
	}else{
		columnNames_.emplace_back(name);
		if(columnNew.size() == 1){
			for(const auto rowPos : iter::range(content_.size())){
				content_[rowPos].emplace_back(columnNew[0]);
			}
		}else{
			for(const auto rowPos : iter::range(content_.size())){
				content_[rowPos].emplace_back(columnNew[rowPos]);
			}
		}
	}
	setColNamePositions();
	return;
}

void table::padWithZeros() {
  addPaddingToEndOfRows("0");
}

table table::getColumns(const VecStr &specificColumnNames) const{
	VecStr missing;
	for(const auto & col : specificColumnNames){
		if(!njh::in(col, columnNames_)){
			missing.emplace_back(col);
		}
	}
	if(!missing.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, table doesn't contain the following columns: " << vectorToString(missing, ",") << std::endl;
		ss << "Options are: " << vectorToString(columnNames_, "," ) << std::endl;
		throw std::runtime_error{ss.str()};
	}
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
  	ss << njh::bashCT::bold
				<< "Can't find col: " << njh::bashCT::red << specifcColumnName
				<< njh::bashCT::reset << "\n";
  	ss << "options are : " << njh::bashCT::bold << vectorToString(columnNames_, ",") << njh::bashCT::reset << "\n";
  	throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
  }
  return getColumn(colPos);
}
VecStr table::getColumn(uint32_t colPos) const {
  VecStr ans;
  if (colPos >= columnNames_.size()) {
  	std::stringstream ss;
  	ss << "positions: " << colPos << " is out of the bounds of "
              << columnNames_.size() << " return nothing" << "\n";
  	throw std::runtime_error{njh::bashCT::boldRed(ss.str())};
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
    setColNamePositions();
  } else {
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "column " << columnIndex << " isn't in table, " << "cols range from 0 to " << columnNames_.size()<< "\n";
		throw std::runtime_error{ss.str()};
  }
}
void table::deleteRow(size_t rowIndex) {
  if (rowIndex < content_.size()) {
    content_.erase(content_.begin() + rowIndex);
  } else {
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "row inex " << rowIndex << " out of range, bottom row: " << content_.size()<< "\n";
		throw std::runtime_error{ss.str()};
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
void table::outPutContents(TableIOOpts options) const {
  if ("tab" ==  options.outDelim_ ) {
    options.outDelim_ = "\t";
  } else if ("whitespace" == options.outDelim_) {
    options.outDelim_ = " ";
  }
  if ("" == options.out_.outFilename_) {
    if (options.outOrganized_) {
      outPutContentOrganized(std::cout);
    } else {
      outPutContents(std::cout, options.outDelim_);
    }
  } else {
  	bool writeHeader = hasHeader_ && !( bfs::exists(options.out_.outName()) && options.out_.append_);
    OutputStream outFile(options.out_);
    if (writeHeader) {
    	outFile << vectorToString(columnNames_, options.outDelim_) << "\n";
    }
    outputVectorOfVectors(content_, options.outDelim_, outFile);
    outFile.flush();
  }
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

void table::sortTable(const std::string &firstColumn,
		const std::string &secondColumn, bool decending) {
	uint32_t colPos1 = getColPos(firstColumn);
	bool col1Numeric = isVecOfDoubleStr(getColumn(firstColumn));
	uint32_t colPos2 = getColPos(secondColumn);
	bool col2Numeric = isVecOfDoubleStr(getColumn(secondColumn));
	if(decending){
		njh::sort(content_,[&colPos1,&col1Numeric,&colPos2,&col2Numeric](const VecStr & vec1, const VecStr & vec2){
				if(vec1[colPos1] == vec2[colPos1]){
					if(col2Numeric){
						return std::stod(vec1[colPos2]) > std::stod(vec2[colPos2]);
					}else{
						return vec1[colPos2] > vec2[colPos2];
					}
				}else{
					if(col1Numeric){
						return std::stod(vec1[colPos1]) > std::stod(vec2[colPos1]);
					}else{
						return vec1[colPos1] > vec2[colPos1];
					}
				}
			});
	}else{
		njh::sort(content_,[&colPos1,&col1Numeric,&colPos2,&col2Numeric](const VecStr & vec1, const VecStr & vec2){
				if(vec1[colPos1] == vec2[colPos1]){
					if(col2Numeric){
						return std::stod(vec1[colPos2]) < std::stod(vec2[colPos2]);
					}else{
						return vec1[colPos2] < vec2[colPos2];
					}
				}else{
					if(col1Numeric){
						return std::stod(vec1[colPos1]) < std::stod(vec2[colPos1]);
					}else{
						return vec1[colPos1] < vec2[colPos1];
					}
				}
			});
	}
}

void table::sortTable(const std::string &firstColumn,
		const std::string &secondColumn, const std::string &thirdColumn,
		bool decending) {

	uint32_t colPos1 = getColPos(firstColumn);
	bool col1Numeric = isVecOfDoubleStr(getColumn(firstColumn));
	uint32_t colPos2 = getColPos(secondColumn);
	bool col2Numeric = isVecOfDoubleStr(getColumn(secondColumn));
	uint32_t colPos3 = getColPos(thirdColumn);
	bool col3Numeric = isVecOfDoubleStr(getColumn(thirdColumn));
	if(decending){
		njh::sort(content_,[&colPos1,&col1Numeric,&colPos2,&col2Numeric,&colPos3,&col3Numeric](const VecStr & vec1, const VecStr & vec2){
				if(vec1[colPos1] == vec2[colPos1]){
					if(vec1[colPos2] == vec2[colPos2]){
						if(col3Numeric){
							return std::stod(vec1[colPos3]) > std::stod(vec2[colPos3]);
						}else{
							return vec1[colPos3] > vec2[colPos3];
						}
					}else{
						if(col2Numeric){
							return std::stod(vec1[colPos2]) > std::stod(vec2[colPos2]);
						}else{
							return vec1[colPos2] > vec2[colPos2];
						}
					}
				}else{
					if(col1Numeric){
						return std::stod(vec1[colPos1]) > std::stod(vec2[colPos1]);
					}else{
						return vec1[colPos1] > vec2[colPos1];
					}
				}
			});
	}else{
		njh::sort(content_,[&colPos1,&col1Numeric,&colPos2,&col2Numeric,&colPos3,&col3Numeric](const VecStr & vec1, const VecStr & vec2){
				if(vec1[colPos1] == vec2[colPos1]){
					if(vec1[colPos2] == vec2[colPos2]){
						if(col3Numeric){
							return std::stod(vec1[colPos3]) < std::stod(vec2[colPos3]);
						}else{
							return vec1[colPos3] < vec2[colPos3];
						}
					}else{
						if(col2Numeric){
							return std::stod(vec1[colPos2]) < std::stod(vec2[colPos2]);
						}else{
							return vec1[colPos2] < vec2[colPos2];
						}
					}
				}else{
					if(col1Numeric){
						return std::stod(vec1[colPos1]) < std::stod(vec2[colPos1]);
					}else{
						return vec1[colPos1] < vec2[colPos1];
					}
				}
			});
	}
}
void table::sortTable(const std::string &firstColumn,
		const std::string & secondColumn, const std::string & thirdColumn,
		const std::string & fourthColumn,
		bool decending){
	uint32_t colPos1 = getColPos(firstColumn);
	bool col1Numeric = isVecOfDoubleStr(getColumn(firstColumn));
	uint32_t colPos2 = getColPos(secondColumn);
	bool col2Numeric = isVecOfDoubleStr(getColumn(secondColumn));
	uint32_t colPos3 = getColPos(thirdColumn);
	bool col3Numeric = isVecOfDoubleStr(getColumn(thirdColumn));
	uint32_t colPos4 = getColPos(fourthColumn);
	bool col4Numeric = isVecOfDoubleStr(getColumn(fourthColumn));
	if(decending){
		njh::sort(content_,[&colPos1,&col1Numeric,&colPos2,&col2Numeric,&colPos3,&col3Numeric,&colPos4,&col4Numeric](const VecStr & vec1, const VecStr & vec2){
				if(vec1[colPos1] == vec2[colPos1]){
					if(vec1[colPos2] == vec2[colPos2]){
						if(vec1[colPos3] == vec2[colPos3]){
							if(col4Numeric){
								return std::stod(vec1[colPos4]) > std::stod(vec2[colPos4]);
							}else{
								return vec1[colPos4] > vec2[colPos4];
							}
						}else{
							if(col3Numeric){
								return std::stod(vec1[colPos3]) > std::stod(vec2[colPos3]);
							}else{
								return vec1[colPos3] > vec2[colPos3];
							}
						}
					}else{
						if(col2Numeric){
							return std::stod(vec1[colPos2]) > std::stod(vec2[colPos2]);
						}else{
							return vec1[colPos2] > vec2[colPos2];
						}
					}
				}else{
					if(col1Numeric){
						return std::stod(vec1[colPos1]) > std::stod(vec2[colPos1]);
					}else{
						return vec1[colPos1] > vec2[colPos1];
					}
				}
			});
	}else{
		njh::sort(content_,[&colPos1,&col1Numeric,&colPos2,&col2Numeric,&colPos3,&col3Numeric,&colPos4,&col4Numeric](const VecStr & vec1, const VecStr & vec2){
				if(vec1[colPos1] == vec2[colPos1]){
					if(vec1[colPos2] == vec2[colPos2]){
						if(vec1[colPos3] == vec2[colPos3]){
							if(col4Numeric){
								return std::stod(vec1[colPos4]) < std::stod(vec2[colPos4]);
							}else{
								return vec1[colPos4] < vec2[colPos4];
							}
						}else{
							if(col3Numeric){
								return std::stod(vec1[colPos3]) < std::stod(vec2[colPos3]);
							}else{
								return vec1[colPos3] < vec2[colPos3];
							}
						}
					}else{
						if(col2Numeric){
							return std::stod(vec1[colPos2]) < std::stod(vec2[colPos2]);
						}else{
							return vec1[colPos2] < vec2[colPos2];
						}
					}
				}else{
					if(col1Numeric){
						return std::stod(vec1[colPos1]) < std::stod(vec2[colPos1]);
					}else{
						return vec1[colPos1] < vec2[colPos1];
					}
				}
			});
	}
}
void table::naturlSortTable(const std::string &byThisColumn, bool decending){
	if (!vectorContains(columnNames_, byThisColumn)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "Table does not contain " << byThisColumn
			 << " not sorting table" << "\n";
		ss << "options are: " << vectorToString(columnNames_, ",") << "\n";
		printVector(columnNames_, ",", ss);
		ss << "\n";
		throw std::runtime_error{ss.str()};
	}

	struct NameWithNameSplit {
		explicit NameWithNameSplit(std::string  name):
						name_(std::move(name)){
			const std::regex regPat_{"([A-Za-z0-9\\.]+)" };
			const std::regex subPat_ {"([A-Za-z]*)([0-9\\.]*)"};

			nameToks_ = njh::tokStrOnMatchRegex(name_, regPat_);
			for(const auto & nameTok : nameToks_) {
				std::smatch nameMatch;
				if(!std::regex_match(nameTok.begin(), nameTok.end(), nameMatch, subPat_)) {
					//					std::stringstream ss;
					//					ss << __PRETTY_FUNCTION__ << ", error " << nameTok << "nameTok didn't match pattern"<< "\n";
					//					throw std::runtime_error{ss.str()};
					subNameToks_.emplace_back(nameTok, std::numeric_limits<double>::min() );
				} else {
					subNameToks_.emplace_back(nameMatch[1], ("" == nameMatch[2] ? std::numeric_limits<double>::min() :std::stod(nameMatch[2]) ) );
				}
			}
		}
		std::string name_;

		std::vector<std::string> nameToks_;
		std::vector<std::pair<std::string, double>> subNameToks_;
	};

	auto colPos = getColPos(byThisColumn);

	njh::sort(content_, [&colPos](const VecStr & row1, const VecStr & row2) {
		const NameWithNameSplit seq1(row1[colPos]);
		const NameWithNameSplit seq2(row2[colPos]);
		auto smallest = std::min(seq1.nameToks_.size(), seq2.nameToks_.size());
		for(uint32_t pos = 0; pos < smallest; ++pos) {
			if(seq1.subNameToks_[pos].first == seq2.subNameToks_[pos].first) {
				if(seq1.subNameToks_[pos].second != seq2.subNameToks_[pos].second) {
					return seq1.subNameToks_[pos].second < seq2.subNameToks_[pos].second;
				}
			} else {
				return seq1.subNameToks_[pos].first < seq2.subNameToks_[pos].first;
			}
		}
		return seq1.subNameToks_.size() < seq2.subNameToks_.size();
	});
	if(decending){
		njh::reverse(content_);
	}
}

void table::sortTable(const std::string &byThisColumn, bool decending) {
	if (!vectorContains(columnNames_, byThisColumn)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "Table does not contain " << byThisColumn
				<< " not sorting table" << "\n";
		ss << "options are: " << vectorToString(columnNames_, ",") << "\n";
		printVector(columnNames_, ",", ss);
		ss << "\n";
		throw std::runtime_error{ss.str()};
	}
	VecStr col = getColumn(byThisColumn);
	uint32_t colPos = getColPos(byThisColumn);
	if (isVecOfDoubleStr(col)) {
		if (isVecOfIntStr(col)) {
			std::sort(content_.begin(), content_.end(),
					[&colPos](const VecStr & vec1, const VecStr & vec2) {
				return std::stoll(vec1[colPos]) < std::stoll(vec2[colPos]);});
		} else {
			std::sort(content_.begin(), content_.end(),
					[&colPos](const VecStr & vec1, const VecStr & vec2) {
				return std::stod(vec1[colPos]) < std::stod(vec2[colPos]);});
		}
	} else {
		std::sort(content_.begin(), content_.end(),
				[&colPos](const VecStr & vec1, const VecStr & vec2) {
			return vec1[colPos] < vec2[colPos];});
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

void table::rbind(const table &otherTable, bool fill) {
	VecStr missingColsFromThis;
	VecStr missingColsFromOther;
	for(const auto & col : otherTable.columnNames_){
		if(!njh::in(col, columnNames_)){
			missingColsFromThis.emplace_back(col);
		}
	}

	for(const auto & col : columnNames_){
		if(!njh::in(col, otherTable.columnNames_)){
			missingColsFromOther.emplace_back(col);
		}
	}
	auto otherTableCopy = otherTable;

	if(!missingColsFromOther.empty() || !missingColsFromThis.empty()){
		if(fill){
			if(!missingColsFromOther.empty()){
				for(const auto & col : missingColsFromOther){
					otherTableCopy.addColumn({"NA"}, col);
				}
			}
			if(!missingColsFromThis.empty()){
				for(const auto & col : missingColsFromThis){
					addColumn({"NA"}, col);
				}
			}
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error" << "\n";
			if(!missingColsFromOther.empty()){
				ss << "Missing the following columns from adding table: " << njh::conToStrEndSpecial(missingColsFromOther, ", ", " and ") << "\n";
			}
			if(!missingColsFromThis.empty()){
				ss << "Missing the following columns from this table  : " << njh::conToStrEndSpecial(missingColsFromThis, ", ", " and ") << "\n";
			}
			throw std::runtime_error{ss.str()};
		}
	}
	otherTableCopy = otherTableCopy.getColumns(columnNames_);
	for(const auto & row : otherTableCopy){
		addRow(row);
	}
}
void table::cbind(const table &otherTable, bool fill) {
	if (!fill && otherTable.content_.size() != content_.size()) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", adding a table that has a different number of rows" << std::endl;
		ss << "This table row number: " << content_.size()
				<< ", other table row number: " << otherTable.content_.size()
				<< std::endl;
		ss << "To fill in the empty rows use fill = true" << std::endl;
		throw std::runtime_error { ss.str() };
	}
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

table table::getColumnsMatchingPattern(const std::regex & pattern) const {
	std::vector<uint32_t> positions = getPositionsMatchingPattern(columnNames_,
			pattern);
	if (positions.size() == 0) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "Couldn't find any columns that matched pattern"<< "\n";
		throw std::runtime_error{ss.str()};
	}
	return getColumns(positions);
}

table table::getColumnsContainingPattern(const std::regex & pattern) const {
	std::vector<uint32_t> positions = getPositionsContainingPattern(columnNames_,
			pattern);
	if (positions.size() == 0) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " Couldn't find any columns that contain pattern"<< "\n";
		throw std::runtime_error{ss.str()};
	}
	return getColumns(positions);
}

table table::getColumnsStartWith(const std::string &startsWith) const{
  std::vector<uint32_t> positions =
      getPositionsOfTargetStartsWith(columnNames_, startsWith);
  if (positions.size() == 0) {
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " Couldn't find any columns that started with: " << startsWith << "\n";
		throw std::runtime_error{ss.str()};
  }
  return getColumns(positions);
}

table table::getColumnsNotMatchingPattern(const std::regex & pattern) const {
	std::vector<uint32_t> positions = getPositionsMatchingPattern(columnNames_,
			pattern);
	if (positions.size() == columnNames_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "Couldn't find any columns that matched pattern"<< "\n";
		throw std::runtime_error{ss.str()};
	}
	std::vector<uint32_t> outPositions;
	for(const auto pos : iter::range(columnNames_.size())){
		if(!njh::in<uint32_t>(pos, positions)){
			outPositions.emplace_back(pos);
		}
	}
	return getColumns(outPositions);
}

table table::getColumnsNotContainingPattern(const std::regex & pattern) const {
	std::vector<uint32_t> positions = getPositionsContainingPattern(columnNames_,
			pattern);
	if (positions.size() == columnNames_.size()) {
    std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " Couldn't find any columns that contain pattern " << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::vector<uint32_t> outPositions;
	for(const auto pos : iter::range(columnNames_.size())){
		if(!njh::in<uint32_t>(pos, positions)){
			outPositions.emplace_back(pos);
		}
	}
	return getColumns(positions);
}

table table::getColumnsNotStartWith(const std::string &startsWith) const{
  std::vector<uint32_t> positions =
      getPositionsOfTargetStartsWith(columnNames_, startsWith);
  if (positions.size() == columnNames_.size()) {
  	std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " Couldn't find any columns that started with: " << startsWith<< "\n";
		throw std::runtime_error{ss.str()};
  }
	std::vector<uint32_t> outPositions;
	for(const auto pos : iter::range(columnNames_.size())){
		if(!njh::in<uint32_t>(pos, positions)){
			outPositions.emplace_back(pos);
		}
	}
  return getColumns(positions);
}




table table::getRowsNotContainingPattern(const std::string &forColumn,
		const std::regex & pattern) const {
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions = getPositionsContainingPattern(col, pattern);
  return table(getTargetsNotAtPositions(content_, positions), columnNames_);
}

table table::getRowsNotMatchingPattern(const std::string &forColumn,
		const std::regex & pattern) const {
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions = getPositionsMatchingPattern(col, pattern);
  return table(getTargetsNotAtPositions(content_, positions), columnNames_);
}

table table::getRowsContainingPattern(const std::string &forColumn,
		const std::regex & pattern) const{
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions = getPositionsContainingPattern(col, pattern);
  return table(getTargetsAtPositions(content_, positions), columnNames_);
}

table table::getRowsMatchingPattern(const std::string &forColumn,
		const std::regex & pattern) const{
  VecStr col = getColumn(forColumn);
  std::vector<uint32_t> positions = getPositionsMatchingPattern(col, pattern);
  return table(getTargetsAtPositions(content_, positions), columnNames_);
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
  for (auto sIter : VecStr{"sum", "mean", "median", "min", "max", "std"}) {
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
      fillWithZeros();
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



std::vector<uint32_t> table::getNumericColumnPositions(bool addZeros) {
  if (addZeros) {
  	fillWithZeros();
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
                              const std::string &function,
															bool addZeros) {
  std::string sep = "______";
  table numericColumns = getAllNumericColumns();
  table nonNumericColumns = getAllNonNumericColumns();
  std::vector<VecStr> collapsedColumns;
  for (const auto & row : nonNumericColumns.content_) {
    collapsedColumns.emplace_back(VecStr{vectorToString(row, sep)});
  }

  table tempTable(collapsedColumns, {"collapsed"});
  tempTable.cbind(numericColumns, true);
  table aggregated = tempTable.aggregateSimple("collapsed", function);
  VecStr collapsedColumn = aggregated.getColumn("collapsed");
  aggregated.deleteColumn("collapsed");
  std::vector<VecStr> expanded;
  for (const auto &element : collapsedColumn) {
    expanded.emplace_back(tokenizeString(element, sep));
  }
  table expandedTable(expanded, nonNumericColumns.columnNames_);
  expandedTable.cbind(aggregated,true);
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
    ans.fillWithZeros();
  }
  return ans;
}
std::map<std::string, uint32_t> table::countColumn(
    const std::string &columnName) {
  auto colPos = getPositionsOfTarget(columnNames_, columnName);
  if (colPos.empty()) {
  	std::stringstream ss;
    ss << "No column named, "
    		<< njh::bashCT::bold << columnName
    		<< njh::bashCT::reset << "\n";
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
	for(const auto rowPos : iter::range(content_.size())){
		if(content_.empty()){
			emptyPositions.emplace_back(rowPos);
		}else {
			bool empty =true;
			for(const auto colPos : iter::range(content_[rowPos].size())){
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
		njh::sort(emptyPositions);
		for(const auto & pos : iter::reversed(emptyPositions)){
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
		ss << njh::bashCT::red << njh::bashCT::bold
				<< "Error in table::countColumn(const std::vector<uint32_t> & colPositions), out of range error\n"
				<< njh::bashCT::blue << vectorToString(colPositions, ",")
				<< njh::bashCT::red << " out of range of " << columnNames_.size()
				<< njh::bashCT::reset << "\n";
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
  auto comp = [&cutOff](const std::string & str){
  	double numValue = std::stod(str);
  	return numValue > cutOff;
  };
  return extractByComp(colName, comp);
}
table table::extractNumColGreater(uint32_t colPos, double cutOff)const{
  auto comp = [&cutOff](const std::string & str){
  	double numValue = std::stod(str);
  	return numValue > cutOff;
  };
  return extractByComp(colPos, comp);
}


void table::checkForColumnsThrow(const VecStr & requiredColumns, const std::string & funcName) const{
	VecStr columnsNotFound;
	for (const auto & col : requiredColumns) {
		if (!njh::in(col, columnNames_)) {
			columnsNotFound.emplace_back(col);
		}
	}
	if (!columnsNotFound.empty()) {
		std::stringstream ss;
		ss << "Need to have " << njh::conToStrEndSpecial(requiredColumns, ",", " and ") << '\n';
		ss << "Did not find " << njh::conToStrEndSpecial(columnsNotFound, ",", " or ") << '\n';
		ss << "Only have " << '\n';
		uint32_t colCount = 1;
		for(const auto & col : columnNames_){
			ss << "\tcolum " << colCount << " is " << col << '\n';
			++colCount;
		}
		throw std::runtime_error { ss.str() };
	}
}


VecStr table::getMissingHeaders(const VecStr requiredColumns) const{
	VecStr columnsNotFound;
	for (const auto & col : requiredColumns) {
		if (!njh::in(col, columnNames_)) {
			columnsNotFound.emplace_back(col);
		}
	}
	return columnsNotFound;
}

table table::splitColWithMeta(const table & inputTab, const splitColWithMetaPars & pars){
	auto inTab = inputTab;
	inTab.checkForColumnsThrow({pars.column_}, __PRETTY_FUNCTION__);
	std::set<std::string> metaFields;
	std::unordered_map<std::string, VecStr> metaValues;

	bool noneContainMeta = true;
	for (const auto & row : inTab.content_) {
		if(MetaDataInName::nameHasMetaData(row[inTab.getColPos(pars.column_)])){
			MetaDataInName rowMeta(row[inTab.getColPos(pars.column_)]);
			auto metas = getVectorOfMapKeys(rowMeta.meta_);
			std::copy(metas.begin(), metas.end(),
					std::inserter(metaFields, metaFields.end()));
			noneContainMeta = false;
		}
	}
	if(!noneContainMeta){
		bool allEmpty = true;
		for (auto & row : inTab.content_) {
			if(MetaDataInName::nameHasMetaData(row[inTab.getColPos(pars.column_)])){
				MetaDataInName rowMeta(row[inTab.getColPos(pars.column_)]);
				for(const auto & m : rowMeta.meta_){
					metaValues[m.first].emplace_back(m.second);
				}
				for(const auto & metaField : metaFields){
					if(!njh::in(metaField, rowMeta.meta_)){
						metaValues[metaField].emplace_back("NA");
					}
				}
				if(!pars.keepMetaInColumn_){
					MetaDataInName::removeMetaDataInName(row[inTab.getColPos(pars.column_)]);
					if("" != row[inTab.getColPos(pars.column_)]){
						allEmpty = false;
					}
				}
			}else{
				for(const auto & metaField : metaFields){
					metaValues[metaField].emplace_back("NA");
				}
			}
		}
		if (allEmpty && pars.removeEmptyColumn_) {
			inTab.deleteColumn(pars.column_);
		}
		for (const auto & m : metaValues) {
			std::string newColName = m.first;
			if (pars.prefixWithColName_) {
				newColName = pars.column_ + "-" + m.first;
			}
			inTab.addColumn(m.second, newColName);
		}
		if(pars.sorting_){
			inTab.sortTable(pars.sortCol_, pars.descending_);
		}
	}else{
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": No values in " << pars.column_ << " contain meta data, doing nothing" << "\n";
		throw std::runtime_error{ss.str()};
	}
	return inTab;
}

}  // namespace njh
