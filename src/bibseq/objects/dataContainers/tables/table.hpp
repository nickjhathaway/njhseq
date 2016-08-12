#pragma once
//
//  table.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/7/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
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
#include "bibseq/utils.h"
#include "bibseq/IO/IOUtils.hpp"
#include "bibseq/objects/dataContainers/tables/TableIOOpts.hpp"


namespace bibseq {



class table {
public:
	table() :
			hasHeader_(false) {
	}
	/**@b Construct with a vec of vec of strings with column name
	 *
	 * @param content The content to be put into the table
	 * @param columnNames The column names for the table
	 */
	table(const std::vector<VecStr> &content, const VecStr &columnNames) :
			content_(content), columnNames_(columnNames), hasHeader_(true) {
		addPaddingToEndOfRows();
		setColNamePositions();
	}
	/**@b Construct with a vector of vector of strings for content, no headers
	 *
	 * @param inContent
	 */
	table(const std::vector<VecStr> &content) :
			content_(content), hasHeader_(false) {
		addPaddingToEndOfRows();
		for (auto i : iter::range(content_[0].size())) {
			columnNames_.emplace_back("col." + leftPadNumStr(i, content_[0].size()));
		}
		setColNamePositions();
	}
	/**@b Construct with just column names, data to be put in latter
	 *
	 * @param columnNames A vector of column names
	 */
	table(const VecStr &columnNames) :
			columnNames_(columnNames), hasHeader_(true) {
		setColNamePositions();
	}
	/**@b Construct with a stream with lines separated by new line characters and each line is delimited
	 *
	 * @param in In stream
	 * @param inDelim The delimiter per line
	 * @param header Whether the first line is a header
	 */
	table(std::istream & in, const std::string &inDelim = "whitespace",
			bool header = false);
	/**@b Construct with a file with lines separated by new line characters and each line is delimited
	 *
	 * @param in In file name
	 * @param inDelim The delimiter per line
	 * @param header Whether the first line is a header
	 */
	table(const std::string &filename, const std::string &inDelim = "whitespace",
			bool header = false);

	/**@b Construct with a file with lines separated by new line characters and each line is delimited
	 *
	 * @param in In file name
	 * @param inDelim The delimiter per line
	 * @param header Whether the first line is a header
	 */
	table(const TableIOOpts & opts);

	/**@b Construct with simple unordered_map
	 *
	 */
	template<typename KEY, typename VALUE>
	table(const std::unordered_map<KEY, VALUE> &inMap, const VecStr &columnNames) :
			columnNames_(columnNames), hasHeader_(true) {
		for (const auto &row : inMap) {
			content_.emplace_back(pairToVecStr(row));
		}
		setColNamePositions();
	}
	template<typename MAP1KEY, typename MAP2KEY, typename MAP2VALUE>
	table(const std::map<MAP1KEY, std::map<MAP2KEY, MAP2VALUE>> &input,
			const VecStr &colNames) :
			columnNames_(colNames), hasHeader_(true) {
		// std::cout <<"at least started" << std::endl;
		for (const auto &m1 : input) {
			// std::cout << "m1.first" << m1.first << std::endl;
			for (const auto &m2 : m1.second) {
				content_.emplace_back(toVecStr(m1.first, m2.first, m2.second));
			}
		}
		setColNamePositions();
	}
	template<typename MAP1KEY, typename MAP2KEY, typename MAP2VALUE>
	table(
			const std::unordered_map<MAP1KEY, std::unordered_map<MAP2KEY, MAP2VALUE>> &input,
			const VecStr &colNames) :
			columnNames_(colNames), hasHeader_(true) {
		// std::cout <<"at least started" << std::endl;
		for (const auto &m1 : input) {
			// std::cout << "m1.first" << m1.first << std::endl;
			for (const auto &m2 : m1.second) {
				content_.emplace_back(toVecStr(m1.first, m2.first, m2.second));
			}
		}
		setColNamePositions();
	}

	template<typename KEY, typename VALUE>
	table(const std::map<KEY, VALUE> &inMap, const VecStr &columnNames) :
			columnNames_(columnNames), hasHeader_(true) {
		for (const auto &row : inMap) {
			content_.emplace_back(pairToVecStr(row));
		}
		setColNamePositions();
	}

	// members
	std::vector<VecStr> content_;
	VecStr columnNames_;
	std::unordered_map<std::string, uint32_t> colNameToPos_;
	bool hasHeader_;
	std::string inDelim_;

	void setRowSize(uint32_t rowSize);
	void setColNamePositions();
	uint32_t getColPos(const std::string & colName) const;
	bool containsColumn(const std::string & colName) const;
	bool containsColumns(const VecStr & colNames) const;
	bool containsColumn(const uint32_t & colPos) const;
	// to ensure all the rows have equal length
	void addPaddingToEndOfRows(const std::string & padding = "");
	void padWithZeros();
	void fillWithZeros();
	//adding coluns
	void addColumn(const VecStr & columnNew, const std::string & name);
	void addSingleValueColumns(const VecStr & columnValues, const VecStr & columnNames);
	template<typename T>
	void addColumn(const std::vector<T> & columnNew, const std::string & name) {
		VecStr add = numVecToVecStr(columnNew);
		addColumn(add, name);
	}
	//adding rows
	void addRows(const std::vector<VecStr> & rows);
	void addRow(const VecStr & row);
	template<typename... T>
	void addRow(const T&... items){
		if(sizeof...(items) != nCol()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error size of adding row doesn't match column number" << std::endl;
			ss << "Row size: " << sizeof...(items) << std::endl;
			ss << "Number of columns " << nCol() << std::endl;
			throw std::runtime_error{ss.str()};
		}
		addRow(toVecStr(items...));
	}
	//get lengths
	uint32_t nCol() const;
	uint32_t nRow() const;
	bool empty() const;

	// outputs
	void outPutContents(TableIOOpts options) const;
	void outPutContents(std::ostream &out, std::string delim) const;
	void outPutContentOrganized(std::ostream &out) const;
	// extracting columns
	table getColumns(const VecStr &specificColumnNames) const;
	table getColumns(const std::vector<uint32_t> &specificColumnPositions) const;
	table getColumnsNotAtPositions(
			const std::vector<uint32_t> &specificColumnPositions) const;
	table getColumnsLoose(const std::string &subStr) const;
	table getColumnsStartWith(const std::string &startsWith) const;
	VecStr getColumn(const std::string &specifcColumnName) const;
	VecStr getColumn(uint32_t pos) const;
	std::vector<std::string *> getColumnPointer(
			const std::string &specificColumnName);
	// deleting columns
	void deleteColumn(const std::string &columnName);
	void deleteColumn(size_t columnIndex);
	// extracting rows
	table getRows(const std::string &forColumn, const std::string &element) const;
	table getRows(const std::vector<uint32_t> &specificRowPositions) const;
	table getRowsLoose(const std::string &forColumn,
			const std::string &subString) const;
	table getRowsStartsWith(const std::string &forColumn,
			const std::string &startsWtih) const;
	// get unique rows only
	table getUniqueRows() const;
	// deleting a row
	void deleteRow(size_t rowIndex);
	// sort the table
	void sortTable(const std::string &byThisColumn, bool decending);
	void sortTable(const std::string &firstColumn,
			const std::string & secondColumn, bool decending);
	void sortTable(const std::string &firstColumn,
			const std::string & secondColumn, const std::string & thirdColumn,
			bool decending);
	void sortTable(const std::string &firstColumn,
			const std::string & secondColumn, const std::string & thirdColumn,
			const std::string & fourthColumn,
			bool decending);

	//
	void trimElementsAtFirstOccurenceOf(const std::string &trimAt);
	// split a table on certain elements in a row
	std::map<std::string, table> splitTableOnColumn(
			const std::string &colName) const;
	std::map<std::string, table> splitTableOnColumnLoose(
			const std::string &colName) const;
	std::map<std::string, table> splitTableOnColumnLoose(
			const std::string &colName, const std::string &firstOccurnceOf) const;
	std::map<std::string, table> splitTableOnColumnLoose(
			const std::string &colName, const VecStr &useTheseNames) const;
	std::vector<uint32_t> getNumericColumnPositions(bool addZeros = true);
	table getAllNumericColumns(bool addZeros = true);
	table getAllNonNumericColumns(bool addZeros = true);
	// adding another table
	void rbind(const table &otherTable, bool fill);
	void cbind(const table &otherTable, bool fill);
	// reversing the table
	void reverseRows();
	void reverseColumns();
	// print out a map of tables
	static void printOutSplitTable(const std::map<std::string, table> &tabSplit,
			std::ostream &out, const std::string &delim, bool outOrganized);
	std::map<std::string, uint32_t> countColumn(const std::string &columnName);
	std::map<std::string, uint32_t> countColumn(uint32_t colPos);

	table countColumn(const VecStr &columnNames);
	table countColumn(const std::vector<uint32_t> & colPositions);

	table getStatsTable() const;
	table aggregateSimple(const std::string &columnName,
			const std::string &function, bool addZeros = true);
	table aggregateAdvance(const std::string &columnName,
			const std::string &function, bool addZeros = true);
	static table cbind(const std::vector<table> &tables,
			const std::string &columnForceMatch, bool addZeros = true);
	static table cbind(const std::map<std::string, table> &tables,
			const std::string &columnForceMatch, bool addZeros = true);
	void removeEmpty(bool addPadding);

	table extractNumColGreater(uint32_t colPos, double cutOff) const;
	table extractNumColGreater(const std::string & colName, double cutOff) const;

	VecStr getColumnLevels(uint32_t colPos)const;
	VecStr getColumnLevels(const std::string & colName)const;

	template<typename UnaryPredicate>
	table extractByComp(uint32_t colPos, UnaryPredicate p) const {
		table out(columnNames_);
		if (colPos < columnNames_.size()) {
			for (const auto & rowPos : iter::range(content_.size())) {
				if (p(content_[rowPos][colPos])) {
					out.content_.emplace_back(content_[rowPos]);
				}
			}
		} else {
			std::stringstream ss;
			ss << "columnPos: " << colPos << "is out of bounds of "
					<< columnNames_.size() << std::endl;
			throw std::out_of_range { ss.str() };
		}
		return out;
	}
	template<typename UnaryPredicate>
	table extractByComp(const std::string & columnName, UnaryPredicate p) const {
		if (bib::in(columnName, columnNames_)) {
			uint32_t pos = getFirstPositionOfTarget(columnNames_, columnName);
			return extractByComp(pos, p);
		} else {
			std::stringstream ss;
			ss << "columnName: " << columnName << "doesn't exist" << std::endl;
			ss << "options are : ";
			printVector(columnNames_, ", ", ss);
			throw std::runtime_error { ss.str() };
			return table();
		}
	}
	auto begin() const {
		return content_.begin();
	}
	auto end() const {
		return content_.end();
	}

	auto begin() {
		return content_.begin();
	}
	auto end() {
		return content_.end();
	}
};
}  // namespace bib


