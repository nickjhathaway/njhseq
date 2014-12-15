#pragma once
//
//  table.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/7/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include <cppitertools/range.hpp>
namespace bibseq {

struct outOptions {
  outOptions()
      : outFilename_(""),
        extention_(".txt"),
        outDelim_("\t"),
        outOrganized_(false),
        overWriteFile_(false),
        exitOnFailureToWrite_(false) {}
  outOptions(const std::string & outFilename,
  		const std::string & extention)
        : outFilename_(outFilename),
          extention_(extention),
          outDelim_("\t"),
          outOrganized_(false),
          overWriteFile_(false),
          exitOnFailureToWrite_(false) {}
  outOptions(const std::string & outFilename,
  		const std::string & extention, const std::string & outDelim_)
        : outFilename_(outFilename),
          extention_(extention),
          outDelim_(outDelim_),
          outOrganized_(false),
          overWriteFile_(false),
          exitOnFailureToWrite_(false) {}
  outOptions(const std::string & outFilename,
  		const std::string & extention, const std::string & outDelim_,
  		bool overWriteFile)
        : outFilename_(outFilename),
          extention_(extention),
          outDelim_(outDelim_),
          outOrganized_(false),
          overWriteFile_(overWriteFile),
          exitOnFailureToWrite_(false) {}

  std::string outFilename_;
  std::string extention_;
  std::string outDelim_;
  bool outOrganized_;
  bool overWriteFile_;
  bool exitOnFailureToWrite_;
};

class table {

 public:
  // constructors
  table() : hasHeader_(false) {}
  table(const std::vector<VecStr> &inContent, const VecStr &inColumnNames)
      : content_(inContent) {
    size_t pos = 0;
    for (const auto &cIter : inColumnNames) {
      colNames_.insert(std::make_pair(cIter, pos));
      ++pos;
    }
    columnNames_ = inColumnNames;
    hasHeader_ = true;
    addPaddingToEndOfRows();
  }
  table(const std::vector<VecStr> &inContent) : content_(inContent) {
    hasHeader_ = false;
    addPaddingToEndOfRows();
    for (auto i : iter::range(0, (int)content_[0].size())) {
      columnNames_.push_back(
          combineStrings({"col.", leftPadNumStr(i, (int)content_[0].size())}));
    }
  }
  table(std::istream & in, const std::string &inDelim = " ",
          bool header = false);

  table(const std::string &filename, const std::string &inDelim = " ",
        bool header = false);
  /*
   * Construct with simple unordered_map
   */
  template <typename KEY, typename VALUE>
  table(const std::unordered_map<KEY, VALUE> &inMap, const VecStr &columnNames)
      : hasHeader_(true), columnNames_(columnNames) {
    for (const auto &row : inMap) {
      content_.emplace_back(pairToVecStr(row));
    }
  }
  template <typename MAP1KEY, typename MAP2KEY, typename MAP2VALUE>
  table(const std::map<MAP1KEY, std::map<MAP2KEY, MAP2VALUE>> &input,
        const VecStr &colNames)
      : hasHeader_(true), columnNames_(colNames) {
    // std::cout <<"at least started" << std::endl;
    for (const auto &m1 : input) {
      // std::cout << "m1.first" << m1.first << std::endl;
      for (const auto &m2 : m1.second) {
        content_.emplace_back(VecStr{to_string(m1.first), to_string(m2.first),
                                     to_string(m2.second)});
      }
    }
  }
  template <typename MAP1KEY, typename MAP2KEY, typename MAP2VALUE>
  table(const std::unordered_map<MAP1KEY,
                                 std::unordered_map<MAP2KEY, MAP2VALUE>> &input,
        const VecStr &colNames)
      : hasHeader_(true), columnNames_(colNames) {
    // std::cout <<"at least started" << std::endl;
    for (const auto &m1 : input) {
      // std::cout << "m1.first" << m1.first << std::endl;
      for (const auto &m2 : m1.second) {
        content_.emplace_back(VecStr{to_string(m1.first), to_string(m2.first),
                                     to_string(m2.second)});
      }
    }
  }
  /*
   * Construct with simple map
   */
  template <typename KEY, typename VALUE>
  table(const std::map<KEY, VALUE> &inMap, const VecStr &columnNames)
      : hasHeader_(true), columnNames_(columnNames) {
    for (const auto &row : inMap) {
      content_.emplace_back(pairToVecStr(row));
    }
  }
  /*
   * Construct with just column names, data to be put in latter
   */
  table(const VecStr &columnNames)
      : hasHeader_(true), columnNames_(columnNames) {
  }

  // memebers
  bool hasHeader_;
  std::vector<VecStr> content_;
  VecStr columnNames_;
  std::map<std::string, size_t> colNames_;
  std::string inDelim_;
  // functions

  // to ensure all the rows have equal length
  void addPaddingToEndOfRows();
  void addPaddingZeros();
  void addZerosToEnds();
  void addColumn(const VecStr & columnNew, const std::string & name);
  template<typename T>
  void addColumn(const std::vector<T> & columnNew, const std::string & name){
  	VecStr add = numVecToVecStr(columnNew);
  	addColumn(add, name);
  }
  // outputs
  void outPutContents(outOptions options) const;
  void outPutContents(std::ostream &out, std::string delim) const;
  void outPutContentOrganized(std::ostream &out) const;
  // extracting columns
  table getColumns(const VecStr &specificColumnNames);
  table getColumns(const std::vector<uint32_t> &specificColumnPositions);
  table getColumnsNotAtPositions(
      const std::vector<uint32_t> &specificColumnPositions);
  table getColumnsLoose(const std::string &subStr);
  table getColumnsStartWith(const std::string &startsWith);
  VecStr getColumn(const std::string &specifcColumnName) const;
  VecStr getColumn(uint32_t pos) const;
  std::vector<std::string *> getColumnPointer(
      const std::string &specificColumnName);
  // deleting columns
  void deleteColumn(const std::string &columnName);
  void deleteColumn(size_t columnIndex);
  // extracting rows
  table getRows(const std::string &forColumn, const std::string &element) const;
  table getRows(const std::vector<uint32_t> &specificRowPositions);
  table getRowsLoose(const std::string &forColumn,
                     const std::string &subString) const;
  table getRowsStartsWith(const std::string &forColumn,
                          const std::string &startsWtih);
  // get unique rows only
  table getUniqueRows();
  // deleting a row
  void deleteRow(size_t rowIndex);
  // sort the table
  void sortTable(const std::string &byThisColumn, bool decending,
                 bool guessFormat = true);
  //
  void trimElementsAtFirstOccurenceOf(const std::string &trimAt);
  // split a table on certain elements in a row
  std::map<std::string, table> splitTableOnColumn(const std::string &colName)
      const;
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
  void rbind(const table &otherTable);
  void cbind(const table &otherTable);
  // reversing the table
  void reverseRows();
  void reverseColumns();
  // print out a map of tables
  static void printOutSplitTable(const std::map<std::string, table> &tabSplit,
                                 std::ostream &out, const std::string &delim,
                                 bool outOrganized);
  std::map<std::string, uint32_t> countColumn(const std::string &columnName);
  std::map<std::string, uint32_t> countColumn(uint32_t colPos);
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

  table extractNumColGreater(uint32_t colPos, double cutOff);
  table extractNumColGreater(const std::string & colName, double cutOff);
  template<typename UnaryPredicate>
  table extractByComp(uint32_t colPos, UnaryPredicate p){
  	table out(columnNames_);
  	if(colPos < len(columnNames_)){
  		for(const auto & rowPos : iter::range(len(content_))){
  			if(p(content_[rowPos][colPos])){
  				out.content_.emplace_back(content_[rowPos]);
  			}
  		}
  	}else{
  		std::cout << "columnPos: " << colPos << "is out of bounds of " << len(columnNames_) << std::endl;
  		exit(1);
  	}
  	return out;
  }
  template<typename UnaryPredicate>
  table extractByComp(const std::string & columnName, UnaryPredicate p){
  	if(in(columnName, columnNames_)){
  		uint32_t pos = getFirstPositionOfTarget(columnNames_,columnName );
  		return extractByComp(pos, p);
  	}else{
  		std::cout << "columnName: " << columnName << "doesn't exist" << std::endl;
  		std::cout << "options are : "; printVector(columnNames_, ", ");
  		exit(1);
  		return table();
  	}
  }
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "table.cpp"
#endif
