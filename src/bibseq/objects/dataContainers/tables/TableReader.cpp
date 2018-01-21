/*
 * TableReader.cpp
 *
 *  Created on: Jan 21, 2018
 *      Author: nick
 */


#include "TableReader.hpp"

namespace bibseq {

TableReader::TableReader(const TableIOOpts & tabOpts): tabOpts_(tabOpts){
	//inital header reader
	bib::files::checkExistenceThrow(tabOpts_.in_.inFilename_);
	std::string currentLine = bib::files::getFirstLine(
			tabOpts_.in_.inFilename_);
	auto toks = tokenizeString(currentLine, tabOpts_.inDelim_, true);
	VecStr columnNames;
	if (!tabOpts_.hasHeader_) {
		for (const auto i : iter::range(toks.size())) {
			columnNames.emplace_back("col." + leftPadNumStr(i, toks.size()));
		}
	} else {
		columnNames = toks;
	}
	header_ = table(columnNames);
	in_ = std::make_unique<InputStream>(tabOpts_.in_);
	if(tabOpts_.hasHeader_){
		bib::files::crossPlatGetline(*in_, currentLine);
	}
}


bool TableReader::getNextRow(VecStr & row){
	std::string currentLine = "";
	row.clear();
	if(bib::files::crossPlatGetline(*in_, currentLine)){
		row = tokenizeString(currentLine, tabOpts_.inDelim_, true);
		if(row.size() != header_.nCol()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error the row has a different number of columns than the first line" << "\n";
			ss << "rowSize: " << row.size() << ", firstLineSize: " << header_.nCol() << "\n";
			ss << "row: " << currentLine << "\n";
			throw std::runtime_error{ss.str()};
		}
		return true;
	}
	return false;
}

VecStr TableReader::extractCols(const VecStr & row, const VecStr & cols) const{
	VecStr ret;
	ret.reserve(cols.size());
	for(const auto & col : cols){
		ret.emplace_back(row[header_.getColPos(col)]);
	}
	return ret;
}

}  // namespace bibseq
