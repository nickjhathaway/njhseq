/*
 * TableWriter.cpp
 *
 *  Created on: Nov 6, 2021
 *      Author: nick
 */




#include "TableWriter.hpp"

namespace njhseq {
TableWriter::TableWriter(const TableIOOpts & tabOpts, const VecStr & header): tabOpts_(tabOpts),header_(header){
	out_ = std::make_unique<OutputStream>(tabOpts_.out_);
	if(tabOpts.hasHeader_){
		header_.outPutContents(*out_, tabOpts_.outDelim_);
	}
}


void TableWriter::writeRow(VecStr & row){
	if(!doNotCheckRowSizes  && row.size() != header_.nCol()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "row size, "<< row.size()	 << " doesn't match header size, " << header_.nCol()	<< "\n";
		ss << "header: " << njh::conToStr(header_.columnNames_, "\t") << "\n";
		ss << "row: " << njh::conToStr(row, "\t") << "\n";
		throw std::runtime_error{ss.str()};
	}
	*out_ << njh::conToStr(row, tabOpts_.outDelim_) << "\n";
}


void TableWriter::writeRowLock(VecStr & row){
	std::lock_guard<std::mutex> lock(mut_);
	writeRow(row);
}






}  // namespace njhseq
