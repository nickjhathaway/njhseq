/*
 * BedRecord.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "BedRecord.hpp"

namespace bibseq {
BedRecord::BedRecord(const std::string & line){
	auto toks = tokenizeString(line, "\t");
	if(toks.size() < 6){
		std::stringstream ss;
		ss << "Error in parsing line: " << line << "\n";
		ss << "Error in read bed file, need to have at least first 6 fields not: " << toks.size() << std::endl;
		ss << "1)chrom,2)chromStart,3)chromEnd,4)name,5)score,6)strand" << std::endl;
		throw std::runtime_error{ss.str()};
	}
	chrom_ = toks[0];
	chromStart_ = bib::lexical_cast<uint32_t>(toks[1]);
	chromEnd_ = bib::lexical_cast<uint32_t>(toks[2]);
	name_ = toks[3];
	score_ = bib::lexical_cast<double>(toks[4]);
	if(toks[5] == "-"){
		strand_ = '-';
	}else if(toks[5] == "+"){
		strand_ = '+';
	}else{
		std::stringstream ss;
		ss << "Error in parsing strand field: " << toks[5] << "\n";
		ss << "should be either + or - not " << toks[5] << "\n";
		throw std::runtime_error{ss.str()};
	}
}

BedRecord::BedRecord():chrom_(""), chromStart_(0), chromEnd_(0), name_(""), score_(0), strand_('+'){

}

bool BedRecord::reverseStrand()const{
	return '-' == strand_;
}

std::string BedRecord::toDelimStr()const{
	return vectorToString(toVecStr(chrom_, chromStart_, chromEnd_, name_, score_, strand_), "\t");
}
} /* namespace bibseq */
