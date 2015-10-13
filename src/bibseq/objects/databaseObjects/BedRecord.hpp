#pragma once
/*
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"

namespace bibseq {




class BedRecord {
public:

	BedRecord(const std::string & line);
	BedRecord();

	std::string chrom_;
	uint32_t chromStart_;
	uint32_t chromEnd_;
	std::string name_;
	double score_;
	char strand_; //either + or -

	bool reverseStrand()const;

	std::string toDelimStr()const;
};


} /* namespace bibseq */



