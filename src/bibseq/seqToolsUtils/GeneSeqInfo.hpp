#pragma once

/*
 * GeneSeqInfo.hpp
 *
 *  Created on: Jan 15, 2017
 *      Author: nick
 */


#include "bibseq/alignment/aligner/aligner.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"


namespace bibseq {


/**@brief class with info on a gene
 *
 * @todo add utilty functions to access info like if queried for amino acid position to return cDna and genomic positions etc
 *
 */
class GeneSeqInfo{
public:
	GeneSeqInfo(const std::string & name,
			const seqInfo & cDna,
			const seqInfo & gDna,
		uint32_t genomicStart, bool oneBasedPos, aligner & alignerObj);

	std::string name_;
	seqInfo cDna_;
	seqInfo gDna_;
	uint32_t genomicStart_;
	bool oneBasedPos_;

	seqInfo cDnaAln_;

	table infoTab_;


};



}  // namespace bibseq
