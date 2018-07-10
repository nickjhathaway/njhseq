#pragma once

/*
 * GeneSeqInfo.hpp
 *
 *  Created on: Jan 15, 2017
 *      Author: nick
 */


#include "bibseq/alignment/aligner/aligner.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"

#include "bibseq/objects/BioDataObject/GenomicRegion.hpp"

namespace bibseq {

/**@brief class with info on a gene
 *
 * @todo add utilty functions to access info like if queried for amino acid position to return cDna and genomic positions etc
 *
 */
class GeneSeqInfo {
public:
	struct GeneSeqInfoPars {

		GenomicRegion region_;
		bool oneBasedPos_ = false;
	};
	GeneSeqInfo(const seqInfo & cDna, const seqInfo & gDna, GeneSeqInfoPars pars);

	void setTable();

	void setCDnaAln(const seqInfo & cDnaAln);
	void setCDnaAlnByAligning(aligner & alignerObj);

	Bed6RecordCore genBedFromAAPositions(uint32_t aaStart, uint32_t aaStop);
	Bed6RecordCore genBedFromCDNAPositions(uint32_t start, uint32_t stop);

	seqInfo cDna_;
	seqInfo gDna_;
	std::shared_ptr<seqInfo> cDnaAln_;
	seqInfo protein_;
	GeneSeqInfoPars pars_;
	table infoTab_;

	struct GenePosInfo {
		uint32_t gDNAPos_ = std::numeric_limits<uint32_t>::max();
		uint32_t cDNAPos_ = std::numeric_limits<uint32_t>::max();
		uint32_t aaPos_ = std::numeric_limits<uint32_t>::max();
		uint32_t codonPos_ = std::numeric_limits<uint32_t>::max();
		char base_;
		char aa_;
	};

	std::unordered_map<uint32_t, GenePosInfo> getInfosByGDNAPos() const;
	std::unordered_map<uint32_t, GenePosInfo> getInfosByCDNAPos() const;
	std::unordered_map<uint32_t, std::tuple<GenePosInfo,GenePosInfo,GenePosInfo>> getInfosByAAPos() const;


};



}  // namespace bibseq
