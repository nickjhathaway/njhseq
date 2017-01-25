/*
 * GeneSeqInfo.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: nick
 */


#include "GeneSeqInfo.hpp"


namespace bibseq {

GeneSeqInfo::GeneSeqInfo(const std::string & name, const seqInfo & cDna,
		const seqInfo & gDna, uint32_t genomicStart, bool oneBasedPos,
		aligner & alignerObj) :
		name_(name), cDna_(cDna), gDna_(gDna), genomicStart_(genomicStart), oneBasedPos_(
				oneBasedPos), infoTab_(VecStr { "cDnaPos", "gDnaPos", "aaPos",
				"codonPos", "base", "aminoAcid" }) {

	alignerObj.alignCacheGlobal(gDna, cDna);
	//check to see if alignment went ok
	for (auto pos : iter::range(len(alignerObj.alignObjectA_))) {
		if (alignerObj.alignObjectA_.seqBase_.seq_[pos] != '-'
				&& alignerObj.alignObjectB_.seqBase_.seq_[pos] != '-') {
			if (alignerObj.alignObjectA_.seqBase_.seq_[pos]
					!= alignerObj.alignObjectB_.seqBase_.seq_[pos]) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error at aln position " << pos
						<< " cDNA doesn't match gDNA, "
								"check alignment parameters or make introns lowercase if they aren't already to aid in proper alignment"
						<< "\n";
				throw std::runtime_error { ss.str() };
			}
		}
	}
	cDnaAln_ = alignerObj.alignObjectB_.seqBase_;
	auto protein = cDna.translateRet(false, false);
	uint32_t posOffset = 0;
	if (oneBasedPos) {
		posOffset = 1;
	}
	for (auto pos : iter::range(len(cDna))) {
		infoTab_.addRow(pos + posOffset,
				alignerObj.getSeqPosForAlnAPos(alignerObj.getAlignPosForSeqBPos(pos))
						+ genomicStart, (pos / 3) + posOffset, (pos % 3) + posOffset,
				cDna.seq_[pos],
				len(protein.seq_) > pos / 3 ?
						std::string(1, protein.seq_[pos / 3]) : "NA");
	}

}




}  // namespace bibseq

