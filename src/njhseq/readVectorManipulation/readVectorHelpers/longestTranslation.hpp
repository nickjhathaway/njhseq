#pragma once

/*
 * longestTranslation.hpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

namespace njhseq {

/**@brief The longest stretch of protein without a stop codon across all six reading frames of a sequence
 *
 */
struct LongestTranslation {
	struct Pars {
		bool removeTrailingStop_ = false; /**< don't add the stop codon that ended the protein*/
		bool forceMStart_ = false; /**< passed to seqInfo::translateRet()*/
	};

	seqInfo protein_;
	uint32_t frame_{0}; /**< the frame start, 0-2*/
	bool isRevComp_{false}; /**< found on the reverse complement*/
	uint32_t seqStart_{0}; /**< start in the (reverse complemented when isRevComp_) sequence*/
	uint32_t seqStop_{0}; /**< stop in the (reverse complemented when isRevComp_) sequence, includes the stop codon when it's added*/
	uint32_t seqLen_{0};
	uint32_t aaStart_{0}; /**< start in the frame's whole translation*/
	uint32_t aaStop_{0}; /**< stop in the frame's whole translation, includes the stop when it's added*/
	uint32_t frameProteinLen_{0}; /**< length of the frame's whole translation*/

	/**@brief Translate all six frames, break each on stops, and keep the longest piece, ties go to the first found (forward frames 0-2, then reverse complement frames 0-2)
	 *
	 * The stop codon that ended the piece is added unless Pars::removeTrailingStop_. This is elucidator guessAProteinFromSeq --getLongest.
	 * A sequence whose every frame is all stops gives an empty protein.
	 *
	 * @param seq the sequence
	 * @param pars stop and start handling
	 * @return the longest piece and where it came from
	 */
	[[nodiscard]] static LongestTranslation translate(const seqInfo & seq, const Pars & pars);

	/**@brief Add frame, seqStart, seqStop, seqLen, aaStart, aaStop, aaLen and isRevComp to the protein's name
	 *
	 * @param overWriteMeta overwrite the fields already in the name, except isRevComp
	 */
	void markName(bool overWriteMeta);
};

}  // namespace njhseq
