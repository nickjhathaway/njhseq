#pragma once
/*
 * FoldbackDetector.hpp
 *
 *  Created on: Nov 5, 2026
 *      Author: nicholashathaway
 */
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include "njhseq/utils.h"

namespace njhseq {

/**@brief Find foldbacks (inverse chimeras) in a read using only the read, no reference or target needed
 *
 * In a foldback the read carries a stretch and then the reverse complement of that stretch. A kmer at position i whose
 * reverse complement is at position j votes for i + j, and all the kmers of two paired arms vote for the same sum
 * whatever is between the arms, so the strongest sum is the fold. The anchors of that sum are chained out along the
 * arms, as insertions and deletions make the sum drift over a long arm, and the arms' ends are then found by extending
 * from the last anchors with errors allowed, since a lost kmer would otherwise put an arm's end too early.
 *
 * What this gives is where the arms are on the read. The stretch between them is the loop, the part of the template
 * that the fold reached back over. Comparing the read to itself can't say which side of the loop the fold happened on,
 * that needs where the arms are on a reference.
 *
 * Kmers that are their own reverse complement can't vote, and two kmers can't vote if they overlap, so an (AT)n tract
 * can't pair with itself to make a call, the arms on either side of it have to pair. The extension can still run into
 * such a tract, which is how a fold at a self complementary tract ends up with no loop.
 *
 * Nothing here tells a foldback made in a PCR from two finished products ligated head to head, or from an inverted
 * repeat that's really in the template, like two Alus facing each other. Those need where the read is on a reference, an
 * inverted repeat in the template maps as one contiguous alignment and a foldback as two on opposite strands.
 *
 * The read should have its adapters and barcodes trimmed. The same adapter on both ends of a molecule is an inverted
 * repeat, and when what's between them is short enough it gets called as a fold with short arms.
 */
class FoldbackDetector{
public:
	struct Pars{
		uint32_t kmerLen_{13}; /**< no more than 31, survives the error between two arms often enough while being unlikely by chance within a read */
		uint32_t seedBand_{30}; /**< the width of the range of sums the strongest fold is looked for within */
		uint32_t maxKmerOccurrences_{4}; /**< kmers seen more than this in the read don't vote, stops repeats from swamping the votes */
		uint32_t minArmLen_{60}; /**< the shorter arm has to be at least this long */
		uint32_t minAnchors_{5}; /**< the fewest kmers that have to pair */
		double minIdentity_{0.70}; /**< the lowest identity between the arms, two arms with 10% error each are around 0.81 */
		double minAnchorsPerBase_{0.03}; /**< the fewest kmers that have to pair per base of the anchored part of the arm, two arms with 10% error each are around 0.06 */
		uint32_t maxChainGap_{150}; /**< the furthest apart two anchors can be and still be on the same arms */
		double maxLoopPerArm_{3}; /**< the loop can be at most this times the longer arm plus maxLoopFlank_, stops sequence shared by the far ends of a read, e.g. adapters, from being a fold */
		uint32_t maxLoopFlank_{100};
		uint32_t extensionBand_{16}; /**< how far apart the two arms can get from insertions and deletions during extension */
		int32_t extensionXDrop_{12}; /**< stop extending when the score has fallen this far below its best */
	};

	struct Fold{
		bool called_{false}; /**< met all the requirements in Pars, the other values are of the best candidate even if false */
		uint32_t anchors_{0}; /**< kmers paired between the two arms */
		uint32_t anchoredLen_{0}; /**< how much of the left arm is between its first and last anchor */
		uint32_t otherVotes_{0}; /**< reverse complement kmer pairs in the read that aren't part of this fold, a high number means more than one fold */
		uint32_t leftStart_{0}; /**< positions on the read, 0 based, starts inclusive and ends exclusive */
		uint32_t leftEnd_{0};
		uint32_t rightStart_{0};
		uint32_t rightEnd_{0};
		double identity_{0}; /**< between the arms over the anchored stretch, a value for how well they pair and not of the reads' quality, two errors can be on the same position */

		[[nodiscard]] uint32_t leftArmLen() const;
		[[nodiscard]] uint32_t rightArmLen() const;
		[[nodiscard]] uint32_t shorterArmLen() const;
		[[nodiscard]] uint32_t loopLen() const;
		[[nodiscard]] double anchorsPerBase() const;
		[[nodiscard]] bool hasCandidate() const;
	};

	FoldbackDetector();
	explicit FoldbackDetector(Pars pars);

	Pars pars_;

	/**@brief Look for a foldback, safe to call from multiple threads
	 *
	 * @param seq the read, either case, anything other than ACGT doesn't vote
	 * @return the best candidate fold, Fold::hasCandidate() will be false if there were no reverse complement kmers at all
	 */
	[[nodiscard]] Fold detect(const std::string & seq) const;
};

}  // namespace njhseq
