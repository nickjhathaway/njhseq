#pragma once
/*
 * AlignmentResults.hpp
 *
 *  Created on: Apr 29, 2017
 *      Author: nick
 */
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

#include "njhseq/seqToolsUtils/seqToolsUtils.hpp"
#include "njhseq/objects/BioDataObject/GenomicRegion.hpp"

namespace njhseq {


/**@brief wrapper class to hold bam alignment record and the sequences from both the reference and the alignment and comparison of the two
 *
 */
class AlignmentResults {
public:

	/**@brief construct with bam alignment record
	 *
	 * @param bAln the alignment record
	 * @param refData the ref data from the alignment file so chromosome name can be recovered
	 * @param keepPlusStrandOrientation whether to keep plus strand orientation or not
	 */
	AlignmentResults(const BamTools::BamAlignment & bAln,
			const BamTools::RefVector & refData, bool keepPlusStrandOrientation = false) ;

	/**@brief extract from a twobit reader the reference sequence from the corresponding alignment
	 *
	 * @param twobitReader
	 */
	void setRefSeq(TwoBit::TwoBitFile & twobitReader) ;

	void setRefSeq(const seqInfo & refSeq) ;

	/**@brief compare the reference sequence and the alignment sequence, AlignmentResults::setRefSeq must be set first
	 *
	 */
	void setComparison(bool keepAlignedObjects);

	/**@brief set the aligned objects
	 *
	 */
	void setAlignedObjects();

	BamTools::BamAlignment bAln_; /**< the alignment record  */
	GenomicRegion gRegion_; /**< the genomic region the seq aligned to*/
	std::shared_ptr<seqInfo> refSeq_; /**< a reference seq record for the genomic region seq is aligned to, must be set with AlignmentResults::setRefSeq*/
	std::shared_ptr<seqInfo> alnSeq_;/**< a seq record of the aligned bases, set on construction*/
	comparison comp_;/**< the comparison between the seq record and the aligned bases, AlignmentResults::refSeq_ must be set first before this is generated */

	std::shared_ptr<seqInfo> refSeqAligned_; /**< the alginment of the reference to this seq, including indels and the such */
	std::shared_ptr<seqInfo> alnSeqAligned_; /**< a seq record the alginment of the aligned seq, set on construction */

	bool keepPlusStrandOrientation_ = false;

	char getAlignedBase(const GenomicRegion & region);


};


/**@brief Get all records from bamFnp that pass the given comparison threshold
 *
 * should only be called on small bam files, no complex checking done, just if the record mapped
 *
 * @param bamFnp the bam file to extract from, can be unsorted and unidexed
 * @param twoBitFnp the twobit file that the sequence was aligned to so the ref seq can be acquired
 * @param allowableErrors the number of errors to allow
 * @return a vector of njhseq::AlignmentResults shared pointers that passed the threshold and that mapped
 */
std::vector<std::shared_ptr<AlignmentResults>> gatherMapResults(
		const bfs::path & bamFnp, const bfs::path & twoBitFnp,
		const comparison & allowableErrors);

/**@brief Get all records from bamFnp that at least mapped
 *
 * should only be called on small bam files, no complex checking done, just if the record mapped
 *
 * @param bamFnp the bam file to extract from, can be unsorted and unidexed
 * @param twoBitFnp the twobit file that the sequence was aligned to so the ref seq can be acquired
 * @return a vector of njhseq::AlignmentResults shared pointers that at least ampped
 */
std::vector<std::shared_ptr<AlignmentResults>> gatherMapResults(
		const bfs::path & bamFnp, const bfs::path & twoBitFnp);

/**@brief get only unique genomic locations, will sort the original vector as well by genomic location
 *
 * @param alnResults the alignment results to only get the unique records from
 * @return the unique genomic locations
 */
std::vector<std::shared_ptr<AlignmentResults>> getUniqueLocationResults(
		std::vector<std::shared_ptr<AlignmentResults>> & alnResults);


} /* namespace njhseq */

