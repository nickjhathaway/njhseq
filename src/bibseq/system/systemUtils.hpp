#pragma once
/*
 * systemUtils.hpp
 *
 *  Created on: Jan 18, 2015
 *      Author: nickhathaway
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "bibseq/objects/seqObjects/readObject.hpp"
namespace bibseq{
/**@brief put padding zeros into the quality to represent gaps in the seq
 *
 * @param info the sequence to readjust the quality for after being externally aligned
 */
void adjustAlnSeqsQual(seqInfo & info);
namespace sys{
/**@brief run muscle on this file
 *
 * @param filename name of the fasta file
 * @return a vector of the algined seqs
 */
std::vector<readObject> muscleSeqs(const std::string & filename);


/**@brief muscle the sequences in seqs, leave the original seqs alone
 *
 * @param seqs the sequence to align
 * @return a vector of the alinged seqs
 */
template<typename T>
std::vector<readObject> muscleSeqsRet(std::vector<T> seqs){
	muscleSeqs(seqs);
	return seqs;
}

/**@brief muscle the sequences in seqs
 *
 * @param seqs the sequence to align
 * @param selected only muscle these selected sequences at these positions, will throw if out of range
 */
template<typename T, typename POSTYPE>
void muscleSeqs(std::vector<T> & seqs, const std::vector<POSTYPE> & selected){
	//create temporary file, the last 6 xs will be randomized characters
	char *tmpname = strdup("/tmp/tmpfileXXXXXX");
	mkstemp(tmpname);
	uint32_t seqsWritten = 0;
	{
		std::ofstream tFile(tmpname);
		if(!tFile){
			throw std::runtime_error{bib::bashCT::boldRed("Error in opening " + std::string(tmpname))};
		}
		//make name the read position as muscle will reorganize the seqs afterwards

		for (const auto & pos : selected) {
			if (pos >= seqs.size()) {
				throw std::out_of_range {
						"Error in bibseq::sys::muscleSeqs, position out of range, pos: "
								+ estd::to_string(pos) + ", size: "
								+ estd::to_string(seqs.size()) };
			}
			//hack because muscle doesn't like stop codons
			if(!bib::containsSubString(seqs[pos].seqBase_.seq_, "*")){
				++seqsWritten;
				tFile << ">" << pos << "\n";
				tFile << seqs[pos].seqBase_.seq_ << "\n";
			}
		}
	}
	//this is for when there are no sequences written due to all having stop codons which muscle won't align;
	if(seqsWritten > 0){
		std::vector<readObject> tempObjs;
		try {
			tempObjs = muscleSeqs(tmpname);
		} catch (std::exception & e) {
			bib::files::bfs::remove(tmpname);
			throw e;
		}

		tempObjs = muscleSeqs(tmpname);
		bib::files::bfs::remove(tmpname);
		//replace the sequences with the aligned sequences and adjust the qual
		for(const auto & read : tempObjs){
			/**@todo to preserve case, create glovalGapInfo from string and rearrange out sequence*/
			auto & currentRead = seqs[std::stoul(read.seqBase_.name_)];
			currentRead.seqBase_.seq_ = read.seqBase_.seq_;
			adjustAlnSeqsQual(currentRead.seqBase_);
		}
	}
}

/**@brief muscle the sequences in seqs
 *
 * @param seqs the sequence to align
 */
template<typename T>
void muscleSeqs(std::vector<T> & seqs){
	std::vector<uint64_t> allSelected(seqs.size());
	bib::iota<uint64_t>(allSelected, 0);
	muscleSeqs(seqs, allSelected);
}





} //namepsace sys
} //namespace bibseq
