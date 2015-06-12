#pragma once
/*
 * systemUtils.hpp
 *
 *  Created on: Jan 18, 2015
 *      Author: nickhathaway
 */

#include "bibseq/objects/seqObjects/readObject.hpp"
namespace bibseq{
/**@b put padding zeros into the quality to represent gaps in the seq
 *
 * @param info the sequence to readjust the quality for after being externally aligned
 */
void adjustAlnSeqsQual(seqInfo & info);
namespace sys{
/**@b run muscle on this file
 *
 * @param filename name of the fasta file
 * @return a vector of the algined seqs
 */
std::vector<readObject> muscleSeqs(const std::string & filename);


/**@b muscle the sequences in seqs, leave the original seqs alone
 *
 * @param seqs the sequence to align
 * @return a vector of the alinged seqs
 */
template<typename T>
std::vector<readObject> muscleSeqsRet(std::vector<T> seqs){
	muscleSeqs(seqs);
	return seqs;
}

/**@b muscle the sequences in seqs
 *
 * @param seqs the sequence to align
 */
template<typename T>
void muscleSeqs(std::vector<T> & seqs){
	//create temporary file, the last 6 xs will be randomized characters
	char *tmpname = strdup("/tmp/tmpfileXXXXXX");
	mkstemp(tmpname);
	{
		std::ofstream tFile(tmpname);
		if(!tFile){
			throw std::runtime_error{bib::bashCT::boldRed("Error in opening " + std::string(tmpname))};
		}
		//make name the read position as muscle will reorganize the seqs aftewards
		uint32_t pos = 0;
		for(const auto & read : seqs){
			if(!containsSubString(read.seqBase_.seq_, "*")){
				tFile << ">" << pos << "\n";
				tFile << read.seqBase_.seq_ << "\n";
			}
			++pos;
		}
	}
	std::vector<readObject> tempObjs;
	tempObjs.reserve(seqs.size());
	try {
		//std::cout << tmpname << std::endl;
		tempObjs = muscleSeqs(tmpname);
	} catch (std::exception & e) {
		//std::cerr << e.what() << std::endl;
		bib::bfs::remove(tmpname);
		throw;
	}
	bib::bfs::remove(tmpname);
	//replace the sequences with the aligned sequences and adjust the qual
	for(const auto & read : tempObjs){
		auto & currentRead = seqs[std::stoul(read.seqBase_.name_)];
		currentRead.seqBase_.seq_ = read.seqBase_.seq_;
		adjustAlnSeqsQual(currentRead.seqBase_);
	}
}



} //namepsace sys
} //namespace bibseq
