#pragma once

/*
 * tandemRepeatUtils.hpp
 *
 *  Created on: May 31, 2020
 *      Author: nick
 */



#include "njhseq/programUtils/seqSetUp.hpp"
#include "njhseq/objects/helperObjects/motif.hpp"
#include "njhseq/IO/SeqIO/SeqInput.hpp"
#include "njhseq/IO/OutputStream.hpp"
#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"




namespace njhseq {

struct SimpleTRFinderLocsPars{
	std::vector<char> alphabet{'A', 'G', 'C', 'T'};
	OutOptions outOpts{bfs::path(""), ".bed"};
	uint32_t minRepeatUnitSize = 1;
	uint32_t maxRepeatUnitSize = 3;
	uint32_t lengthCutOff = 12;
	uint32_t minNumRepeats = 2;
	bool searchAllUnits = false;
	bool doNotAddFlankingSeq = false;
	bool addFullSeqToOuput = false;
	uint32_t numThreads = 1;

	bool verbose = false;
	bool debug = false;

	void setDefaultOpts(seqSetUp & setUp);

};

void runSimpleTRFinderLocsPars(const SeqIOOptions & seqInput, const SimpleTRFinderLocsPars & pars);


}  // namespace njhseq
