/*
 * BamToolsUtils.cpp
 *
 *  Created on: May 28, 2015
 *      Author: nick
 */

#include "BamToolsUtils.hpp"
#include "bibseq/IO/readObjectIO.hpp"
#include <api/BamReader.h>

namespace bibseq {


seqInfo bamAlnToSeqInfo(const BamTools::BamAlignment & aln){
	return seqInfo(aln.Name, aln.QueryBases,
			aln.Qualities, readObjectIO::SangerQualOffset);
}


void bamWriteFasta(const BamTools::BamAlignment & aln,
		std::ostream & out){
	if(aln.IsReverseStrand()){
		auto tempInfo = bamAlnToSeqInfo(aln);
		tempInfo.reverseComplementRead(false);
		tempInfo.outPutFastq(out);
	}else{
		out << ">" << aln.Name << "\n";
		out << aln.QueryBases << "\n";
	}
}

void bamWriteFastq(const BamTools::BamAlignment & aln,
		std::ostream & out){
	if(aln.IsReverseStrand()){
		auto tempInfo = bamAlnToSeqInfo(aln);
		tempInfo.reverseComplementRead(false);
		tempInfo.outPutFastq(out);
	}else{
		out << "@" << aln.Name << "\n";
		out << aln.QueryBases << "\n";
		out << "+\n";
		out << aln.Qualities << "\n";
	}
}


bool IsBamSorted(const std::string & filename){
	BamTools::BamReader bReader;
	bReader.Open(filename);
	BamTools::BamAlignment aln;
	int32_t previousPos = 0;
	while(bReader.GetNextAlignmentCore(aln)){
		if(aln.IsMapped()){
			//std::cout << previousPos << std::endl;
			if(previousPos > aln.Position){
				return false;
			}
			previousPos = aln.Position;
		}
	}
	return true;
}

} /* namespace bibseq */
