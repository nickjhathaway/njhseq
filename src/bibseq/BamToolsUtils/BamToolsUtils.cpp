
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//
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
