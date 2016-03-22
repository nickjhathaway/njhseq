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
/*
 * BamToolsUtils.cpp
 *
 *  Created on: May 28, 2015
 *      Author: nick
 */

#include "BamToolsUtils.hpp"

#include <api/BamReader.h>

namespace bibseq {

std::unordered_map<size_t,alnInfoLocal> bamAlnToAlnInfoLocal(const BamTools::BamAlignment & bAln){
	std::unordered_map<size_t,alnInfoLocal> ret;
	//check for N cigar data
	bool hasN = false;
	for(const auto & cigar : bAln.CigarData){
		if(cigar.Type == 'N'){
			hasN = true;
			break;
		}
	}
	if(hasN){
		size_t refOffSet = 0;
		size_t seqOffSet = 0;
		size_t currentPosition = bAln.Position;
		std::vector<BamTools::CigarOp> currentOps;
		for(const auto & cData : bAln.CigarData){
			switch (cData.Type) {
				case 'M':
					refOffSet += cData.Length;
					seqOffSet += cData.Length;
					currentOps.emplace_back(cData);
					break;
				case 'D':
					refOffSet += cData.Length;
					currentOps.emplace_back(cData);
					break;
				case 'N':
					ret[currentPosition] = bamAlnToAlnInfoLocal(currentOps);
					refOffSet += cData.Length;
					currentOps.clear();
					currentOps.emplace_back(BamTools::CigarOp('S', seqOffSet));
					currentPosition = bAln.Position + refOffSet;
					break;
				case 'I':
					//do nothing
					seqOffSet += cData.Length;
					currentOps.emplace_back(cData);
					break;
				case 'S':
					//do nothing
					currentOps.emplace_back(cData);
					seqOffSet += cData.Length;
					break;
				case 'H':
					//do nothing
					currentOps.emplace_back(cData);
					//seqOffSet += cData.Length;
					break;
				default:
					std::cerr << __PRETTY_FUNCTION__ << ": Unhandled case: " << cData.Type << std::endl;
					break;
			}
		}
		if(currentOps.empty()){
			std::stringstream ss;
			ss << "Error in " << __PRETTY_FUNCTION__ << " in processing cigar operations, the options were empty in the loop" << std::endl;
			throw std::runtime_error{ss.str()};
		}
		ret[currentPosition] = bamAlnToAlnInfoLocal(currentOps);
	}else{
		ret[bAln.Position] = bamAlnToAlnInfoLocal(bAln.CigarData);
	}
	return ret;
}

alnInfoLocal bamAlnToAlnInfoLocal(const std::vector<BamTools::CigarOp> & cigarData){
	alnInfoLocal ret;
	ret.localASize_ = 0;
	ret.localAStart_ = 0;
	ret.addFromFile_ = false;
	ret.score_ = std::numeric_limits<double>::max();
	//first check for soft clips
	bool softClipFront = false;
	bool softClipBack = false;
	bool hardClipFront = false;
	bool hardClipBack = false;
	if (cigarData.front().Type == 'H') {
		hardClipFront = true;
		if (cigarData.size() > 1 && (cigarData.begin() + 1)->Type == 'S') {
			softClipFront = true;
			ret.localBStart_ = cigarData.front().Length;
		}
	} else if (cigarData.front().Type == 'S') {
		softClipFront = true;
		ret.localBStart_ = cigarData.front().Length;
	}
	if (cigarData.back().Type == 'H') {
		hardClipBack = true;
		if (cigarData.size() > 1 && (cigarData.end() - 2)->Type == 'S') {
			softClipBack = true;
		}
	} else if (cigarData.back().Type == 'S') {
		softClipBack = true;
	}

	for(const auto & cData : cigarData){
		switch (cData.Type) {
			case 'M':
				ret.localASize_ += cData.Length;
				ret.localBSize_ += cData.Length;
				break;
			case 'I':
				ret.localBSize_ += cData.Length;
				//gap position is position in the unclipped string so the start position needs to be added as well
				ret.gapInfos_.emplace_back(gapInfo(ret.localASize_ + ret.localAStart_,cData.Length, true));
				break;
			case 'D':
				ret.localASize_ += cData.Length;
				//gap position is position in the unclipped string so the start position needs to be added as well
				ret.gapInfos_.emplace_back(gapInfo(ret.localBSize_ + ret.localBStart_,cData.Length, false));
				break;
			case 'S':
				//do nothing
				break;
			case 'H':
				//do nothing
				break;
			default:
				std::cerr << __PRETTY_FUNCTION__ << ": Unhandled case: " << cData.Type << std::endl;
				break;
		}
	}
	std::sort(ret.gapInfos_.begin(), ret.gapInfos_.end(),[](const gapInfo & gap1, const gapInfo & gap2){ return gap1 > gap2;});
	return ret;
}


seqInfo bamAlnToSeqInfo(const BamTools::BamAlignment & aln, bool complement) {
	auto ret = seqInfo(aln.Name, aln.QueryBases, aln.Qualities,
			SangerQualOffset);
	if (complement) {
		ret.reverseComplementRead(false, false);
	}
	return ret;
}


void bamWriteFasta(const BamTools::BamAlignment & aln,
		std::ostream & out, bool complement){
	if(complement && aln.IsReverseStrand()){
		auto tempInfo = bamAlnToSeqInfo(aln, complement);
		tempInfo.outPutFastq(out);
	}else{
		out << ">" << aln.Name << "\n";
		out << aln.QueryBases << "\n";
	}
}

void bamWriteFastq(const BamTools::BamAlignment & aln,
		std::ostream & out, bool complement){
	if(complement && aln.IsReverseStrand()){
		auto tempInfo = bamAlnToSeqInfo(aln, complement);
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
