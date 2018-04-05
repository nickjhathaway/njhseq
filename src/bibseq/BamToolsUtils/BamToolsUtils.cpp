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

#include <api/BamReader.h>
#include "BamToolsUtils.hpp"

namespace bibseq {


std::vector<GenomicRegion> genGenRegionsFromRefData(const BamTools::RefVector & rData){
	std::vector<GenomicRegion> ret;
	for(const auto & r : rData){
		ret.emplace_back(r.RefName, r.RefName, 0, r.RefLength, false);
	}
	return ret;
}

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
	//bool softClipFront = false;
	//bool softClipBack = false;
	//bool hardClipFront = false;
	//bool hardClipBack = false;
	if (cigarData.front().Type == 'H') {
		//hardClipFront = true;
		if (cigarData.size() > 1 && (cigarData.begin() + 1)->Type == 'S') {
			//softClipFront = true;
			ret.localBStart_ = cigarData.front().Length;
		}
	} else if (cigarData.front().Type == 'S') {
		//softClipFront = true;
		ret.localBStart_ = cigarData.front().Length;
	}
	/*
	if(cigarData.back().Type == 'H') {
		hardClipBack = true;
		if (cigarData.size() > 1 && (cigarData.end() - 2)->Type == 'S') {
			softClipBack = true;
		}
	} else if (cigarData.back().Type == 'S') {
		softClipBack = true;
	}*/

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


seqInfo bamAlnToSeqInfo(const BamTools::BamAlignment & aln, bool keepPlusStrandOrientation) {
	auto ret = seqInfo(aln.Name, aln.QueryBases, aln.Qualities,
			SangerQualOffset);
	if (aln.IsReverseStrand() && !keepPlusStrandOrientation) {
		ret.reverseComplementRead(false, true);
	}
	return ret;
}


void bamWriteFasta(const BamTools::BamAlignment & aln,
		std::ostream & out){
	if(aln.IsReverseStrand()){
		auto tempInfo = bamAlnToSeqInfo(aln);
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
		tempInfo.outPutFastq(out);
	}else{
		out << "@" << aln.Name << "\n";
		out << aln.QueryBases << "\n";
		out << "+\n";
		out << aln.Qualities << "\n";
	}
}


bool IsBamSorted(const std::string & filename, bool verbose) {
	BamTools::BamReader bReader;
	bReader.Open(filename);
	BamTools::BamAlignment aln;
	BamTools::BamAlignment previousAln;
	std::vector<int32_t> alreadySeenRefIds;
	bReader.GetNextAlignmentCore(previousAln);
	uint32_t count = 1;
	while (bReader.GetNextAlignmentCore(aln)) {
		if(verbose && count % 1000 == 0){
			std::cout << '\r' << count;
			std::cout.flush();
		}
		if (previousAln.RefID == aln.RefID) {
			if (previousAln.Position > aln.Position) {
				return false;
			}
		} else {
			alreadySeenRefIds.emplace_back(previousAln.RefID);
			if (bib::in(aln.RefID, alreadySeenRefIds)) {
				return false;
			}
		}
		previousAln = aln;
		++count;
	}
	if(verbose){
		std::cout << std::endl;
	}
	return true;
}


void checkBamOpenThrow(BamTools::BamReader & bReader,
		const bfs::path & bamFnp) {
	if (!bReader.IsOpen()) {
		std::stringstream ss;
		if (!bfs::exists(bamFnp)) {
			ss << "Error " << bamFnp << " doesn't exist" << std::endl;
		} else {
			ss << "Error in opening " << bamFnp << std::endl;
		}
		throw std::runtime_error { ss.str() };
	}
}

void loadBamIndexThrow(BamTools::BamReader & bReader){
	if(!bReader.LocateIndex()){
		std::stringstream ss;
		ss << "Error: can't find index for " << bReader.GetFilename() << std::endl;
		ss << "Should be " << bReader.GetFilename() << ".bai" << std::endl;
		throw std::runtime_error{ss.str()};
	}
}


std::vector<BamTools::RefData> tabToRefDataVec(const table & refDataTab){
	std::vector<BamTools::RefData> ret;
	VecStr neededCols{"refId", "refName", "refLength"};
	VecStr missing;
	for(const auto & col : neededCols){
		if(!bib::in(col, refDataTab.columnNames_)){
			missing.emplace_back(col);
		}
	}

	if (!missing.empty()) {
		std::stringstream ss;
		ss << "In: " << __PRETTY_FUNCTION__ << "Missing the following columns: "
				<< bib::bashCT::boldRed(vectorToString(missing, ",")) << std::endl;
		throw std::runtime_error { ss.str() };
	}
	uint32_t id = 0;
	for (const auto & row : refDataTab.content_) {
		auto currentRefId = estd::stou(row[refDataTab.getColPos("refId")]);
		if (currentRefId != id) {
			std::stringstream ss;
			ss << "In: " << bib::bashCT::boldRed(__PRETTY_FUNCTION__)
					<< ", RefIds are not in order, expected " << id << " but found "
					<< currentRefId << " indstead" << std::endl;
			throw std::runtime_error { ss.str() };
		}
		++id;
		ret.emplace_back(row[refDataTab.getColPos("refName")],
				bib::lexical_cast<int32_t>(row[refDataTab.getColPos("refLength")]));
	}
	return ret;
}

table refDataVecToTab(const std::vector<BamTools::RefData> & refInfos){
	table ret(VecStr{"refId", "refName", "refLength"});
	for(const auto pos : iter::range(refInfos.size())){
		ret.addRow(pos, refInfos[pos].RefName, refInfos[pos].RefLength);
	}
	return ret;
}

void logAlnInfo(std::ostream & out, BamTools::RefVector & refInfo,
		const BamTools::BamAlignment & aln, std::string indent ) {

	out << indent << "aln.Name:                       " << aln.Name << std::endl;
	out << indent << "aln.MapQuality:                 " << aln.MapQuality<< std::endl;
	out << indent << "aln.IsFirstMate():              " << bib::colorBool(aln.IsFirstMate()) << std::endl;
	out << indent << "aln.Position:                   " << aln.Position << std::endl;
	out << indent << "refInfo[aln.RefID].RefName:     " << refInfo[aln.RefID].RefName << std::endl;
	out << indent << "aln.GetEndPosition():           " << aln.GetEndPosition() << std::endl;
	out << indent << "aln.IsReverseStrand():          " << bib::colorBool(aln.IsReverseStrand()) << std::endl;
	out << indent << "aln.AlignedBases.size():        " << aln.AlignedBases.size() << std::endl;
	out << indent << "aln.QueryBases.size():          " << aln.QueryBases.size() << std::endl;
	out << indent << "aln.MatePosition:               " << aln.MatePosition << std::endl;
	out << indent << "refInfo[aln.MateRefID].RefName: " << refInfo[aln.MateRefID].RefName << std::endl;
	out << indent << "aln.IsMateReverseStrand():      " << bib::colorBool(aln.IsMateReverseStrand()) << std::endl;
	out << indent << "aln.InsertSize:                 " << aln.InsertSize << std::endl << std::endl;

}

void setBamFileRegionThrow(BamTools::BamReader & bReader, const GenomicRegion & region){
	auto refData = bReader.GetReferenceData();
	std::unordered_map<std::string, uint32_t> refNameToId;
	for (auto pos : iter::range(refData.size())) {
		refNameToId[refData[pos].RefName] = pos;
	}
	if(!bib::in(region.chrom_,refNameToId)){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
		ss << "Region: " << region.chrom_ << " not found in bam file" << std::endl;
		ss << "Options are: " << bib::conToStr(getVectorOfMapKeys(refNameToId), ",") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (!bReader.SetRegion(refNameToId.at(region.chrom_), region.start_,
			refNameToId.at(region.chrom_), region.end_)) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
		ss << bReader.GetErrorString() << std::endl;
		ss << "Failed to set region" << std::endl;
		ss << "Region: " << refNameToId.at(region.chrom_) << std::endl;
		ss << "Start: " << region.start_ << std::endl;
		ss << "Stop: " << region.end_ << std::endl;
		throw std::runtime_error { ss.str() };
	}
}

void checkBamFilesForIndexesAndAbilityToOpen(const std::vector<bfs::path> & bamFnps){
	bool fail = false;
	std::stringstream ss;
	for(const auto & bamFnp : bamFnps){
		if(!bfs::exists(bamFnp)){
			ss << "Error " << bamFnp  << " doesn't exist"<< "\n";
			fail = true;
		}else{
			BamTools::BamReader bReader;
			bReader.Open(bamFnp.string());
			if (!bReader.IsOpen()) {
				ss << "Error in opening " << bamFnp << "\n";
				fail = true;
			}else if(!bReader.LocateIndex()){
				ss << "Error: can't find index for " << bReader.GetFilename() << "\n";
				ss << "Should be " << bReader.GetFilename() << ".bai" << "\n";
				fail = true;
			}
		}
	}
	if(fail){
		throw std::runtime_error{ss.str()};
	}
}

} /* namespace bibseq */