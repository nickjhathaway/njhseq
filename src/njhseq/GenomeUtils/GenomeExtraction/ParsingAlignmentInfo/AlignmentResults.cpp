/*
 * AlignmentResults.cpp
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
#include "AlignmentResults.hpp"
#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"


namespace njhseq {

AlignmentResults::AlignmentResults(const BamTools::BamAlignment & bAln,
		const BamTools::RefVector & refData,
		bool keepPlusStrandOrientation) :
		bAln_(bAln),
		gRegion_(bAln, refData),
		alnSeq_(std::make_shared<seqInfo>(bamAlnToSeqInfo(bAln, keepPlusStrandOrientation))),
		keepPlusStrandOrientation_(keepPlusStrandOrientation){
	if(keepPlusStrandOrientation){
		gRegion_.reverseSrand_ = false;
	}
	gRegion_.setUidWtihCoordsStrand();
}

void AlignmentResults::setRefSeq(TwoBit::TwoBitFile & twobitReader) {
	if ("*" == gRegion_.chrom_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error setting ref seq, alignment "
				<< bAln_.Name << " didn't map" << "\n";
		throw std::runtime_error { ss.str() };
	}
	refSeq_ = std::make_shared<seqInfo>(gRegion_.extractSeq(twobitReader));
	//make seq is upper case, other wise lower case masked sequence will appear as mismatches
	stringToUpper(refSeq_->seq_);
}


void AlignmentResults::setRefSeq(const seqInfo & refSeq){
	refSeq_ = std::make_shared<seqInfo>(refSeq);
	//make seq is upper case, other wise lower case masked sequence will appear as mismatches
	stringToUpper(refSeq_->seq_);
}

void AlignmentResults::setComparison(bool keepAlignedObjects){
	if(nullptr == refSeq_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error setting comparison, refSeq hasn't been added yet"<< "\n";
		throw std::runtime_error { ss.str() };
	}
	aligner alignerObj(10, gapScoringParameters(5,1), substituteMatrix::createDegenScoreMatrix(1,-1));
	alignerObj.weighHomopolymers_ = false;
	alignerObj.countEndGaps_ = false;
	auto alnInfo = bamAlnToAlnInfoLocal(bAln_);


	alignerObj.alignObjectA_.seqBase_ = *refSeq_;
	alignerObj.alignObjectB_.seqBase_ = *alnSeq_;
	alignerObj.parts_.lHolder_ = alnInfo.begin()->second;
	if(bAln_.IsReverseStrand() && !keepPlusStrandOrientation_){
		alignerObj.alignObjectA_.seqBase_.reverseComplementRead(false, true);
		alignerObj.alignObjectB_.seqBase_.reverseComplementRead(false, true);
	}
	alignerObj.rearrangeObjsLocal(alignerObj.alignObjectA_, alignerObj.alignObjectB_);
	if(bAln_.IsReverseStrand() && !keepPlusStrandOrientation_){
		alignerObj.alignObjectA_.seqBase_.reverseComplementRead(false, true);
		alignerObj.alignObjectB_.seqBase_.reverseComplementRead(false, true);
	}
	alignerObj.scoreAlignment(false);
	if(keepAlignedObjects){
		refSeqAligned_ = std::make_shared<seqInfo>(alignerObj.alignObjectA_.seqBase_);
		alnSeqAligned_ = std::make_shared<seqInfo>(alignerObj.alignObjectB_.seqBase_);
	}
	alignerObj.profileAlignment(refSeq_, alnSeq_, false, true, false);
	for(auto & m : alignerObj.comp_.distances_.mismatches_){
		if(!gRegion_.reverseSrand_){
			m.second.refBasePos = gRegion_.start_ + m.second.refBasePos;
		}else{
			m.second.refBasePos = gRegion_.end_ - 1 - m.second.refBasePos;
		}
	}
	for(auto & g : alignerObj.comp_.distances_.alignmentGaps_){
		if(!gRegion_.reverseSrand_){
			g.second.refPos_ = gRegion_.start_ + g.second.refPos_;
		}else{
			g.second.refPos_ = gRegion_.end_ - 1 - g.second.refPos_;
		}
	}
	comp_ = alignerObj.comp_;
}

void AlignmentResults::setAlignedObjects(){
	if(nullptr == refSeq_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error setting comparison, refSeq hasn't been added yet"<< "\n";
		throw std::runtime_error { ss.str() };
	}
	aligner alignerObj(500, gapScoringParameters(5,1), substituteMatrix(2,-2));
	alignerObj.weighHomopolymers_ = false;
	alignerObj.countEndGaps_ = false;
	auto alnInfo = bamAlnToAlnInfoLocal(bAln_);
	alignerObj.alignObjectA_.seqBase_ = *refSeq_;
	alignerObj.alignObjectB_.seqBase_ = *alnSeq_;
	alignerObj.parts_.lHolder_ = alnInfo.begin()->second;
	if(bAln_.IsReverseStrand() && !keepPlusStrandOrientation_){
		alignerObj.alignObjectA_.seqBase_.reverseComplementRead(false, true);
		alignerObj.alignObjectB_.seqBase_.reverseComplementRead(false, true);
	}
	alignerObj.rearrangeObjsLocal(alignerObj.alignObjectA_, alignerObj.alignObjectB_);
	if(bAln_.IsReverseStrand() && !keepPlusStrandOrientation_){
		alignerObj.alignObjectA_.seqBase_.reverseComplementRead(false, true);
		alignerObj.alignObjectB_.seqBase_.reverseComplementRead(false, true);
	}
	refSeqAligned_ = std::make_shared<seqInfo>(alignerObj.alignObjectA_.seqBase_);
	alnSeqAligned_ = std::make_shared<seqInfo>(alignerObj.alignObjectB_.seqBase_);
}

char AlignmentResults::getAlignedBase(const GenomicRegion & region){
	if(region.chrom_ != gRegion_.chrom_ ||
			region.reverseSrand_ !=gRegion_.reverseSrand_  ||
			region.start_ < gRegion_.start_
			|| region.start_ >= gRegion_.end_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error region and determine region don't match up" << "\n";
		ss << "Query Region: " << njh::json::writeAsOneLine(region.toJson()) << '\n';
		ss << "Determined Region: " << njh::json::writeAsOneLine(gRegion_.toJson()) << '\n';
	}
	if(nullptr == refSeqAligned_ || nullptr == alnSeqAligned_){
		setAlignedObjects();
	}
	uint32_t normalizedPos = region.start_ - gRegion_.start_;
	return alnSeqAligned_->seq_[getAlnPosForRealPos(refSeqAligned_->seq_, normalizedPos)];

}


std::vector<std::shared_ptr<AlignmentResults>> gatherMapResults(
		const bfs::path & bamFnp, const bfs::path & twoBitFnp,
		const comparison & allowableErrors) {

	auto initialRes = gatherMapResults(bamFnp, twoBitFnp);
	std::vector<std::shared_ptr<AlignmentResults>> ret;

	for(const auto & res : initialRes){
		if(allowableErrors.passErrorProfile(res->comp_)){
			ret.emplace_back(res);
		}
	}


	return ret;
}

std::vector<std::shared_ptr<AlignmentResults>> gatherMapResults(
		const bfs::path & bamFnp, const bfs::path & twoBitFnp){
	BamTools::BamReader bReader;
	bReader.Open(bamFnp.string());
	checkBamOpenThrow(bReader, bamFnp.string());
	TwoBit::TwoBitFile twobitReader(twoBitFnp);
	auto twobitSeqNames = twobitReader.sequenceNames();
	auto bRefData = bReader.GetReferenceData();
	BamTools::BamAlignment bAln;
	std::vector<std::shared_ptr<AlignmentResults>> ret;
	while (bReader.GetNextAlignment(bAln)) {
		if (bAln.IsMapped()) {
			auto results = std::make_shared<AlignmentResults>(bAln, bRefData);
			results->setRefSeq(twobitReader);
			results->setComparison(false);
			ret.emplace_back(results);
		}
	}
	return ret;
}

std::vector<std::shared_ptr<AlignmentResults>> getUniqueLocationResults(
		std::vector<std::shared_ptr<AlignmentResults>> & alnResults) {
	if(alnResults.size() <= 1){
		return alnResults;
	}
	njh::sort(alnResults,
			[](const std::shared_ptr<AlignmentResults> & results1,
					const std::shared_ptr<AlignmentResults> & results2) {
				return results1->gRegion_.createUidFromCoords() < results2->gRegion_.createUidFromCoords();
			});
	std::vector<std::shared_ptr<AlignmentResults>>  ret;
	ret.emplace_back(alnResults.front());
	for(const auto pos : iter::range<uint32_t>(1, alnResults.size())){
		if(ret.back()->gRegion_.createUidFromCoords() != alnResults[pos]->gRegion_.createUidFromCoords()){
			ret.emplace_back(alnResults[pos]);
		}
	}
	return ret;
}

} /* namespace njhseq */
