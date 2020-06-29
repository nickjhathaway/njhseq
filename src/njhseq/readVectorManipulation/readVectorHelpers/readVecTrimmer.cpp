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
#include "readVecTrimmer.hpp"
#include "njhseq/IO/SeqIO/SeqOutput.hpp"



namespace njhseq {


std::vector<seqInfo> readVecTrimmer::trimCircularGenomeToRef(seqInfo seq,
		const trimCircularGenomeToRefPars & pars,
		aligner & alignerObj){
	std::vector<seqInfo> outSeqs;

	if(!pars.doNotReOrientDirection_){
		kmerInfo refKInfo(pars.refSeq_.seq_, pars.kmerLength_, false);
		kmerInfo seqKInfo(seq.seq_, pars.kmerLength_, true);
		auto compForward = refKInfo.compareKmers(seqKInfo);
		auto compReverse = refKInfo.compareKmersRevComp(seqKInfo);
		if(compReverse.second > compForward.second){
			seq.reverseComplementRead(false, true);
		}
	}

	alignerObj.alignCacheGlobal(pars.refSeq_, seq);
	alignerObj.rearrangeObjsGlobal(pars.refSeq_, seq);
//	auto debugOutOpts = SeqIOOptions::genFastaOut("test.fasta");
//	debugOutOpts.out_.overWriteFile_ =true;
//	SeqOutput debugWriter(debugOutOpts);
//	debugWriter.openOut();
//	debugWriter.write(alignerObj.alignObjectA_);
//	debugWriter.write(alignerObj.alignObjectB_);
	if(
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' == alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' == alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()      ){
		//1
		//b starts in a with back overhang
		auto lastRefAlnPos =    alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
		auto firstAlnInputPos = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');

		auto refStartPosition = alignerObj.getSeqPosForAlnAPos(firstAlnInputPos);

		auto originalAlign1 = alignerObj.alignObjectA_;
		auto originalAlign2 = alignerObj.alignObjectB_;

		auto inputPos = alignerObj.getSeqPosForAlnBPos(lastRefAlnPos) + 1;
		auto backSeq = seq.getSubRead(inputPos);
		auto forwardSeqRef = pars.refSeq_.getSubRead(0, len(backSeq)+ pars.padding_);

		alignerObj.alignCacheGlobal(forwardSeqRef, backSeq);
		alignerObj.rearrangeObjsGlobal(forwardSeqRef, backSeq);
//		debugWriter.write(alignerObj.alignObjectA_);
//		debugWriter.write(alignerObj.alignObjectB_);

		auto lastRealignedPos = alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
		auto lastRefPositionForBackSeqInclusive = alignerObj.getSeqPosForAlnAPos(lastRealignedPos);
		auto lastRefPositionForBackSeqEnd = lastRefPositionForBackSeqInclusive + 1;

		uint32_t endPositionForBackSeqForOutput = std::numeric_limits<uint32_t>::max();
		uint32_t originalInputSeqStartForOutput = std::numeric_limits<uint32_t>::max();

		if(pars.preferHeader_){
			if(lastRefPositionForBackSeqEnd > refStartPosition){
				endPositionForBackSeqForOutput = alignerObj.getSeqPosForAlnBPos(alignerObj.getAlignPosForSeqAPos(refStartPosition));
			}
			originalInputSeqStartForOutput = 0;
		} else {
			if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
				//back seq went further than front of ref seq
				auto lastRefAlnPositionInRealalignment = alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
				endPositionForBackSeqForOutput = alignerObj.getSeqPosForAlnBPos(lastRefAlnPositionInRealalignment) + 1;
			}else{
				endPositionForBackSeqForOutput = len(backSeq);
				if(firstAlnInputPos <= lastRefPositionForBackSeqInclusive){
					originalInputSeqStartForOutput = getRealPosForAlnPos(originalAlign2.seqBase_.seq_, getAlnPosForRealPos(originalAlign1.seqBase_.seq_, lastRefPositionForBackSeqInclusive)) + 1;
				}else{
					originalInputSeqStartForOutput = 0;
				}
			}
		}
		auto startPositionForBackSeqForOutput = 0;
		if('-' != alignerObj.alignObjectA_.seqBase_.seq_.front() ){
			startPositionForBackSeqForOutput = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-'));
		}
		auto outSeq = backSeq.getSubRead(startPositionForBackSeqForOutput, endPositionForBackSeqForOutput - startPositionForBackSeqForOutput);
		if(lastRefPositionForBackSeqEnd >= refStartPosition){
			if(std::numeric_limits<uint32_t>::max() != originalInputSeqStartForOutput ){
				outSeq.append(seq.getSubRead(originalInputSeqStartForOutput, inputPos - originalInputSeqStartForOutput));
				if(pars.mark_){
					MetaDataInName seqMeta;
					if(MetaDataInName::nameHasMetaData(outSeq.name_)){
						seqMeta= MetaDataInName(outSeq.name_);
					}
					seqMeta.addMeta("length", len(outSeq), true);
					seqMeta.addMeta("inputLength", len(seq), true);
					seqMeta.addMeta("refTransitionPoint", refStartPosition);
					seqMeta.addMeta("inputPosition", njh::pasteAsStr(
							inputPos + startPositionForBackSeqForOutput, "-", inputPos + endPositionForBackSeqForOutput - startPositionForBackSeqForOutput,
							",",
							originalInputSeqStartForOutput, "-", inputPos - originalInputSeqStartForOutput));
					seqMeta.resetMetaInName(outSeq.name_);
				}
			}else{
				if(pars.mark_){
					MetaDataInName seqMeta;
					if(MetaDataInName::nameHasMetaData(outSeq.name_)){
						seqMeta= MetaDataInName(outSeq.name_);
					}
					seqMeta.addMeta("length", len(outSeq), true);
					seqMeta.addMeta("inputLength", len(seq), true);
					seqMeta.addMeta("refTransitionPoint", refStartPosition);
					seqMeta.addMeta("inputPosition", njh::pasteAsStr(
							inputPos + startPositionForBackSeqForOutput, "-", inputPos + endPositionForBackSeqForOutput - startPositionForBackSeqForOutput));
					seqMeta.resetMetaInName(outSeq.name_);
				}
			}
			outSeqs.emplace_back(outSeq);
		} else {
			auto backSeqTrimmed = seq.getSubRead(originalInputSeqStartForOutput, inputPos - originalInputSeqStartForOutput);
			outSeq.on_ = false;
			backSeqTrimmed.on_ = false;
			if(pars.mark_){
				{
					MetaDataInName seqMeta;
					if(MetaDataInName::nameHasMetaData(outSeq.name_)){
						seqMeta= MetaDataInName(outSeq.name_);
					}
					seqMeta.addMeta("length", len(outSeq), true);
					seqMeta.addMeta("inputLength", len(seq), true);
					seqMeta.addMeta("refTransitionPoint", refStartPosition);
					seqMeta.addMeta("inputPosition", njh::pasteAsStr(
							inputPos + startPositionForBackSeqForOutput, "-", inputPos + endPositionForBackSeqForOutput - startPositionForBackSeqForOutput));
					seqMeta.resetMetaInName(outSeq.name_);
				}
				{
					MetaDataInName seqMeta;
					if(MetaDataInName::nameHasMetaData(backSeqTrimmed.name_)){
						seqMeta= MetaDataInName(backSeqTrimmed.name_);
					}
					seqMeta.addMeta("length", len(backSeqTrimmed), true);
					seqMeta.addMeta("inputLength", len(seq), true);
					seqMeta.addMeta("refTransitionPoint", refStartPosition);
					seqMeta.addMeta("inputPosition", njh::pasteAsStr(originalInputSeqStartForOutput, "-", inputPos - originalInputSeqStartForOutput));
					seqMeta.resetMetaInName(backSeqTrimmed.name_);
				}
			}
			outSeqs.emplace_back(outSeq);
			outSeqs.emplace_back(backSeqTrimmed);
		}
	} else if (
			 '-' == alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' == alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//2
		//b ends in a with front overhang
		auto originalAlign1 = alignerObj.alignObjectA_;
		auto originalAlign2 = alignerObj.alignObjectB_;

		auto inputFrontEndPos = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-'));
		auto backRefPosRaw =    alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-')) + 1;
		auto backRefPosAdjusted = backRefPosRaw;

		if(backRefPosAdjusted > pars.padding_){
			backRefPosAdjusted -= pars.padding_;
		}else{
			backRefPosAdjusted = 0;
		}

		auto backRefPosRawInBack = backRefPosRaw - backRefPosAdjusted;

		auto refBackSeq = pars.refSeq_.getSubRead(backRefPosAdjusted);
		auto frontInputSeq = seq.getSubRead(0, inputFrontEndPos);

		alignerObj.alignCacheGlobal(refBackSeq, frontInputSeq);
		alignerObj.rearrangeObjsGlobal(refBackSeq, frontInputSeq);

		auto remainderInputSeq = seq.getSubRead(inputFrontEndPos);

		uint32_t positionOfFrontInputSeqInBackRef = 0;
		if('-' == alignerObj.alignObjectB_.seqBase_.seq_.front()){
			positionOfFrontInputSeqInBackRef = alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-'));
		}else{
			positionOfFrontInputSeqInBackRef = 0;
		}
		uint32_t positionOfFrontSeqInOriginalRef = positionOfFrontInputSeqInBackRef + backRefPosAdjusted;

		uint32_t remainderEndPosForOutput = std::numeric_limits<uint32_t>::max();
		uint32_t frontSeqStartPosForOutput = std::numeric_limits<uint32_t>::max();

		if(pars.preferHeader_){
			//front seq has overhang at the beginning
			if('-' != alignerObj.alignObjectB_.seqBase_.seq_.front()){
				frontSeqStartPosForOutput = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-'));
			}else{
				frontSeqStartPosForOutput = 0;
			}
			if(positionOfFrontSeqInOriginalRef < backRefPosRaw){
				remainderEndPosForOutput = getRealPosForAlnPos(originalAlign2.seqBase_.seq_, getAlnPosForRealPos(originalAlign1.seqBase_.seq_, positionOfFrontSeqInOriginalRef)) - inputFrontEndPos;
			}else{
				remainderEndPosForOutput = len(remainderInputSeq);
			}
		} else {
			if(backRefPosRawInBack > positionOfFrontInputSeqInBackRef){
				frontSeqStartPosForOutput = alignerObj.getSeqPosForAlnBPos(alignerObj.getAlignPosForSeqAPos(backRefPosRawInBack));
			}else{
				frontSeqStartPosForOutput = 0;
			}
			remainderEndPosForOutput = len(remainderInputSeq);
		}
		uint32_t frontSeqEndPosForOutput = len(frontInputSeq);
		if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back()){
			frontSeqEndPosForOutput = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-'));
		}

		auto outSeq = remainderInputSeq.getSubRead(0, remainderEndPosForOutput);
		if (positionOfFrontSeqInOriginalRef <= backRefPosRaw) {
			outSeq.append(frontInputSeq.getSubRead(frontSeqStartPosForOutput, frontSeqEndPosForOutput - frontSeqStartPosForOutput));
			if(pars.mark_){
				MetaDataInName seqMeta;
				if(MetaDataInName::nameHasMetaData(outSeq.name_)){
					seqMeta= MetaDataInName(outSeq.name_);
				}
				seqMeta.addMeta("length", len(outSeq), true);
				seqMeta.addMeta("inputLength", len(seq), true);
				seqMeta.addMeta("refTransitionPoint", backRefPosRaw);
				seqMeta.addMeta("inputPosition", njh::pasteAsStr(
						inputFrontEndPos, "-", inputFrontEndPos + remainderEndPosForOutput,
						",",
						frontSeqStartPosForOutput, "-", frontSeqEndPosForOutput));
				seqMeta.resetMetaInName(outSeq.name_);
			}
			outSeqs.emplace_back(outSeq);
		} else {
			auto trimmedFrontSeq = frontInputSeq.getSubRead(frontSeqStartPosForOutput, frontSeqEndPosForOutput - frontSeqStartPosForOutput);
			outSeq.on_ = false;
			trimmedFrontSeq.on_ = false;
			if(pars.mark_){
				{
					MetaDataInName seqMeta;
					if(MetaDataInName::nameHasMetaData(outSeq.name_)){
						seqMeta= MetaDataInName(outSeq.name_);
					}
					seqMeta.addMeta("length", len(outSeq), true);
					seqMeta.addMeta("inputLength", len(seq), true);
					seqMeta.addMeta("refTransitionPoint", backRefPosRaw);
					seqMeta.addMeta("inputPosition", njh::pasteAsStr(
							inputFrontEndPos, "-", inputFrontEndPos + remainderEndPosForOutput));
					seqMeta.resetMetaInName(outSeq.name_);
				}
				{
					MetaDataInName seqMeta;
					if(MetaDataInName::nameHasMetaData(trimmedFrontSeq.name_)){
						seqMeta= MetaDataInName(trimmedFrontSeq.name_);
					}
					seqMeta.addMeta("length", len(trimmedFrontSeq), true);
					seqMeta.addMeta("inputLength", len(seq), true);
					seqMeta.addMeta("refTransitionPoint", backRefPosRaw);
					seqMeta.addMeta("inputPosition", njh::pasteAsStr(frontSeqStartPosForOutput, "-", frontSeqEndPosForOutput));
					seqMeta.resetMetaInName(trimmedFrontSeq.name_);
				}
			}
			outSeqs.emplace_back(outSeq);
			outSeqs.emplace_back(trimmedFrontSeq);
		}
	} else if (
						 '-' == alignerObj.alignObjectA_.seqBase_.seq_.front() &&
				 	   '-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
						 '-' == alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
						 '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//3
		//b has overlap on front and back
		auto inputStart = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-'));
		auto inputEnd = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-')) + 1;
		auto outSeq  = seq.getSubRead(inputStart, inputEnd - inputStart);
		if(pars.mark_){
			MetaDataInName seqMeta;
			if(MetaDataInName::nameHasMetaData(outSeq.name_)){
				seqMeta = MetaDataInName(outSeq.name_);
			}
			seqMeta.addMeta("length", len(outSeq), true);
			seqMeta.addMeta("inputLength", len(seq), true);
			seqMeta.addMeta("inputPosition", njh::pasteAsStr(inputStart, "-", inputEnd - inputStart));
			seqMeta.resetMetaInName(outSeq.name_);
		}
		outSeqs.emplace_back(outSeq);
	}   else if (
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' == alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' == alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//4
		//b found completely within a
		auto outSeq = seq;
		outSeq.on_ = false;
		outSeqs.emplace_back(outSeq);
	}   else if (
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' == alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//5
		//starts together but a has overhang at end
		auto outSeq = seq;
		outSeq.on_ = false;
		outSeqs.emplace_back(outSeq);
	}   else if (
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' == alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//6
		//a has overhang at front but both end together
		auto outSeq = seq;
		outSeq.on_ = false;
		outSeqs.emplace_back(outSeq);
	}   else if (
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//7
		//perfect alignment nothing to do
		auto outSeq = seq;
		outSeq.on_ = false;
		outSeqs.emplace_back(outSeq);
	}   else if (
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' == alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//8
		//start together b has overhang at end
		auto inputStart = 0;
		auto inputEnd = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-')) + 1;
		auto outSeq = seq.getSubRead(inputStart, inputEnd - inputStart);
		if(pars.mark_){
			MetaDataInName seqMeta;
			if(MetaDataInName::nameHasMetaData(outSeq.name_)){
				seqMeta = MetaDataInName(outSeq.name_);
			}
			seqMeta.addMeta("length", len(outSeq), true);
			seqMeta.addMeta("inputLength", len(seq), true);
			seqMeta.addMeta("inputPosition", njh::pasteAsStr(inputStart, "-", inputEnd - inputStart));
			seqMeta.resetMetaInName(outSeq.name_);
		}
		outSeqs.emplace_back(outSeq);
	}     else if (
			 '-' == alignerObj.alignObjectA_.seqBase_.seq_.front() &&
	 	   '-' != alignerObj.alignObjectB_.seqBase_.seq_.front() &&
			 '-' != alignerObj.alignObjectA_.seqBase_.seq_.back()  &&
			 '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()     ) {
		//9
		//end together b has overhang at front
		auto inputStart = alignerObj.getSeqPosForAlnBPos(alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-'));
		auto inputEnd = len(seq);
		auto outSeq = seq.getSubRead(inputStart, inputEnd - inputStart);
		if(pars.mark_){
			MetaDataInName seqMeta;
			if(MetaDataInName::nameHasMetaData(outSeq.name_)){
				seqMeta = MetaDataInName(outSeq.name_);
			}
			seqMeta.addMeta("length", len(outSeq), true);
			seqMeta.addMeta("inputLength", len(seq), true);
			seqMeta.addMeta("inputPosition", njh::pasteAsStr(inputStart, "-", inputEnd - inputStart));
			seqMeta.resetMetaInName(outSeq.name_);
		}
		outSeqs.emplace_back(outSeq);
	}

	return outSeqs;
}


readVecTrimmer::TrimEdgesByLowEntropyRes::TrimEdgesByLowEntropyRes(uint32_t start, uint32_t end): start_(start), end_(end){

}


readVecTrimmer::TrimEdgesByLowEntropyRes readVecTrimmer::determineTrimPostionsByLowEntropy(const std::string & seq, const readVecTrimmer::TrimEdgesByLowEntropyPars & pars){
	readVecTrimmer::TrimEdgesByLowEntropyRes ret(0, seq.size());

	if (len(seq) >= pars.windowSize) {
		for (auto pos : iter::range<uint32_t>(0, len(seq) - pars.windowSize + 1,
				pars.windowStep)) {
			kmerInfo kInfo(seq.substr(pos, pars.windowSize), pars.kLen, false);
			if (kInfo.computeKmerEntropy() < pars.entropyCutOff) {
				ret.start_ = pos + pars.windowSize;
			} else {
				break;
			}
		}
		for (auto pos : iter::range<uint32_t>(0, len(seq) - pars.windowSize + 1,
				pars.windowStep)) {
			kmerInfo kInfo(
					seq.substr(seq.size() - pos - pars.windowSize, pars.windowSize),
					pars.kLen, false);
			if (kInfo.computeKmerEntropy() < pars.entropyCutOff) {
				ret.end_ = seq.size() - pos - pars.windowSize;
			} else {
				break;
			}
		}
	} else {
		kmerInfo kInfo(seq, pars.kLen, false);
		if (kInfo.computeKmerEntropy() < pars.entropyCutOff) {
			ret.end_ = 0;
		}
	}
	return ret;
}

readVecTrimmer::BreakUpRes::BreakUpRes(const seqInfo & seqBase, uint32_t start,
		uint32_t end, const std::string & pat) :
		seqBase_(seqBase), start_(start), end_(end), pat_(pat) {

}

std::vector<readVecTrimmer::BreakUpRes> readVecTrimmer::breakUpSeqOnPat(const seqInfo & seq,
		const std::string & pattern) {
	njh::PatPosFinder pFinder(pattern);
	return breakUpSeqOnPat(seq, pFinder);
}
std::vector<readVecTrimmer::BreakUpRes> readVecTrimmer::breakUpSeqOnPat(const seqInfo & seq,
		const njh::PatPosFinder & pFinder) {
	std::vector<BreakUpRes> ret;
	auto pats = pFinder.getPatPositions(seq.seq_);

	if (!pats.empty()) {
		if (0 != pats.front().pos_) {
			size_t start = 0;
			size_t end = pats.front().pos_;
			auto subSeq = seq.getSubRead(start, end - start);
			ret.emplace_back(subSeq, start, end, pats.front().pat_);
		}
		if (pats.size() > 1) {
			for (const auto & patPos : iter::range(pats.size() - 1)) {
				const auto & p = pats[patPos];
				size_t start = p.pos_ + p.pat_.size();
				size_t end = pats[patPos + 1].pos_;
				auto subSeq = seq.getSubRead(start, end - start);
				ret.emplace_back(subSeq, start, end, p.pat_);
			}
		}
		if (seq.seq_.size() != pats.back().end()) {
			size_t start = pats.back().end();
			size_t end = seq.seq_.size();
			auto subSeq = seq.getSubRead(start, end - start);
			ret.emplace_back(subSeq, start, end, pats.back().pat_);
		}
	} else {
		ret.emplace_back(seq, 0, len(seq), "");
	}
	return ret;
}


void readVecTrimmer::trimAtRstripQualScore(seqInfo &seq, const uint32_t qualCutOff){
	if(!seq.qual_.empty()){
		if(qualCutOff == seq.qual_.back()){
			uint32_t pos = seq.qual_.size() - 1;
			while(pos != 0 && qualCutOff == seq.qual_[pos - 1]){
				--pos;
			}
			seq.trimBack(pos);
		}
	}
}

void readVecTrimmer::trimAtRstripBase(seqInfo &seq, const char base){
	if(!seq.seq_.empty()){
		if(base == seq.seq_.back()){
			uint32_t pos = seq.seq_.size() - 1;
			while(pos != 0 && base == seq.seq_[pos - 1]){
				--pos;
			}
			seq.trimBack(pos);
		}
	}
}



void readVecTrimmer::trimAtLstripQualScore(seqInfo &seq, const uint32_t qualCutOff){
	if(!seq.qual_.empty()){
		if(qualCutOff == seq.qual_.front()){
			uint32_t pos = 0;
			while(pos + 1 != seq.qual_.size() && qualCutOff == seq.qual_[pos + 1]){
				++pos;
			}
			seq.trimFront(pos + 1);
		}
	}
}

void readVecTrimmer::trimAtLstripBase(seqInfo &seq, const char base){
	if(!seq.seq_.empty()){
		if(base == seq.seq_.front()){
			uint32_t pos = 0;
			while(pos + 1 != seq.seq_.size() && base == seq.seq_[pos + 1]){
				++pos;
			}
			seq.trimFront(pos + 1);
		}
	}
}




void readVecTrimmer::trimAtFirstBase(seqInfo &seq, const char base){
	auto pos = seq.seq_.find(base);
	if(std::string::npos != pos){
		seq.trimBack(pos);
	}
}

void readVecTrimmer::trimEdgesForLowEntropy(seqInfo &seq,
		const TrimEdgesByLowEntropyPars& pars) {
	auto positions = determineTrimPostionsByLowEntropy(seq.seq_, pars);
	if (positions.start_ < positions.end_) {
		seq = seq.getSubRead(positions.start_, positions.end_ - positions.start_);
		if (pars.mark) {
			MetaDataInName seqMeta;
			if (MetaDataInName::nameHasMetaData(seq.name_)) {
				seqMeta = MetaDataInName(seq.name_);
			}
			seqMeta.addMeta("trimStart", positions.start_, true);
			seqMeta.addMeta("trimEnd", positions.end_, true);
			seqMeta.resetMetaInName(seq.name_);
		}
		seq.on_ = true;
	} else {
		seq.on_ = false;
	}
}

void readVecTrimmer::trimToMaxLength(seqInfo &seq, size_t maxLength) {
	if (maxLength != 0 && len(seq) > maxLength - 1) {
		seq.trimBack(maxLength);
	} else if (len(seq) < maxLength) {
		seq.on_ = false;
	}
}

void readVecTrimmer::trimAtFirstQualScore(seqInfo &seq,
		const uint32_t qualCutOff) {
	auto iter = std::find_if(seq.qual_.begin(), seq.qual_.end(),
			[&qualCutOff](uint32_t qual) {return qual <= qualCutOff;});
	if (seq.qual_.end() != iter) {
		seq.trimBack(iter - seq.qual_.begin());
	}
}

void readVecTrimmer::trimToLastQualScore(seqInfo &seq,
		const uint32_t qualCutOff) {
	if(0 == len(seq)){
		return;
	}
	auto iter = std::find_if(seq.qual_.rbegin(), seq.qual_.rend(),
			[&qualCutOff](uint32_t qual) {return qual <= qualCutOff;});

	if (seq.qual_.rend() != iter) {
		//position is 1 plus the actual found quality because seq.trimFront position is not inclusive
		uint32_t position = seq.qual_.size() - (iter - seq.qual_.rbegin());
		seq.trimFront(position);
	}
}



void readVecTrimmer::trimAtLastBase(seqInfo &seq, const char base){
	auto pos = seq.seq_.rfind(base);
	if(std::string::npos != pos){
		seq.trimBack(pos);
	}
}

void readVecTrimmer::trimAtFirstSeq(seqInfo &seq, const std::string & str){
	auto pos = seq.seq_.find(str);
	if(std::string::npos != pos){
		seq.trimBack(pos);
	}
}

void readVecTrimmer::trimAtLastSeq(seqInfo &seq, const std::string & str){
	auto pos = seq.seq_.rfind(str);
	if(std::string::npos != pos){
		seq.trimBack(pos);
	}
}

void readVecTrimmer::trimOffEndBases(seqInfo &seq, size_t endBases){
	if(len(seq) < endBases){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " : error, length of seq is smaller than the number of bases to trim" << "\n";
		ss << "len of seq: " << len(seq) << ", bases to trim: " << endBases << "\n";
		throw std::runtime_error{ss.str()};
	}
	seq.trimBack(seq.seq_.size() - endBases);
}

void readVecTrimmer::trimOffForwardBases(seqInfo &seq, size_t forwardBases){
	if(len(seq) < forwardBases){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " : error, length of seq is smaller than the number of bases to trim" << "\n";
		ss << "len of seq: " << len(seq) << ", bases to trim: " << forwardBases << "\n";
		throw std::runtime_error{ss.str()};
	}
	seq.trimFront(forwardBases);
}

void readVecTrimmer::trimEnds(seqInfo &seq, size_t forwardBases, size_t endBases){
	trimOffForwardBases(seq, forwardBases);
	trimOffEndBases(seq, endBases);
}

void readVecTrimmer::trimAtSequence(seqInfo &seq, const seqInfo &reversePrimer,
		aligner &alignObj, const comparison & allowableErrors, const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	//first is position of the first base, second is the position of the last base of the match
	std::pair<int, int> rPos;
	if (tSeqPars.within_ >= len(seq)) {
		if (tSeqPars.local_) {
			rPos = alignObj.findReversePrimer(seq, reversePrimer);
			alignObj.rearrangeObjs(seq, reversePrimer, tSeqPars.local_);
			alignObj.profilePrimerAlignment(seq, reversePrimer);
		} else {
			alignObj.alignCacheGlobal(seq, reversePrimer);
			alignObj.rearrangeObjs(seq, reversePrimer, tSeqPars.local_);
			alignObj.profilePrimerAlignment(seq, reversePrimer);
			if('-' == alignObj.alignObjectA_.seqBase_.seq_.front() ||
					'-' != alignObj.alignObjectB_.seqBase_.seq_.front()){
				rPos.first = 0;
			}else{
				//not if if reversePrimer is all - to begin with
				rPos.first = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-'));
			}
			if ('-' == alignObj.alignObjectA_.seqBase_.seq_.back()
					|| '-' != alignObj.alignObjectB_.seqBase_.seq_.back()) {
				rPos.second = len(seq) - 1;
			}else{
				rPos.second = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-'));
			}
		}
	}else{
		auto fromPos = len(seq) - tSeqPars.within_;
		auto backSeq = seq.getSubRead(fromPos);
		if(tSeqPars.local_){
			rPos = alignObj.findReversePrimer(backSeq, reversePrimer);
			alignObj.rearrangeObjs(backSeq, reversePrimer, tSeqPars.local_);
			alignObj.profilePrimerAlignment(backSeq, reversePrimer);
		}else{
			alignObj.alignCacheGlobal(backSeq, reversePrimer);
			alignObj.rearrangeObjs(backSeq, reversePrimer, tSeqPars.local_);
			alignObj.profilePrimerAlignment(backSeq, reversePrimer);
			if('-' == alignObj.alignObjectA_.seqBase_.seq_.front() ||
					'-' != alignObj.alignObjectB_.seqBase_.seq_.front()){
				rPos.first = 0;
			}else{
				//not if if reversePrimer is all - to begin with
				rPos.first = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-'));
			}
			if ('-' == alignObj.alignObjectA_.seqBase_.seq_.back()
					|| '-' != alignObj.alignObjectB_.seqBase_.seq_.back()) {
				rPos.second = len(backSeq) - 1;
			}else{
				rPos.second = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-'));
			}
		}
		rPos.first += fromPos;
		rPos.second += fromPos;
		/*debugging
		alignObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
		alignObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		std::cout << rPos.first << " : " << rPos.second << std::endl;*/
	}



	bool passInspection = allowableErrors.passErrorProfile(alignObj.comp_);
	if (alignObj.comp_.distances_.query_.coverage_
			< allowableErrors.distances_.query_.coverage_) {
		passInspection = false;
	}
	/*debugging
	if(seq.name_ == "contig1[length=273;region=ama1_subRegion;sample=QG0183-C]_t210"){
		std::cout << seq.toJsonJustInfo() << std::endl;
		std::cout << alignObj.comp_.toJson() << std::endl;
		std::cout << "passInspection: " << njh::colorBool(passInspection) << std::endl;

	}*/

	seq.on_ = passInspection;
	if(seq.on_ || tSeqPars.alwaysTrim){
		if (tSeqPars.includeSequence_) {
			if (tSeqPars.sequenceToLowerCase_) {
				if (tSeqPars.removePreviousSameBases_) {
					while (seq.seq_[rPos.first] == seq.seq_[rPos.first - 1]) {
						--rPos.first;
					}
				}
				changeSubStrToLowerToEnd(seq.seq_, rPos.first);
			}
			seq.trimBack(rPos.second + 1);
			//seq.setClip(0, rPos.second);
		} else {
			if (tSeqPars.removePreviousSameBases_) {
				bool foundPrevious = false;
				while (seq.seq_[rPos.first] == seq.seq_[rPos.first - 1]) {
					foundPrevious = true;
					--rPos.first;
				}

				//set clip will keep the second position you give it so you have to minus one more if you found any
				if (foundPrevious) {
					--rPos.first;
				}
			}
			seq.setClip(0, rPos.first - 1);
		}
	}
}

void readVecTrimmer::trimBeforeSequence(seqInfo &seq, const seqInfo &forwardSeq,
		aligner &alignObj, const comparison & allowableErrors, const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	std::pair<int, int> fPos;
	if (tSeqPars.within_ >= len(seq)) {
		if (tSeqPars.local_) {
			fPos = alignObj.findReversePrimer(seq, forwardSeq);
			alignObj.rearrangeObjs(seq, forwardSeq, true);
			alignObj.profilePrimerAlignment(seq, forwardSeq);
		} else {
			alignObj.alignCacheGlobal(seq, forwardSeq);
			alignObj.profilePrimerAlignment(seq, forwardSeq);
			if('-' == alignObj.alignObjectA_.seqBase_.seq_.front() ||
					'-' != alignObj.alignObjectB_.seqBase_.seq_.front()){
				fPos.first = 0;
			}else{
				//not if if forwardSeq is all - to begin with
				fPos.first = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-'));
			}
			if ('-' == alignObj.alignObjectA_.seqBase_.seq_.back()
					|| '-' != alignObj.alignObjectB_.seqBase_.seq_.back()) {
				fPos.second = len(seq) - 1;
			}else{
				fPos.second = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-'));
			}
		}
		//debugging
		/*
		alignObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
		alignObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		std::cout << fPos.first << " : " << fPos.second << std::endl;
	  */
	}else{
		auto frontSeq = seq.getSubRead(0, tSeqPars.within_);
		if(tSeqPars.local_){
			fPos = alignObj.findReversePrimer(frontSeq, forwardSeq);
			alignObj.rearrangeObjs(frontSeq, forwardSeq, tSeqPars.local_);
			alignObj.profilePrimerAlignment(frontSeq, forwardSeq);
		}else{
			alignObj.alignCacheGlobal(frontSeq, forwardSeq);
			alignObj.profilePrimerAlignment(frontSeq, forwardSeq);
			if('-' == alignObj.alignObjectA_.seqBase_.seq_.front() ||
					'-' != alignObj.alignObjectB_.seqBase_.seq_.front()){
				fPos.first = 0;
			}else{
				//not if if forwardSeq is all - to begin with
				fPos.first = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-'));
			}
			if ('-' == alignObj.alignObjectA_.seqBase_.seq_.back()
					|| '-' != alignObj.alignObjectB_.seqBase_.seq_.back()) {
				fPos.second = len(frontSeq) - 1;
			}else{
				fPos.second = alignObj.getSeqPosForAlnAPos(alignObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-'));
			}
		}
		//debugging
		/*
		alignObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
		alignObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
		std::cout << fPos.first << " : " << fPos.second << std::endl;
		*/
	}
	/*
	std::pair<int, int> fPos = alignObj.findReversePrimer(seq, forwardSeq);
	alignObj.rearrangeObjs(seq, forwardSeq, true);
	alignObj.profilePrimerAlignment(seq, forwardSeq);
*/


	bool passInspection = allowableErrors.passErrorProfile(alignObj.comp_);
	if (alignObj.comp_.distances_.query_.coverage_
			< allowableErrors.distances_.query_.coverage_) {
		passInspection = false;
	}
	//debugging
	/*
	if(!passInspection){
		std::cout << seq.toJsonJustInfo() << std::endl;
		std::cout << alignObj.comp_.toJson() << std::endl;
		std::cout << "passInspection: " << njh::colorBool(passInspection) << std::endl;
	}*/
	seq.on_ = passInspection;

	if(seq.on_ || tSeqPars.alwaysTrim){
		if (tSeqPars.includeSequence_) {
			if (tSeqPars.sequenceToLowerCase_) {
				if (tSeqPars.removePreviousSameBases_) {
					while (seq.seq_[fPos.second] == seq.seq_[fPos.second + 1]) {
						++fPos.second;
					}
				}
				changeSubStrToLowerFromBegining(seq.seq_, fPos.second);
			}
			seq.trimFront(fPos.first);
		} else {
			if (tSeqPars.removePreviousSameBases_) {
				while (seq.seq_[fPos.second] == seq.seq_[fPos.second + 1]) {
					++fPos.second;
				}
			}
			//debugging
			/*
			std::cout << "Trimming to " << fPos.second + 1 << std::endl;
			*/
			seq.trimFront(fPos.second + 1);
		}
	}
}


void readVecTrimmer::trimBetweenSequences(seqInfo &seq, const seqInfo &forwardSeq,
		const seqInfo &backSeq, aligner &alignObj, const comparison & allowableErrors,
		const FullTrimReadsPars::trimSeqPars & tSeqPars) {
	trimAtSequence(seq, backSeq, alignObj, allowableErrors, tSeqPars);
	if(seq.on_){
		trimBeforeSequence(seq, forwardSeq, alignObj, allowableErrors, tSeqPars);
	}
}




}  // namespace njhseq
