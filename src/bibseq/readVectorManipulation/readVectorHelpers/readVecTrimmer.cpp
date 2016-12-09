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
#include "readVecTrimmer.hpp"




namespace bibseq {



void readVecTrimmer::trimAtFirstBase(seqInfo &seq, const char base){
	auto pos = seq.seq_.find(base);
	if(std::string::npos != pos){
		seq.trimBack(pos);
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

void readVecTrimmer::trimAtSequence(seqInfo &seq, seqInfo &reversePrimer,
		aligner &alignObj, comparison allowableErrors, trimSeqPars tSeqPars) {
	std::pair<int, int> rPos = alignObj.findReversePrimer(seq, reversePrimer);
	alignObj.rearrangeObjs(seq, reversePrimer, true);

	alignObj.profilePrimerAlignment(seq, reversePrimer);
	bool passInspection = allowableErrors.passErrorProfile(alignObj.comp_);
	if (alignObj.comp_.distances_.query_.coverage_
			< allowableErrors.distances_.query_.coverage_) {
		passInspection = false;
	}
	seq.on_ = passInspection;
	if(seq.on_){
		if (tSeqPars.includeSequence_) {
			seq.setClip(0, rPos.second);
			if (tSeqPars.sequenceToLowerCase_) {
				if (tSeqPars.removePreviousSameBases_) {
					while (seq.seq_[rPos.first] == seq.seq_[rPos.first - 1]) {
						--rPos.first;
					}
				}
				changeSubStrToLowerToEnd(seq.seq_, rPos.first);
			}
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

void readVecTrimmer::trimBeforeSequence(seqInfo &seq, seqInfo &forwardSeq,
		aligner &alignObj, comparison allowableErrors, trimSeqPars tSeqPars) {

	std::pair<int, int> fPos = alignObj.findReversePrimer(seq, forwardSeq);
	alignObj.rearrangeObjs(seq, forwardSeq, true);
	alignObj.profilePrimerAlignment(seq, forwardSeq);



	bool passInspection = allowableErrors.passErrorProfile(alignObj.comp_);
	if (alignObj.comp_.distances_.query_.coverage_
			< allowableErrors.distances_.query_.coverage_) {
		passInspection = false;
	}
	seq.on_ = passInspection;
	if(seq.on_){
		if (tSeqPars.includeSequence_) {
			seq.trimFront(fPos.first);
			if (tSeqPars.sequenceToLowerCase_) {
				if (tSeqPars.removePreviousSameBases_) {
					while (seq.seq_[fPos.second] == seq.seq_[fPos.second + 1]) {
						++fPos.second;
					}
				}
				changeSubStrToLowerFromBegining(seq.seq_, fPos.second);
			}
		} else {
			if (tSeqPars.removePreviousSameBases_) {
				while (seq.seq_[fPos.second] == seq.seq_[fPos.second + 1]) {
					++fPos.second;
				}
			}
			seq.trimFront(fPos.second + 1);
		}
	}

}


void readVecTrimmer::trimBetweenSequences(seqInfo &seq, seqInfo &forwardSeq,
		seqInfo &backSeq, aligner &alignObj, comparison allowableErrors,
		trimSeqPars tSeqPars) {
	trimAtSequence(seq, backSeq, alignObj, allowableErrors, tSeqPars);
	if(seq.on_){
		trimBeforeSequence(seq, forwardSeq, alignObj, allowableErrors, tSeqPars);
	}
}




}  // namespace bibseq
