/*
 * PairedRead.cpp
 *
 *  Created on: Jul 26, 2015
 *      Author: nick
 */
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
#include "PairedRead.hpp"

namespace njhseq {
PairedRead::PairedRead() :
		readObject(), mateSeqBase_(seqInfo()) {

}

PairedRead::PairedRead(const seqInfo & seqFirst, const seqInfo & seqSecond,
		bool processed) :
		readObject(seqFirst, processed), mateSeqBase_(seqSecond) {

}

void PairedRead::toggleReverseComplement(){
	mateSeqBase_.reverseComplementRead(true, true);
	mateRComplemented_ = !mateRComplemented_;
}

double PairedRead::getQualCheck(uint32_t qualCutOff)const{
	uint32_t count = njh::count_if(seqBase_.qual_, [&qualCutOff](uint32_t qual){ return qual >=qualCutOff;});
	count += njh::count_if(mateSeqBase_.qual_, [&qualCutOff](uint32_t qual){ return qual >=qualCutOff;});
	return static_cast<double>(count)/(mateSeqBase_.qual_.size() + seqBase_.qual_.size());
}

void PairedRead::setBaseCountOnQualCheck(uint32_t qualCheck){
	fractionAboveQualCheck_ =  getQualCheck(qualCheck);
}


void PairedRead::setLetterCount() {
	counter_.reset();
	counter_.increaseCountByString(seqBase_.seq_, seqBase_.cnt_);
	if(mateRComplemented_){
		counter_.increaseCountByString(mateSeqBase_.seq_, seqBase_.cnt_);
	}else{
		counter_.increaseCountByString(seqUtil::reverseComplement(mateSeqBase_.seq_, "DNA"), seqBase_.cnt_);
	}
	counter_.resetAlphabet(true);
	counter_.setFractions();
}

void PairedRead::setLetterCount(const std::vector<char> & alph){
	counter_ = charCounter(alph);
	setLetterCount();
}


double PairedRead::getAverageErrorRate() const{
  double sum = 0;
  for (const auto& q : seqBase_.qual_) {
    sum += pow(10.0, -(q / 10.0));
  }
  for (const auto& q : mateSeqBase_.qual_) {
    sum += pow(10.0, -(q / 10.0));
  }
  return sum / len(*this);
}

void PairedRead::createCondensedSeq(){
  condensedSeq = "";
  condensedSeqQual.clear();
  condensedSeqCount.clear();
  condensedSeqQualPos.clear();
  {
    int currentCount = 1;
    std::vector<uint32_t> currentQuals;
    currentQuals.push_back(seqBase_.qual_[0]);
    std::pair<uint32_t, uint32_t> currentQualPos {0,1};
    uint32_t i = 1;
    for (; i < seqBase_.seq_.length(); i++) {
      if (seqBase_.seq_[i] == seqBase_.seq_[i - 1]) {
        currentQuals.push_back(seqBase_.qual_[i]);
        ++currentCount;
      } else {
        condensedSeq.push_back(seqBase_.seq_[i - 1]);
        condensedSeqQual.push_back(vectorMean(currentQuals));
        currentQualPos.second = currentQuals.size();
        condensedSeqQualPos.emplace_back(currentQualPos);
        currentQualPos.first = i;
        condensedSeqCount.push_back(currentCount);
        currentCount = 1;
        currentQuals.clear();
        currentQuals.push_back(seqBase_.qual_[i]);
      }
    }
    condensedSeq.push_back(seqBase_.seq_[i - 1]);
    condensedSeqQual.push_back(vectorMean(currentQuals));
    currentQualPos.second = currentQuals.size();
    condensedSeqQualPos.emplace_back(currentQualPos);
    condensedSeqCount.push_back(currentCount);
  }
  {
    int currentCount = 1;
    std::vector<uint32_t> currentQuals;
    currentQuals.push_back(mateSeqBase_.qual_[0]);
    std::pair<uint32_t, uint32_t> currentQualPos {0,1};
    uint32_t i = 1;
    for (; i < mateSeqBase_.seq_.length(); i++) {
      if (mateSeqBase_.seq_[i] == mateSeqBase_.seq_[i - 1]) {
        currentQuals.push_back(mateSeqBase_.qual_[i]);
        ++currentCount;
      } else {
        condensedSeq.push_back(mateSeqBase_.seq_[i - 1]);
        condensedSeqQual.push_back(vectorMean(currentQuals));
        currentQualPos.second = currentQuals.size();
        condensedSeqQualPos.emplace_back(currentQualPos);
        currentQualPos.first = i;
        condensedSeqCount.push_back(currentCount);
        currentCount = 1;
        currentQuals.clear();
        currentQuals.push_back(mateSeqBase_.qual_[i]);
      }
    }
    condensedSeq.push_back(mateSeqBase_.seq_[i - 1]);
    condensedSeqQual.push_back(vectorMean(currentQuals));
    currentQualPos.second = currentQuals.size();
    condensedSeqQualPos.emplace_back(currentQualPos);
    condensedSeqCount.push_back(currentCount);
  }
}

void PairedRead::outFastq(std::ostream & firstOut, std::ostream & secondOut) const {
	seqBase_.outPutFastq(firstOut);
	mateSeqBase_.outPutFastq(secondOut);
}

} /* namespace njhseq */
