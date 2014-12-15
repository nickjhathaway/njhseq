//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "seqInfo.hpp"
#include "bibseq/helpers/seqUtil.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {

seqInfo seqInfo::getSubRead(uint32_t pos, uint32_t size) const{
  return seqInfo(name_, seq_.substr(pos, size), getSubVector(qual_, pos, size),
                 cnt_, frac_);
}
seqInfo seqInfo::getSubRead(uint32_t pos) const{
  return seqInfo(name_, seq_.substr(pos, seq_.size() - pos),
                 getSubVector(qual_, pos, len(seq_) - pos), cnt_, frac_);
}

// description
void seqInfo::printDescription(std::ostream& out, bool deep) const {
  out << "seqInfo{" << std::endl << "name_:" << name_ << std::endl
      << "seq_:" << seq_ << std::endl << "qual_:" << qual_ << std::endl
      << "cnt_:" << cnt_ << std::endl << "frac_:" << frac_ << std::endl << "}"
      << std::endl;
}

Json::Value seqInfo::toJson()const{
	Json::Value ret;
	ret["seq"] = bib::json::toJson(seq_);
	ret["qual"] = bib::json::toJson(qual_);
	ret["cnt"] = bib::json::toJson(cnt_);
	ret["frac"] = bib::json::toJson(frac_);
	ret["name"] = bib::json::toJson(name_);
	return ret;
}

const std::unordered_map<char, uint32_t> seqInfo::ansiBaseColor = std::unordered_map<char, uint32_t>{
		{'A',203},
		{'a',203},
		{'C',83},
		{'c',83},
		{'G',227},
		{'g',227},
		{'T',69},
		{'t',69},
		{'n',145},
		{'N',145},
		{'-',102}
};

void seqInfo::outPutSeqAnsi(std::ostream& fastaFile) const{
	fastaFile << bib::bashCT::addBGColor(145) << ">" << name_ << bib::bashCT::reset << std::endl;
	for(const auto & c : seq_){
		fastaFile << bib::bashCT::addBGColor(ansiBaseColor.at(c)) << c << bib::bashCT::reset;
	}
	fastaFile << bib::bashCT::reset << std::endl;
}
void seqInfo::prepend(const std::string& seq, uint32_t defaultQuality){
	prepend(seq, std::vector<uint32_t> (1,defaultQuality));
}
void seqInfo::append(const std::string& seq, uint32_t defaultQuality){
	append(seq, std::vector<uint32_t> (1,defaultQuality));
}
void seqInfo::append(const std::string& seq,
                     const std::vector<uint32_t>& qual) {
  if (qual.size() == 1) {
    seq_.append(seq);
    addOtherVec(qual_, std::vector<uint32_t>(seq.size(), qual.front()));
  } else if (qual.size() == seq.size()) {
    seq_.append(seq);
    addOtherVec(qual_, qual);
  } else {
    std::cout << "Need to supply either single a quality or same amount of "
                 "quality score for seq length" << std::endl;
    std::cout << "trying to add " << qual.size()
              << " qualities and a seq of length " << seq.length() << std::endl;
    exit(1);
  }
}

void seqInfo::reverseComplementRead(bool mark) {
  std::vector<std::vector<uint32_t>> quals;
  std::vector<uint32_t> currentQuals;
  currentQuals.push_back(qual_[0]);
  for (uint32_t i = 1; i < seq_.length(); ++i) {

    if (seq_[i] == seq_[i - 1]) {
      currentQuals.push_back(qual_[i]);
    } else {
      quals.push_back(currentQuals);
      currentQuals.clear();
      currentQuals.push_back(qual_[i]);
    }
  }
  quals.push_back(currentQuals);
  qual_.clear();
  seq_ = seqUtil::reverseComplement(seq_, "DNA");
  for (auto iter = quals.rbegin(); iter != quals.rend(); ++iter) {
    addOtherVec(qual_, *iter);
  }
  if(mark){
  	name_.append("_Comp");
  }
}
void seqInfo::prepend(const std::string& seq,
                      const std::vector<uint32_t>& qual) {
  if (qual.size() == 1) {
    seq_.insert(seq_.begin(), seq.begin(), seq.end());
    prependVec(qual_, std::vector<uint32_t>(seq.size(), qual.front()));
  } else if (qual.size() == seq.size()) {
    seq_.insert(seq_.begin(), seq.begin(), seq.end());
    prependVec(qual_, qual);
  } else {
    std::cout << "Need to supply either single a quality or same amount of "
                 "quality score for seq length" << std::endl;
    std::cout << "trying to add " << qual.size()
              << " qualities and a seq of length " << seq.length() << std::endl;
    exit(1);
  }
}
const std::vector<uint32_t> seqInfo::getLeadQual(int posA, int out) const {
  std::vector<uint32_t> ans;
  int lowerBound = 0;
  if (posA - out > lowerBound) {
    lowerBound = posA - out;
  }
  for (auto i : iter::range(lowerBound, posA)) {
    ans.push_back(qual_[i]);
  }
  return ans;
}
const std::vector<uint32_t> seqInfo::getTrailQual(int posA, int out) const {
  std::vector<uint32_t> ans;
  int higherBound = (int)qual_.size() - 1;
  if (posA + out + 1 < higherBound) {
    higherBound = posA + out + 1;
  }
  for (auto i : iter::range(posA + 1, higherBound)) {
    ans.push_back(qual_[i]);
  }
  return ans;
}
bool seqInfo::checkLeadQual(int pos, uint32_t secondayQual, int out) const {
  int lowerBound = 0;
  if (pos - out > lowerBound) {
    lowerBound = pos - out;
  }
  for (auto i : iter::range(lowerBound, pos)) {
    if (qual_[i] <= secondayQual) {
      return false;
    }
  }
  return true;
}
bool seqInfo::checkTrailQual(int pos, uint32_t secondayQual, int out) const {
  int higherBound = (int)qual_.size() - 1;
  if (pos + out + 1 < higherBound) {
    higherBound = pos + out + 1;
  }
  for (auto i : iter::range(pos + 1, higherBound)) {
    if (qual_[i] <= secondayQual) {
      return false;
    }
  }
  return true;
}
bool seqInfo::checkPrimaryQual(int pos, uint32_t primaryQual) const {
  if (qual_[pos] <= primaryQual) {
    return false;
  }
  return true;
}
bool seqInfo::checkQual(int pos, uint32_t primaryQual, uint32_t secondayQual,
                        int out) const {
  if (!checkPrimaryQual(pos, primaryQual)) {
    return false;
  }
  if (!checkLeadQual(pos, secondayQual, out) ||
      !checkTrailQual(pos, secondayQual, out)) {
    return false;
  }
  return true;
}
int seqInfo::findLowestNeighborhoodQual(int posA, int out) const {
  int lowerBound = 0;
  int higherBound = (int)qual_.size() - 1;
  if (posA - out > lowerBound) {
    lowerBound = posA - out;
  }
  if (posA + out + 1 < higherBound) {
    higherBound = posA + out + 1;
  }
  uint32_t lowestQual = UINT32_MAX;
  for (auto i : iter::range(lowerBound, higherBound)) {
    if (posA != i) {
      if (qual_[i] < lowestQual) {
        lowestQual = qual_[i];
      }
    }
  }
  return lowestQual;
}
void seqInfo::removeBase(size_t pos) {
  if (pos >= seq_.size()) {
    std::cout << "pos: " << pos << " out of bounds of " << seq_.size()
              << std::endl;
  } else {
    seq_.erase(seq_.begin() + pos);
    qual_.erase(qual_.begin() + pos);
  }
}

void seqInfo::removeLowQualityBases(uint32_t qualCutOff) {
  std::vector<size_t> positions;
  size_t pos = 0;
  for (auto qIter : qual_) {
    if (qIter < qualCutOff) {
      positions.push_back(pos);
    }
    ++pos;
  }
  std::reverse(positions.begin(), positions.end());
  for (const auto& pIter : positions) {
    removeBase(pIter);
  }
}
void seqInfo::removeGaps() {
  for (auto i : iter::range((int)seq_.size(), -1, -1)) {
    if (seq_[i] == '-') {
      removeBase(i);
    }
  }
  return;
}

// quality strings for printing
std::string seqInfo::getQualString() const { return vectorToString(qual_); }
std::string seqInfo::getFastqString(const std::vector<uint32_t>& quals,
                                    uint32_t offset) {
  std::string convertedQuals = "";
  for (const auto& q : quals) {
    if (q <= 93) {
      convertedQuals.push_back((char)(q + offset));
    } else {
      convertedQuals.push_back((char)(93 + offset));
    }
  }
  return convertedQuals;
}
std::string seqInfo::getFastqQualString(int offset) const {
  return getFastqString(qual_, offset);
}

//
void seqInfo::markAsChimeric() {
  if (name_.find("CHI") == std::string::npos) {
    name_ = "CHI_" + name_;
  }
}

// outputs
void seqInfo::outPutFastq(std::ostream& fastqFile) const {
  fastqFile << "@" << name_ << std::endl;
  fastqFile << seq_ << std::endl;
  fastqFile << "+" << std::endl;
  fastqFile << getFastqQualString(33) << std::endl;
}

void seqInfo::outPutSeq(std::ostream& fastaFile) const {
  fastaFile << ">" << name_ << std::endl;
  fastaFile << seq_ << std::endl;
}

void seqInfo::outPutQual(std::ostream& qualFile) const {
  qualFile << ">" << name_ << std::endl;
  qualFile << getQualString() << std::endl;
}

bool seqInfo::degenCompare(const seqInfo & otherInfo, const substituteMatrix & compareScores)const{
	if(seq_.length() != otherInfo.seq_.length()){
		return false;
	}
	auto comp = [&](const char & a, const char & b){
		return compareScores.mat_[a][b] > 0;
	};
	return std::equal(seq_.begin(), seq_.end(), otherInfo.seq_.begin(), comp);
}

}  // namespace bib
