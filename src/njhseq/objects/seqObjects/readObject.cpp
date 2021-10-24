#include "readObject.hpp"

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
namespace njhseq {


readObject::readObject(const seqInfo& seqBase, bool processed)
    : baseReadObject(seqBase) {
  processRead(processed);
  initializeOthers();
  //std::cout << "readObject constructor: " << std::endl;
  //std::cout << seqBase_.name_ << std::endl;
  //std::cout << seqBase_.cnt_ << std::endl;
  //std::cout << seqBase_.frac_ << std::endl;
}

readObject::readObject() :
		baseReadObject() {
	initializeOthers();
}


void readObject::resetMetaInName() {
	seqBase_.resetMetaInName(meta_);
}

bool readObject::nameHasMetaData() const {
	return seqBase_.nameHasMetaData();
}

void readObject::processNameForMeta() {
	meta_.meta_.clear();
	meta_.processNameForMeta(seqBase_.name_, false);
}

bool readObject::containsMeta(const std::string & key) const {
	return meta_.containsMeta(key);
}

std::string readObject::getMeta(const std::string & key) const {
	return meta_.getMeta(key);
}


Json::Value readObject::toJson()const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["super"] = baseReadObject::toJson();
	ret["sampName"] = njh::json::toJson(sampName);
	ret["expectsString"] = njh::json::toJson(expectsString);
	ret["meta_"] = njh::json::toJson(meta_);
	ret["averageErrorRate"] = njh::json::toJson(averageErrorRate);
	ret["fractionAboveQualCheck_"] = njh::json::toJson(fractionAboveQualCheck_);
	ret["remove"] = njh::json::toJson(remove);
	ret["counter_"] = njh::json::toJson(counter_);
	ret["condensedSeq"] = njh::json::toJson(condensedSeq);
	ret["condensedSeqQual"] = njh::json::toJson(condensedSeqQual);
	ret["condensedSeqQualPos"] = njh::json::toJson(condensedSeqQualPos);
	ret["condensedSeqCount"] = njh::json::toJson(condensedSeqCount);
	return ret;
}



std::string readObject::getOwnSampName() const {
  return seqBase_.getOwnSampName();
}

/// chagening the name of the cluster
void readObject::updateName() {
	seqBase_.updateName();
}

void readObject::setName(const std::string& newName) {
  seqBase_.name_ = newName + "_t" + estd::to_string(seqBase_.cnt_);
  sampName = getOwnSampName();
}


std::string readObject::getReadId() const {
  return seqBase_.getReadId();
}
std::string readObject::getStubName(bool removeChiFlag) const {
  return seqBase_.getStubName(removeChiFlag);
}




void readObject::appendName(const std::string& add) {
  size_t totalPos = seqBase_.name_.rfind("_t");
  if (totalPos != std::string::npos) {
    seqBase_.name_ = seqBase_.name_.substr(0, totalPos);
  }
  seqBase_.name_ += add + "_t" + estd::to_string(seqBase_.cnt_);
}

void readObject::processRead(bool processed) {
	seqBase_.processRead(processed);
}


void readObject::initializeOthers() {

  condensedSeq = "";
  condensedSeqCount.clear();
  averageErrorRate = getAverageErrorRate();
  remove = false;

  //basesAboveQualCheck_ = 0;
  fractionAboveQualCheck_ = 0;
  sampName = getOwnSampName();
}

void readObject::addQual(const std::string & qualString) {
	seqBase_.addQual(qualString);
  averageErrorRate = getAverageErrorRate();
}

void readObject::addQual(const std::string & qualString, uint32_t offSet) {
	seqBase_.addQual(qualString, offSet);
  averageErrorRate = getAverageErrorRate();
}

void readObject::addQual(const std::vector<uint8_t> & quals){
	seqBase_.addQual(quals);
  averageErrorRate = getAverageErrorRate();
}

void readObject::clearObject() {
  seqBase_.name_ = "";
  seqBase_.seq_ = "";
  seqBase_.qual_.clear();
  seqBase_.cnt_ = 0;
  condensedSeq = "";
  condensedSeqCount.clear();
}

void readObject::createCondensedSeq() {
  condensedSeq = "";
  condensedSeqQual.clear();
  condensedSeqCount.clear();
  condensedSeqQualPos.clear();
  int currentCount = 1;
  std::vector<uint32_t> currentQuals;
  currentQuals.push_back(seqBase_.qual_[0]);
  std::pair<uint32_t, uint32_t> currentQualPos {0,1};
//  bool print = false;
//  if (seqBase_.name_.find("WGHHT:02045:02871_t46.000000") !=
//      std::string::npos) {
//    print = false;
//  }
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
//      if (print) {
//        std::cout << "Base: " << seqBase_.seq_[i - 1] << std::endl;
//        std::cout << "\t" << vectorToString(currentQuals) << std::endl;
//        std::cout << "\t" << currentCount << std::endl;
//      }
      currentCount = 1;
      currentQuals.clear();
      currentQuals.push_back(seqBase_.qual_[i]);
    }
  }
//  if (print) {
//    std::cout << "Base: " << seqBase_.seq_[i - 1] << std::endl;
//    std::cout << "\t" << vectorToString(currentQuals) << std::endl;
//    std::cout << "\t" << currentCount << std::endl;
//  }
  condensedSeq.push_back(seqBase_.seq_[i - 1]);
  condensedSeqQual.push_back(vectorMean(currentQuals));
  currentQualPos.second = currentQuals.size();
  condensedSeqQualPos.emplace_back(currentQualPos);
  condensedSeqCount.push_back(currentCount);
//  if (print) {
//    std::cout << condensedSeq << std::endl;
//    std::cout << vectorToString(condensedSeqQual) << std::endl;
//    std::cout << vectorToString(condensedSeqCount) << std::endl;
//    std::cout << seqBase_.seq_ << std::endl;
//    std::cout << vectorToString(seqBase_.qual_) << std::endl;
//    std::cout << seqBase_.seq_.size() << std::endl;
//    std::cout << seqBase_.qual_.size() << std::endl;
//    checkSeqQual(std::cout);
//    exit(1);
//  }
}

void readObject::setClip(size_t leftPos, size_t rightPos) {
	seqBase_.setClip(leftPos, rightPos);
	averageErrorRate = getAverageErrorRate();
}
void readObject::setClip(size_t rightPos) {
	setClip(0, rightPos);
}
void readObject::setClip(const std::pair<int, int>& positions) {
	setClip(positions.first, positions.second);
}

void readObject::trimFront(size_t upToPosNotIncluding) {
  setClip(upToPosNotIncluding, seqBase_.seq_.size() - 1);
}
void readObject::trimBack(size_t fromPositionIncluding) {
  setClip(0, fromPositionIncluding - 1);

}

double readObject::getAverageQual() const {
	return seqBase_.getAverageQual();
}

double readObject::getAverageErrorRate() const {
	return seqBase_.getAverageErrorRate();
}

uint32_t readObject::getSumQual() const {
	return seqBase_.getSumQual();
}

void readObject::setFractionByCount(double totalNumberOfReads) {
  seqBase_.setFractionByCount(totalNumberOfReads);
}


std::vector<size_t> readObject::findSubsequenceOccurences(
    const std::string& search) const {
  return findOccurences(seqBase_.seq_, search);
}





void readObject::outPutCondensedSeq(std::ostream& out) const {
  out << ">" << seqBase_.name_ << std::endl;
  out << condensedSeq << std::endl;
}
void readObject::outPutCondensedQual(std::ostream& out) const {
  out << ">" << seqBase_.name_ << std::endl;
  out << vectorToString(condensedSeqQual) << std::endl;
}

// special output
void readObject::checkSeqQual(std::ostream& outFile) const {
  for (uint32_t i = 0; i < seqBase_.seq_.length(); ++i) {
    outFile << seqBase_.seq_[i] << ":" << seqBase_.qual_[i] << " ";
  }
  outFile << std::endl;
}

// counter setters
void readObject::setLetterCount() {
	counter_.reset();
	counter_.increaseCountByString(seqBase_.seq_, seqBase_.cnt_);
	counter_.resetAlphabet(true);
	counter_.setFractions();
}

void readObject::setLetterCount(const std::vector<char> & alph){
	counter_ = charCounter(alph);
	counter_.increaseCountByString(seqBase_.seq_, seqBase_.cnt_);
	counter_.setFractions();
}



// comparisons
bool readObject::operator<(const readObject& otherRead) const {
	/*
	if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
		if (averageErrorRate < otherRead.averageErrorRate) {
			return true;
		} else {
			return false;
		}
	} else {
		return seqBase_.cnt_ > otherRead.seqBase_.cnt_;
	}*/
	if (roundDecPlaces(seqBase_.cnt_, 2) == roundDecPlaces(otherRead.seqBase_.cnt_, 2) ) {
		if (roundDecPlaces(averageErrorRate, 2)  < roundDecPlaces(otherRead.averageErrorRate, 2) ) {
			return true;
		} else {
			return false;
		}
	} else {
		return roundDecPlaces(seqBase_.cnt_, 2)  > roundDecPlaces(otherRead.seqBase_.cnt_, 2) ;
	}
}

bool readObject::operator>(const readObject& otherRead) const {
	/*
	if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
		if (averageErrorRate > otherRead.averageErrorRate) {
			return true;
		} else {
			return false;
		}
	} else {
		return seqBase_.cnt_ < otherRead.seqBase_.cnt_;
	}*/
	if (roundDecPlaces(seqBase_.cnt_, 2) == roundDecPlaces(otherRead.seqBase_.cnt_, 2) ) {
		if (roundDecPlaces(averageErrorRate, 2)  > roundDecPlaces(otherRead.averageErrorRate, 2) ) {
			return true;
		} else {
			return false;
		}
	} else {
		return roundDecPlaces(seqBase_.cnt_, 2)  < roundDecPlaces(otherRead.seqBase_.cnt_, 2) ;
	}
}

/*
bool readObject::operator==(const readObject& otherRead) const {
  if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
    return true;
  } else {
    return false;
  }
}*/

bool readObject::operator>=(const readObject& otherRead) const {
	/*if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
		if (averageErrorRate >= otherRead.averageErrorRate) {
			return true;
		} else {
			return false;
		}
	} else {
		return seqBase_.cnt_ <= otherRead.seqBase_.cnt_;
	}*/
	if (roundDecPlaces(seqBase_.cnt_, 2) == roundDecPlaces(otherRead.seqBase_.cnt_, 2) ) {
		if (roundDecPlaces(averageErrorRate, 2)  >= roundDecPlaces(otherRead.averageErrorRate, 2) ) {
			return true;
		} else {
			return false;
		}
	} else {
		return roundDecPlaces(seqBase_.cnt_, 2)  <= roundDecPlaces(otherRead.seqBase_.cnt_, 2) ;
	}
}

bool readObject::operator<=(const readObject& otherRead) const {
	/*if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
		if (averageErrorRate <= otherRead.averageErrorRate) {
			return true;
		} else {
			return false;
		}
	} else {
		return seqBase_.cnt_ >= otherRead.seqBase_.cnt_;
	}*/
	if (roundDecPlaces(seqBase_.cnt_, 2) == roundDecPlaces(otherRead.seqBase_.cnt_, 2) ) {
		if (roundDecPlaces(averageErrorRate, 2)  <= roundDecPlaces(otherRead.averageErrorRate, 2) ) {
			return true;
		} else {
			return false;
		}
	} else {
		return roundDecPlaces(seqBase_.cnt_, 2)  >= roundDecPlaces(otherRead.seqBase_.cnt_, 2) ;
	}
}

double readObject::getQualCheck(uint32_t qualCutOff) const{
  auto basesAboveQualCheck = njh::count_if(seqBase_.qual_, [&qualCutOff](const uint32_t & q){ return q>=qualCutOff;});
  return static_cast<double>(basesAboveQualCheck) / seqBase_.qual_.size();
}

void readObject::setBaseCountOnQualCheck(uint32_t qualCheck) {
  fractionAboveQualCheck_ = getQualCheck(qualCheck);
}

void readObject::replace(const std::string& toBeReplaced,
                         const std::string& replaceWith, bool allOccurences) {
  std::vector<size_t> occurences = findSubsequenceOccurences(toBeReplaced);
  std::reverse(occurences.begin(), occurences.end());
  for (const auto pos : occurences) {
    std::vector<uint32_t> currentQuals;
    for (const auto subPos : iter::range(pos, pos + toBeReplaced.size())) {
      currentQuals.push_back(seqBase_.qual_[subPos]);
    }
    for (size_t i = 0; i < toBeReplaced.size(); ++i) {
      seqBase_.seq_.erase(seqBase_.seq_.begin() + pos);
      seqBase_.qual_.erase(seqBase_.qual_.begin() + pos);
    }
    double meanQual = vectorMean(currentQuals);
    seqBase_.qual_.insert(seqBase_.qual_.begin() + pos, replaceWith.size(),
                          meanQual);
    seqBase_.seq_.insert(pos, replaceWith.c_str());
  }
}




void readObject::translate(bool complement, bool reverse, size_t start ){
	seqBase_.translate(complement, reverse, start);
}
readObject readObject::translateRet(bool complement, bool reverse, size_t start ) const{
	return readObject(seqBase_.translateRet(complement, reverse, start), false);
}


double readObject::getGCContent() {
  setLetterCount();
  counter_.calcGcContent();
  return counter_.gcContent_;
}

void readObject::adjustHomopolyerRunQualities() {
  createCondensedSeq();
  seqBase_.qual_.clear();
  for (const auto i : iter::range<uint64_t>(0, condensedSeq.length())) {
    addOtherVec(seqBase_.qual_, std::vector<uint8_t>(condensedSeqCount[i],
                                                      condensedSeqQual[i]));
  }
}
void readObject::updateQualCounts(std::map<uint32_t, uint32_t>& qualCounts)
    const {
  for (const auto& q : seqBase_.qual_) {
    ++qualCounts[q];
  }
  return;
}
void readObject::updateQualCounts(
    std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize, const std::array<double, 100>& qualErrorLookUp) const {

  for (const auto pos : iter::range(seqBase_.qual_.size())) {
    updateQaulCountsAtPos(pos, counts, qualWindowSize, qualErrorLookUp);
  }
}
void readObject::updateQualWindowCounts(
    uint32_t pos, std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize) const {
  std::vector<uint8_t> currentQuals;
  currentQuals.push_back(seqBase_.qual_[pos]);
  addOtherVec(currentQuals, seqBase_.getLeadQual(pos, qualWindowSize));
  addOtherVec(currentQuals, seqBase_.getTrailQual(pos, qualWindowSize));
  counts["mean"][roundDecPlaces(vectorMean(currentQuals), 2)] += seqBase_.cnt_;
  counts["median"][roundDecPlaces(vectorMedianRef(currentQuals), 2)] +=
      seqBase_.cnt_;
  counts["min"][roundDecPlaces(vectorMinimum(currentQuals), 2)] +=
      seqBase_.cnt_;
  return;
}
void readObject::updateQualWindowCounts(
    std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize) const {
  for (const auto pos : iter::range(seqBase_.qual_.size())) {
    updateQualWindowCounts(pos, counts, qualWindowSize);
  }
}

void readObject::updateBaseQualCounts(std::map<double, uint32_t>& baseCounts,
                                      uint32_t pos) const {
  baseCounts[roundDecPlaces(seqBase_.qual_[pos], 2)] += seqBase_.cnt_;
}

void readObject::updateBaseQualCounts(std::map<double, uint32_t>& baseCounts)
    const {
  for (const auto pos : iter::range(seqBase_.qual_.size())) {
    updateBaseQualCounts(baseCounts, pos);
  }
}
void readObject::updateQaulCountsAtPos(
    uint32_t pos, std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize, const std::array<double, 100>& qualErrorLookUp) const {

  std::vector<uint8_t> currentQuals;
  currentQuals.push_back(seqBase_.qual_[pos]);
  addOtherVec(currentQuals, seqBase_.getLeadQual(pos, qualWindowSize));
  addOtherVec(currentQuals, seqBase_.getTrailQual(pos, qualWindowSize));
  auto stats = getStatsOnVec(currentQuals);
  ++counts["base"][seqBase_.qual_[pos]];
  for (const auto& stat : stats) {
    if (stat.first == "std") {
      continue;
    } else if ("mean" == stat.first) {
      ++counts[stat.first][roundDecPlaces(stat.second, 2)];
    } else {
      ++counts[stat.first][stat.second];
    }
  }
  double currentProbSum = 0.00;
  for (const auto& q : currentQuals) {
    currentProbSum += qualErrorLookUp[q];
  }
  ++counts["probSum"][roundDecPlaces(currentProbSum, 4)];
  ++counts["probMean"]
          [roundDecPlaces(currentProbSum / (2 * qualWindowSize + 1), 4)];
}





readObject::~readObject(){}

}  // namespace njh
