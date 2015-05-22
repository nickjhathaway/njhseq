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

#include "readObject.hpp"
namespace bibseq {

std::string readObject::getOtherReadSampName(const readObject& cr) const {
  VecStr toks = tokenizeString(cr.seqBase_.name_, ".");
  return replaceString(toks[0], "CHI_", "");
}
std::string readObject::getOwnSampName() const {
  VecStr toks = tokenizeString(seqBase_.name_, ".");
  return replaceString(toks[0], "CHI_", "");
}

/// chagening the name of the cluster
void readObject::updateName() {
  size_t totalPos = seqBase_.name_.rfind("_t");
  if (totalPos != std::string::npos) {
    seqBase_.name_ = seqBase_.name_.substr(0, totalPos);
  }
  seqBase_.name_ += "_t" + to_string(seqBase_.cnt_);
}

void readObject::setName(const std::string& newName) {
  seqBase_.name_ = newName + "_t" + to_string(seqBase_.cnt_);
  sampName = getOwnSampName();
}
void readObject::setFractionName() {
  size_t tPos = seqBase_.name_.rfind("_t");
  size_t fPos = seqBase_.name_.rfind("_f");
  if (tPos == std::string::npos && fPos == std::string::npos) {

  } else if (tPos == std::string::npos) {
    seqBase_.name_ = seqBase_.name_.substr(0, fPos);
  } else {
    seqBase_.name_ = seqBase_.name_.substr(0, tPos);
  }
  seqBase_.name_ += "_f" + std::to_string(seqBase_.frac_);
}
std::string readObject::getReadId() const {
  size_t firstPeriod = seqBase_.name_.find(".");
  size_t firstUnder = seqBase_.name_.find("_", firstPeriod);
  return seqBase_.name_.substr(firstPeriod + 1, firstUnder - firstPeriod - 1);
}
std::string readObject::getStubName(bool removeChiFlag) const {
  size_t tPos = seqBase_.name_.rfind("_t");
  size_t fPos = seqBase_.name_.rfind("_f");
  std::string outString = seqBase_.name_;
  if (tPos == std::string::npos && fPos == std::string::npos) {
    outString = seqBase_.name_;
  } else if (tPos != std::string::npos && fPos != std::string::npos) {
    if (tPos < fPos) {
      outString = seqBase_.name_.substr(0, tPos);
    } else {
      outString = seqBase_.name_.substr(0, fPos);
    }
  } else if (tPos == std::string::npos) {
    outString = seqBase_.name_.substr(0, fPos);
  } else {
    outString = seqBase_.name_.substr(0, tPos);
  }

  if (removeChiFlag) {
    outString = replaceString(outString, "CHI_", "");
  }
  return outString;
}

void readObject::setCumulativeFractionName() {
  size_t tPos = seqBase_.name_.rfind("_t");
  size_t fPos = seqBase_.name_.rfind("_f");
  if (tPos == std::string::npos && fPos == std::string::npos) {

  } else if (tPos == std::string::npos) {
    seqBase_.name_ = seqBase_.name_.substr(0, fPos);
  } else {
    seqBase_.name_ = seqBase_.name_.substr(0, tPos);
  }
  seqBase_.name_ += "_f" + std::to_string(cumulativeFraction);
}
void readObject::setNormalizedFractionName() {
  size_t tPos = seqBase_.name_.rfind("_t");
  size_t fPos = seqBase_.name_.rfind("_f");
  if (tPos == std::string::npos && fPos == std::string::npos) {

  } else if (tPos == std::string::npos) {
    seqBase_.name_ = seqBase_.name_.substr(0, fPos);
  } else {
    seqBase_.name_ = seqBase_.name_.substr(0, tPos);
  }
  seqBase_.name_ += "_f" + std::to_string(normalizedFraction);
}

void readObject::appendName(const std::string& add) {
  size_t totalPos = seqBase_.name_.rfind("_t");
  if (totalPos != std::string::npos) {
    seqBase_.name_ = seqBase_.name_.substr(0, totalPos);
  }
  seqBase_.name_ += add + "_t" + to_string(seqBase_.cnt_);
}

void readObject::processRead(bool processed) {
  bool setFraction = false;
  if (processed) {
    VecStr toks;
    bool containsAllNumbers = true;
    if (seqBase_.name_.find("_t") == std::string::npos &&
        seqBase_.name_.find("_f") == std::string::npos) {
      if (seqBase_.name_.rfind("_") == std::string::npos) {
        std::cout << "Improper name format for processed read, should have a "
                     "_# or _t# where # is the number of reads the sequence "
                     "represents" << std::endl;
        std::cout << "failed due to name not have a _ or _t, " << seqBase_.name_
                  << std::endl;
        exit(1);
      } else {
        toks = tokenizeString(seqBase_.name_, "_");
      }
    } else if (seqBase_.name_.find("_t") != std::string::npos) {
      toks = tokenizeString(seqBase_.name_, "_t");
    } else {
      toks = tokenizeString(seqBase_.name_, "_f");
      setFraction = true;
    }
    containsAllNumbers = stringContainsAllDigitsDouble(toks[toks.size() - 1]);
    if (containsAllNumbers) {
      seqBase_.cnt_ = std::stod(toks[toks.size() - 1]);
      if (setFraction) {
        seqBase_.frac_ = seqBase_.cnt_;
        seqBase_.cnt_ = seqBase_.frac_ * 1000;
      }
      seqBase_.name_ = seqBase_.name_.substr(0, seqBase_.name_.rfind("_")) +
                       "_t" + to_string(seqBase_.cnt_);
    } else {
      std::cout << "Improper name format for processed read, should have a _# "
                   "or _t# where # is the number of reads the sequence "
                   "represents" << std::endl;
      std::cout << "failed due to # containing a non-digit character, "
                << toks[toks.size() - 1] << std::endl;
      exit(1);
      // error
    }
  } else {
    // seqBase_.cnt_ = seqBase.cnt_;
  }
}
void readObject::initializeOthers() {
  flowValues.clear();
  seqClip = "";
  qualityClip.clear();
  condensedSeq = "";
  condensedSeqCount.clear();
  averageErrorRate = getAverageErrorRate();
  remove = false;
  cumulativeFraction = 0;
  normalizedFraction = 0;
  numberOfFlows = 0;
  basesAboveQualCheck_ = 0;
  fractionAboveQualCheck_ = 0;
  sampName = getOwnSampName();
}

void readObject::addQual(std::string qualString) {
  seqBase_.qual_ = stringToVector<uint32_t>(qualString);
  averageErrorRate = getAverageErrorRate();
}

void readObject::addQual(std::string qualString, int offSet) {
  seqBase_.qual_.clear();
  for (size_t i = 0; i < qualString.size(); ++i) {
    seqBase_.qual_.push_back(qualString[i] - offSet);
  }
  averageErrorRate = getAverageErrorRate();
}

void readObject::clearObject() {
  seqBase_.name_ = "";
  seqBase_.seq_ = "";
  seqBase_.qual_.clear();
  flowValues.clear();
  seqClip = "";
  qualityClip.clear();
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
  bool print = false;
  if (seqBase_.name_.find("WGHHT:02045:02871_t46.000000") !=
      std::string::npos) {
    print = false;
  }
  uint32_t i = 1;
  for (; i < seqBase_.seq_.length(); i++) {
    if (seqBase_.seq_[i] == seqBase_.seq_[i - 1]) {
      currentQuals.push_back(seqBase_.qual_[i]);
      ++currentCount;
    } else {
      condensedSeq.push_back(seqBase_.seq_[i - 1]);
      condensedSeqQual.push_back(vectorMean(currentQuals));
      currentQualPos.second = len(currentQuals);
      condensedSeqQualPos.emplace_back(currentQualPos);
      currentQualPos.first = i;
      condensedSeqCount.push_back(currentCount);
      if (print) {
        std::cout << "Base: " << seqBase_.seq_[i - 1] << std::endl;
        std::cout << "\t" << vectorToString(currentQuals) << std::endl;
        std::cout << "\t" << currentCount << std::endl;
      }
      currentCount = 1;
      currentQuals.clear();
      currentQuals.push_back(seqBase_.qual_[i]);
    }
  }
  if (print) {
    std::cout << "Base: " << seqBase_.seq_[i - 1] << std::endl;
    std::cout << "\t" << vectorToString(currentQuals) << std::endl;
    std::cout << "\t" << currentCount << std::endl;
  }
  condensedSeq.push_back(seqBase_.seq_[i - 1]);
  condensedSeqQual.push_back(vectorMean(currentQuals));
  currentQualPos.second = len(currentQuals);
  condensedSeqQualPos.emplace_back(currentQualPos);
  condensedSeqCount.push_back(currentCount);
  if (print) {
    std::cout << condensedSeq << std::endl;
    std::cout << vectorToString(condensedSeqQual) << std::endl;
    std::cout << vectorToString(condensedSeqCount) << std::endl;
    std::cout << seqBase_.seq_ << std::endl;
    std::cout << vectorToString(seqBase_.qual_) << std::endl;
    std::cout << seqBase_.seq_.size() << std::endl;
    std::cout << seqBase_.qual_.size() << std::endl;
    checkSeqQual(std::cout);
    exit(1);
  }
}

void readObject::setClip(size_t leftPos, size_t rightPos) {
  seqBase_.seq_ = seqBase_.seq_.substr(leftPos, rightPos - leftPos + 1);
  seqBase_.qual_.erase(seqBase_.qual_.begin() + rightPos + 1,
                       seqBase_.qual_.end());
  seqBase_.qual_.erase(seqBase_.qual_.begin(),
                       seqBase_.qual_.begin() + leftPos);
  averageErrorRate = getAverageErrorRate();
}
void readObject::setClip(size_t rightPos) { setClip(0, rightPos); }
void readObject::setClip(const std::pair<int, int>& positions) {
  setClip(positions.first, positions.second);
}

void readObject::trimFront(size_t upToPosNotIncluding) {
  setClip(upToPosNotIncluding, seqBase_.seq_.size() - 1);
  /*
  seq=seqBase_.seq_.substr(upToPosNotIncluding);
  seqBase_.qual_.erase(qual.begin(),qual.begin()+upToPosNotIncluding-1);
  */
  averageErrorRate = getAverageErrorRate();
}
void readObject::trimBack(size_t fromPositionIncluding) {
  setClip(0, fromPositionIncluding - 1);
  /*
  seq=seqBase_.seq_.substr(0,fromPositionIncluding);
  seqBase_.qual_.erase(qual.begin()+fromPositionIncluding,qual.end());*/
  averageErrorRate = getAverageErrorRate();
}

double readObject::getAverageQual() const {

  double sum = 0;
  for (const auto& q : seqBase_.qual_) {
    sum += q;
  }
  return sum / seqBase_.qual_.size();
}

double readObject::getAverageErrorRate() const {
  double sum = 0;
  for (const auto& q : seqBase_.qual_) {
    sum += pow(10.0, -(q / 10.0));
  }

  return sum / seqBase_.qual_.size();
}

double readObject::getSumQual() const {
  double sum = 0;
  for (const auto& q : seqBase_.qual_) {
    sum += q;
  }
  return sum;
}

void readObject::setFractionByCount(size_t totalNumberOfReads) {
  seqBase_.frac_ = (double)seqBase_.cnt_ / totalNumberOfReads;
  cumulativeFraction = seqBase_.frac_;
}

void readObject::setNormalizedFraction(size_t totalNumberOfSamples) {
  normalizedFraction = seqBase_.frac_ / totalNumberOfSamples;
}

void readObject::setNormalizedFractionByCumulativeFraction(
    size_t totalNumberOfSamples) {
  normalizedFraction = cumulativeFraction / totalNumberOfSamples;
}

std::vector<size_t> readObject::findSubsequenceOccurences(
    const std::string& search) const {
  return findOccurences(seqBase_.seq_, search);
}

void readObject::outPutFlows(std::ostream& flowsFile) const {
  flowsFile << ">" << seqBase_.name_ << std::endl;
  flowsFile << vectorToString(flowValues) << std::endl;
}

void readObject::outPutFlowsRaw(std::ostream& out) const {
  out << ">" << seqBase_.name_ << std::endl;
  out << numberOfFlows << " " << vectorToString(flowValues) << std::endl;
}

void readObject::outPutPyroData(std::ostream& pyroNoiseFile) const {
  pyroNoiseFile << seqBase_.name_ << " " << numberOfFlows << " "
                << vectorToString(flowValues) << std::endl;
}

void readObject::outPutSeqClip(std::ostream& out) const {
  out << ">" << seqBase_.name_ << std::endl;
  out << seqClip << std::endl;
}

void readObject::outPutQualClip(std::ostream& out) const {
  out << ">" << seqBase_.name_ << std::endl;
  out << vectorToString(qualityClip) << std::endl;
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
  for (int i = 0; i < (int)seqBase_.seq_.length(); i++) {
    outFile << seqBase_.seq_[i] << ":" << seqBase_.qual_[i] << " ";
  }
  outFile << std::endl;
}

// counter setters
void readObject::setLetterCount() {
	counter_.reset();
	counter_.increaseCountByString(seqBase_.seq_, seqBase_.cnt_);
	counter_.setFractions();
  //counter_ = letterCounter(seqBase_.seq_, seqBase_.qual_);
}
void readObject::setLetterCount(const std::vector<char> & alph){
	counter_ = charCounterArray(alph);
	counter_.increaseCountByString(seqBase_.seq_, seqBase_.cnt_);
	counter_.setFractions();
}
void readObject::setCondensedCounter() {
  condensedCounter = letterCounter(condensedSeq);
}

// comparisons
bool readObject::operator<(const readObject& otherRead) const {
  if (seqBase_.cnt_ > otherRead.seqBase_.cnt_) {
    return true;
  } else if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
    if (averageErrorRate < otherRead.averageErrorRate) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

bool readObject::operator>(const readObject& otherRead) const {
  if (seqBase_.cnt_ < otherRead.seqBase_.cnt_) {
    return true;
  } else if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
    if (averageErrorRate > otherRead.averageErrorRate) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

bool readObject::operator==(const readObject& otherRead) const {
  if (seqBase_.cnt_ == otherRead.seqBase_.cnt_) {
    return true;
  } else {
    return false;
  }
}
bool readObject::operator>=(const readObject& otherRead) const {
  if (seqBase_.cnt_ <= otherRead.seqBase_.cnt_) {
    return true;
  } else {
    return false;
  }
}

bool readObject::operator<=(const readObject& otherRead) const {
  if (seqBase_.cnt_ >= otherRead.seqBase_.cnt_) {
    return true;
  } else {
    return false;
  }
}

void readObject::setBaseCountOnQualCheck(uint32_t qualCheck) {
  basesAboveQualCheck_ = count_if(seqBase_.qual_, [&](const uint32_t & q){ return q>=qualCheck;});
  fractionAboveQualCheck_ = (double)basesAboveQualCheck_ / seqBase_.qual_.size();
}

void readObject::replace(const std::string& toBeReplaced,
                         const std::string& replaceWith, bool allOccurences) {
  std::vector<size_t> occurences = findSubsequenceOccurences(toBeReplaced);
  std::reverse(occurences.begin(), occurences.end());
  for (std::vector<size_t>::iterator oIter = occurences.begin();
       oIter != occurences.end(); ++oIter) {
    std::vector<uint32_t> currentQuals;
    for (size_t i = *oIter; i < *oIter + toBeReplaced.size(); ++i) {
      currentQuals.push_back(seqBase_.qual_[i]);
    }
    for (size_t i = 0; i < toBeReplaced.size(); ++i) {
      seqBase_.seq_.erase(seqBase_.seq_.begin() + *oIter);
      seqBase_.qual_.erase(seqBase_.qual_.begin() + *oIter);
    }
    double meanQual = vectorMean(currentQuals);
    seqBase_.qual_.insert(seqBase_.qual_.begin() + *oIter, replaceWith.size(),
                          meanQual);
    seqBase_.seq_.insert(*oIter, replaceWith.c_str());
  }
}

uint32_t readObject::getNumberOfBasesFromFlow(
    const std::vector<double>& inFlows) {
  uint32_t numOfBases = 0;
  for (auto fIter = inFlows.begin(); fIter != inFlows.end(); ++fIter) {
    if (fIter == inFlows.begin() && *fIter > 20) {
      continue;
    }
    if ((*fIter) > 0.7) {
      // non noisy signal increase bases count by nearest integer
      numOfBases += std::round(*fIter);
    }
  }
  return numOfBases;
}
int readObject::clipToNinetyPercentOfFlows(size_t cutOff) {
  processedFlowValues = flowValues;
  if (processedFlowValues.size() > cutOff) {
    processedFlowValues.erase(processedFlowValues.begin() + cutOff,
                              processedFlowValues.end());
  }
  std::cout << "Clipped to : " << processedFlowValues.size() << std::endl;
  int clipTo = getNumberOfBasesFromFlow(flowValues);
  // setClip(0, clipTo);
  return clipTo;
}

size_t readObject::findFirstNoisyFlow(const std::vector<double>& inFlows) {
  int noise = 0;
  int signal = 0;
  int basesSinceSignal = 0;
  size_t pos = 0;
  for (std::vector<double>::const_iterator fIter = inFlows.begin();
       fIter != inFlows.end(); ++fIter) {

    if (*fIter > 0.5) {
      if (*fIter < 0.7) {
        ++noise;
      } else {
        ++signal;
        basesSinceSignal = 0;
      }
    } else {
      ++basesSinceSignal;
    }
    if (basesSinceSignal >= 4 || noise >= 1) {
      break;
    }
    ++pos;
  }
  return pos;
}

bool readObject::flowNoiseProcess(size_t cutoff) {
  size_t firstNoiseReadPos = findFirstNoisyFlow(flowValues);
  if (firstNoiseReadPos < 360) {
    return false;
  }
  if (firstNoiseReadPos < cutoff) {
    cutoff = firstNoiseReadPos;
  }
  processedFlowValues = flowValues;
  processedFlowValues.erase(processedFlowValues.begin() + cutoff,
                            processedFlowValues.end());

  auto numberOfBases = getNumberOfBasesFromFlow(processedFlowValues);
  if (numberOfBases < seqBase_.seq_.size()) {
    setClip(0, numberOfBases);
  }

  return true;
}

void readObject::convertToProteinFromcDNA(bool transcribeToRNAFirst,
                                          size_t start, bool forceStartM) {
  if (transcribeToRNAFirst) {
    seqUtil::convertToProteinFromcDNA(seqBase_.seq_, start, forceStartM);
  } else {
    seqBase_.seq_ =
        seqUtil::convertToProtein(seqBase_.seq_, start, forceStartM);
  }
}
std::string readObject::getProteinFromcDNA(bool transcribeToRNAFirst,
                                           size_t start,
                                           bool forceStartM) const {
  if (transcribeToRNAFirst) {
    return seqUtil::convertToProteinFromcDNAReturn(seqBase_.seq_, start,
                                                   forceStartM);
  } else {
    return seqUtil::convertToProtein(seqBase_.seq_, start, forceStartM);
  }
}
double readObject::getGCContent() {
  setLetterCount();
  counter_.calcGcContent();
  return counter_.gcContent;
}

void readObject::adjustHomopolyerRunQualities() {
  createCondensedSeq();
  seqBase_.qual_.clear();
  for (const auto& i : iter::range<uint64_t>(0, condensedSeq.length())) {
    addOtherVec(seqBase_.qual_, std::vector<uint32_t>(condensedSeqCount[i],
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

  for (const auto& pos : iter::range(len(seqBase_.qual_))) {
    updateQaulCountsAtPos(pos, counts, qualWindowSize, qualErrorLookUp);
  }
}
void readObject::updateQualWindowCounts(
    uint32_t pos, std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize) const {
  std::vector<uint32_t> currentQuals;
  currentQuals.push_back(seqBase_.qual_[pos]);
  addOtherVec(currentQuals, seqBase_.getLeadQual(pos, qualWindowSize));
  addOtherVec(currentQuals, seqBase_.getTrailQual(pos, qualWindowSize));
  counts["mean"][roundDecPlaces(vectorMean(currentQuals), 2)] += seqBase_.cnt_;
  counts["median"][roundDecPlaces(vectorMedian(currentQuals), 2)] +=
      seqBase_.cnt_;
  counts["min"][roundDecPlaces(vectorMinimum(currentQuals), 2)] +=
      seqBase_.cnt_;
  return;
}
void readObject::updateQualWindowCounts(
    std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize) const {
  for (const auto& pos : iter::range(seqBase_.qual_.size())) {
    updateQualWindowCounts(pos, counts, qualWindowSize);
  }
}

void readObject::updateBaseQualCounts(std::map<double, uint32_t>& baseCounts,
                                      uint32_t pos) const {
  baseCounts[roundDecPlaces(seqBase_.qual_[pos], 2)] += seqBase_.cnt_;
}

void readObject::updateBaseQualCounts(std::map<double, uint32_t>& baseCounts)
    const {
  for (const auto& pos : iter::range(seqBase_.qual_.size())) {
    updateBaseQualCounts(baseCounts, pos);
  }
}
void readObject::updateQaulCountsAtPos(
    uint32_t pos, std::map<std::string, std::map<double, uint32_t>>& counts,
    int qualWindowSize, const std::array<double, 100>& qualErrorLookUp) const {

  std::vector<uint32_t> currentQuals;
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

void readObject::setQualMeans(uint32_t windowSize) {
  qualMeans_.clear();
  for (auto i : iter::range(len(seqBase_.qual_))) {
    uint32_t lowerBound = 0;
    uint32_t higherBound = len(seqBase_.qual_);
    if (i > windowSize) {
      lowerBound = i - windowSize;
    }
    if (i + windowSize + 1 < higherBound) {
      higherBound = i + windowSize + 1;
    }
    double sum = 0.0;
    for (auto cursor : iter::range(lowerBound, higherBound)) {
      sum += seqBase_.qual_[cursor];
    }
    qualMeans_.push_back(sum / (higherBound - lowerBound));
  }
}
void readObject::setQualProbMeans(uint32_t windowSize,
                                  const std::array<double, 100>& errorLookUp) {
  qualProbMeans_.clear();
  for (auto i : iter::range(len(seqBase_.qual_))) {
    uint32_t lowerBound = 0;
    uint32_t higherBound = len(seqBase_.qual_);
    if (i > windowSize) {
      lowerBound = i - windowSize;
    }
    if (i + windowSize + 1 < higherBound) {
      higherBound = i + windowSize + 1;
    }
    double sum = 0.0;
    for (auto cursor : iter::range(lowerBound, higherBound)) {
      sum += errorLookUp[seqBase_.qual_[cursor]];
    }
    qualProbMeans_.push_back(sum / (higherBound - lowerBound));
  }
}

void readObject::printDescription(std::ostream& out, bool deep) const {
  baseReadObject::printDescription(out);
  out << "readObject{" << std::endl << "qualMeans_:" << qualMeans_ << std::endl
      << "qualProbMeans_:" << qualProbMeans_ << std::endl
      << "flowValues:" << flowValues << std::endl
      << "processedFlowValues:" << processedFlowValues << std::endl
      << "seqClip:" << seqClip << std::endl << "qualityClip:" << qualityClip
      << std::endl << "sampName:" << sampName << std::endl
      << "expectsString:" << expectsString << std::endl
      << "numberOfFlows:" << numberOfFlows << std::endl
      << "cumulativeFraction:" << cumulativeFraction << std::endl
      << "normalizedFraction:" << normalizedFraction << std::endl
      << "averageErrorRate:" << averageErrorRate << std::endl
      << "basesAboveQualCheck:" << basesAboveQualCheck_ << std::endl
      << "fractionAboveQualCheck:" << fractionAboveQualCheck_ << std::endl
      << "remove:" << remove << std::endl;
  if (deep) {
    counter_.printDescription(out, deep);
  }
  out << "condensedSeq:" << condensedSeq << std::endl
      << "condensedSeqQual:" << condensedSeqQual << std::endl
      << "condensedSeqCount:" << condensedSeqCount << std::endl;
  if (deep) {
    condensedCounter.printDescription(out, deep);
  }
  out << "}" << std::endl;
}

}  // namespace bib