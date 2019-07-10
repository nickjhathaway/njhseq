#include "seqInfo.hpp"
#include "njhseq/helpers/seqUtil.hpp"
#include <njhcpp/bashUtils.h>

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

seqInfo::seqInfo() :
		name_(""), seq_(""), cnt_(1), frac_(0) {
}
seqInfo::seqInfo(const std::string & name) :
		name_(name), seq_(""), cnt_(1), frac_(0) {
}
seqInfo::seqInfo(const std::string& name, const std::string& seq,
		const std::vector<uint32_t>& qual) :
		name_(name), seq_(seq), qual_(qual), cnt_(1), frac_(0) {
	if(qual_.size() != seq.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error seq and qual have to be the same size; seq: " << seq.size() << " qual: " << qual_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}
}
seqInfo::seqInfo(const std::string& name, const std::string& seq,
		const std::vector<uint32_t>& qual, double cnt) :
		name_(name), seq_(seq), qual_(qual), cnt_(cnt), frac_(0) {
	if(qual_.size() != seq.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error seq and qual have to be the same size; seq: " << seq.size() << " qual: " << qual_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}
}
seqInfo::seqInfo(const std::string& name, const std::string& seq) :
		name_(name), seq_(seq), qual_(std::vector<uint32_t>(seq.size(), 40)), cnt_(
				1), frac_(0) {

}

//seqInfo::seqInfo(const std::string& name, const std::string& seq,
//		const std::string& stringQual) :
//		name_(name), seq_(seq), qual_(stringToVector<uint32_t>(stringQual)), cnt_(
//				1), frac_(0) {
//}

seqInfo::seqInfo(const std::string& name, const std::string& seq,
		const std::string& stringQual, uint32_t off_set) :
		name_(name), seq_(seq), qual_(std::vector<uint32_t>(0)), cnt_(1), frac_(0) {
	for (const auto & c : stringQual) {
		qual_.emplace_back(c - off_set);
	}
	if(qual_.size() != seq.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error seq and qual have to be the same size; seq: " << seq.size() << " qual: " << qual_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}
}
seqInfo::seqInfo(const std::string& name, const std::string& seq,
		const std::vector<uint32_t>& qual, double cnt, double frac) :
		name_(name), seq_(seq), qual_(qual), cnt_(cnt), frac_(frac) {
	if(qual_.size() != seq.size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error seq and qual have to be the same size; seq: " << seq.size() << " qual: " << qual_.size()<< "\n";
		throw std::runtime_error{ss.str()};
	}
}

seqInfo seqInfo::getSubRead(uint32_t pos, uint32_t size) const {
	return seqInfo(name_, seq_.substr(pos, size), getSubVector(qual_, pos, size),
			cnt_, frac_);
}
seqInfo seqInfo::getSubRead(uint32_t pos) const {
	return seqInfo(name_, seq_.substr(pos, seq_.size() - pos),
			getSubVector(qual_, pos, seq_.size() - pos), cnt_, frac_);
}



void seqInfo::updateName() {
	size_t totalPos = name_.rfind("_t");
	if (totalPos != std::string::npos) {
		name_ = name_.substr(0, totalPos);
	}
	name_ += "_t" + estd::to_string(cnt_);
}


Json::Value seqInfo::toJsonJustInfo() const {
	Json::Value ret;
	ret["cnt"] = njh::json::toJson(cnt_);
	ret["frac"] = njh::json::toJson(frac_);
	ret["name"] = njh::json::toJson(name_);
	ret["on"] = njh::json::toJson(on_);
	return ret;
}
Json::Value seqInfo::toJson() const {
	Json::Value ret;
	ret["seq"] = njh::json::toJson(seq_);
	ret["qual"] = njh::json::toJson(qual_);
	ret["cnt"] = njh::json::toJson(cnt_);
	ret["frac"] = njh::json::toJson(frac_);
	ret["name"] = njh::json::toJson(name_);
	ret["on"] = njh::json::toJson(on_);
	return ret;
}

void seqInfo::translate(bool complement, bool reverse, size_t start) {
	if (complement && reverse) {
		reverseComplementRead(true, true);
		seq_ = seqUtil::convertToProtein(seq_, start, false);
	} else if (complement) {
		seqUtil::convertToProteinFromcDNA(seq_, start, false);
	} else {
		seq_ = seqUtil::convertToProtein(seq_, start, false);
	}
	qual_ = std::vector<uint32_t>(seq_.size(), 40);
}

seqInfo seqInfo::translateRet(bool complement, bool reverse,
		size_t start) const {
	seqInfo ret(*this);
	ret.translate(complement, reverse, start);
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << "\t" << ret.seq_ << std::endl;
	return ret;
}

void seqInfo::processRead(bool processed) {
	bool setFraction = false;
	if (processed) {
		VecStr toks;
		bool containsAllNumbers = true;
		if (name_.rfind("_t") == std::string::npos
				&& name_.rfind("_f") == std::string::npos) {
			if (name_.rfind("_") == std::string::npos) {
				std::stringstream ss;
				ss << "Improper name format for processed read, should have a "
						"_#, _t#, or _f# where # is the number of reads the sequence "
						"represents" << "\n";
				ss << "failed due to name not have a _ or _t, " << name_ << "\n";
				throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
			} else {
				toks = tokenizeString(name_, "_");
			}
		}else if (name_.rfind("_t") != std::string::npos
				&& name_.rfind("_f") != std::string::npos){
			if(name_.rfind("_t") > name_.rfind("_f")){
				toks = tokenizeString(name_, "_t");
			}else{
				toks = tokenizeString(name_, "_f");
				setFraction = true;
			}
		} else if (name_.rfind("_t") != std::string::npos) {
			toks = tokenizeString(name_, "_t");
		} else {
			toks = tokenizeString(name_, "_f");
			setFraction = true;
		}
		containsAllNumbers = isDoubleStr(toks[toks.size() - 1]);
		if (containsAllNumbers) {
			cnt_ = std::stod(toks[toks.size() - 1]);
			if (setFraction) {
				frac_ = cnt_;
				cnt_ = frac_ * 1000;
			}
			//name_ = name_.substr(0, name_.rfind("_")) +
			//                 "_t" + estd::to_string(cnt_);
		} else {
			std::stringstream ss;
			ss << "Improper name format for processed read, should have a _# "
					"or _t# where # is the number of reads the sequence "
					"represents" << "\n";
			ss << "failed due to # containing a non-digit character, "
					<< toks[toks.size() - 1] << " from " << name_<< "\n";
			throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
		}
	} else {
		//leave the count that it was originally when constructed
		//cnt_ = 1;
	}
}

const std::unordered_map<char, uint32_t> seqInfo::ansiBaseColor =
		std::unordered_map<char, uint32_t> { { 'A', 203 }, { 'a', 203 },
				{ 'C', 83 }, { 'c', 83 }, { 'G', 227 }, { 'g', 227 }, { 'T', 69 }, {
						't', 69 }, { 'n', 145 }, { 'N', 145 }, { '-', 102 } };

void seqInfo::outPutSeqAnsi(std::ostream& fastaFile) const {
	fastaFile << njh::bashCT::addBGColor(145) << ">" << name_
			<< njh::bashCT::reset << "\n";
	for (const auto & c : seq_) {
		if(ansiBaseColor.find(c) == ansiBaseColor.end()){
			fastaFile << njh::bashCT::addBGColor(16) << c
					<< njh::bashCT::reset;
		}else{
			fastaFile << njh::bashCT::addBGColor(ansiBaseColor.at(c)) << c
					<< njh::bashCT::reset;
		}
	}
	fastaFile << njh::bashCT::reset << "\n";
}
void seqInfo::prepend(const std::string& seq, uint32_t defaultQuality) {
	prepend(seq, std::vector<uint32_t>(1, defaultQuality));
}
void seqInfo::append(const std::string& seq, uint32_t defaultQuality) {
	append(seq, std::vector<uint32_t>(1, defaultQuality));
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
		std::stringstream ss;
		ss << "Need to supply either single a quality or same amount of "
				"quality score for seq length" << "\n";
		ss << "trying to add " << qual.size() << " qualities and a seq of length "
				<< seq.length() << "\n";
		throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
	}
}

void seqInfo::prepend(const char & base, uint32_t quality) {
	seq_.insert(seq_.begin(), base);
	qual_.insert(qual_.begin(), quality);
}

void seqInfo::append(const char & base, uint32_t quality) {
	seq_.push_back(base);
	qual_.emplace_back(quality);
}

void seqInfo::insert(uint32_t pos, const seqInfo & otherInfo){
	if(pos >= len(*this)){
		std::stringstream ss;
		ss << __FILE__ << ": " << __LINE__ << " : " << __PRETTY_FUNCTION__ << ", error out of range insert" << "\n";
		ss << "Position " << pos << " is greater than or equal to the length of the current sequence " << len(*this) << "\n";
		throw std::runtime_error{ss.str()};
	}
	seq_.insert(seq_.begin() + pos, otherInfo.seq_.begin(), otherInfo.seq_.end());
	qual_.insert(qual_.begin() + pos, otherInfo.qual_.begin(), otherInfo.qual_.end());
}

void seqInfo::reverseHRunsQuals() {
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
	for (auto hrQuals : quals) {
		njh::reverse(hrQuals);
		addOtherVec(qual_, hrQuals);
	}
}

void seqInfo::reverseComplementRead(bool mark, bool regQualReverse) {
	if (regQualReverse) {
		njh::reverse(qual_);
	} else {
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
		for (auto iter = quals.rbegin(); iter != quals.rend(); ++iter) {
			addOtherVec(qual_, *iter);
		}
	}
	seq_ = seqUtil::reverseComplement(seq_, "DNA");

	if (mark) {
		//if name already contains _Comp and then remove _Comp
		if (name_.find("_Comp") != std::string::npos) {
			name_ = njh::replaceString(name_, "_Comp", "");
		} else {
			name_.append("_Comp");
		}
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
		std::stringstream ss;
		ss << "Need to supply either single a quality or same amount of "
				"quality score for seq length" << "\n";
		ss << "trying to add " << qual.size() << " qualities and a seq of length "
				<< seq.length() << "\n";
		throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
	}
}

void seqInfo::prepend(const seqInfo & other){
	prepend(other.seq_, other.qual_);
}

void seqInfo::append(const seqInfo & other){
	append(other.seq_, other.qual_);
}


const std::vector<uint32_t> seqInfo::getLeadQual(uint32_t posA,
		uint32_t out) const {
	std::vector<uint32_t> ans;
	uint32_t lowerBound = 0;
	if (posA - out > lowerBound) {
		lowerBound = posA - out;
	}
	for (auto i : iter::range(lowerBound, posA)) {
		ans.emplace_back(qual_[i]);
	}
	return ans;
}
const std::vector<uint32_t> seqInfo::getTrailQual(uint32_t posA,
		uint32_t out) const {
	std::vector<uint32_t> ret;
	uint32_t higherBound = qual_.size() - 1;
	if (posA + out + 1 < higherBound) {
		higherBound = posA + out + 1;
	}
	for (auto i : iter::range(posA + 1, higherBound)) {
		ret.emplace_back(qual_[i]);
	}
	return ret;
}
bool seqInfo::checkLeadQual(uint32_t pos, uint32_t secondayQual,
		uint32_t out) const {
	uint32_t lowerBound = 0;
	if (static_cast<int32_t>(pos) - out > lowerBound) {
		lowerBound = pos - out;
	}
	for (auto i : iter::range(lowerBound, pos)) {
		if (qual_[i] <= secondayQual) {
			return false;
		}
	}
	return true;
}
bool seqInfo::checkTrailQual(uint32_t pos, uint32_t secondayQual,
		uint32_t out) const {
	uint32_t higherBound = qual_.size() - 1;
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
bool seqInfo::checkPrimaryQual(uint32_t pos, uint32_t primaryQual) const {
	if (qual_[pos] <= primaryQual) {
		return false;
	}
	return true;
}
bool seqInfo::checkQual(uint32_t pos, const QualScorePars & qScorePars) const {
	if (!checkPrimaryQual(pos, qScorePars.primaryQual_)) {
		return false;
	}
	if(0 != qScorePars.qualThresWindow_){
		if (!checkLeadQual(pos, qScorePars.secondaryQual_, qScorePars.qualThresWindow_)
				|| !checkTrailQual(pos, qScorePars.secondaryQual_, qScorePars.qualThresWindow_)) {
			return false;
		}
	}
	return true;
}
uint32_t seqInfo::findLowestNeighborhoodQual(uint32_t posA,
		uint32_t out) const {
	uint32_t lowerBound = 0;
	uint32_t higherBound = qual_.size() - 1;
	if (static_cast<int32_t>(posA) - out > lowerBound) {
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
		std::stringstream ss;
		ss << "pos: " << pos << " out of bounds of seq " << seq_.size() << "\n";
		throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
	} else if (pos >= qual_.size()) {
		std::stringstream ss;
		ss << "pos: " << pos << " out of bounds of qual " << qual_.size() << "\n";
		throw std::runtime_error { njh::bashCT::boldRed(ss.str()) };
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
	for (auto pos : iter::range<int32_t>(seq_.size(), -1, -1)) {
		if ('-' == seq_[pos] || '.' == seq_[pos]) {
			removeBase(pos);
		}
	}
	return;
}

// quality strings for printing
std::string seqInfo::getQualString() const {
	return vectorToString(qual_);
}

void seqInfo::setFractionByCount(double totalNumberOfReads){
	frac_ = cnt_/totalNumberOfReads;
}

std::string seqInfo::getFastqString(const std::vector<uint32_t>& quals,
		uint32_t offset) {
	std::string convertedQuals = "";
	for (const auto& q : quals) {
		if (q <= 93) {
			convertedQuals.push_back(static_cast<char>(q + offset));
		} else {
			convertedQuals.push_back(static_cast<char>(93 + offset));
		}
	}
	return convertedQuals;
}
std::string seqInfo::getFastqQualString(uint32_t offset) const {
	return getFastqString(qual_, offset);
}

//
void seqInfo::markAsChimeric() {
	if (name_.find("CHI") == std::string::npos) {
		name_ = "CHI_" + name_;
	}
}

void seqInfo::unmarkAsChimeric() {
	if (name_.find("CHI") != std::string::npos) {
		name_ = njh::replaceString(name_, "CHI_", "");
	}
}

// outputs
void seqInfo::outPutFastq(std::ostream& fastqFile) const {
	fastqFile << "@" << name_ << "\n";
	fastqFile << seq_ << "\n";
	fastqFile << "+" << "\n";
	fastqFile << getFastqQualString(SangerQualOffset) << "\n";
}

void seqInfo::outPutSeq(std::ostream& fastaFile) const {
	fastaFile << ">" << name_ << "\n";
	fastaFile << seq_ << "\n";
}

void seqInfo::outPutQual(std::ostream& qualFile) const {
	qualFile << ">" << name_ << "\n";
	qualFile << getQualString() << "\n";
}
/*
void seqInfo::outPut(std::ostream& outFile,
		const SeqIOOptions & options) const {
	if (options.out_.outFormat_ == "fastq") {
		outPutFastq(outFile);
	} else if (options.out_.outFormat_ == "fasta") {
		outPutSeq(outFile);
	} else {
		throw std::runtime_error { njh::bashCT::boldRed(
				"in " + std::string(__PRETTY_FUNCTION__) + " : unrecognized option: "
						+ options.out_.outFormat_) };
	}
}

void seqInfo::outPut(std::ostream& outFile, std::ostream& outFile2,
		const SeqIOOptions & options) const {
	if (options.out_.outFormat_ == "fastqQual") {
		outPutSeq(outFile);
		outPutQual(outFile2);
	} else {
		throw std::runtime_error { njh::bashCT::boldRed(
				"in " + std::string(__PRETTY_FUNCTION__) + " : unrecognized option: "
						+ options.out_.outFormat_) };
	}
}*/

bool seqInfo::degenCompare(const seqInfo & otherInfo,
		const substituteMatrix & compareScores) const {
	if (seq_.length() != otherInfo.seq_.length()) {
		return false;
	}
	auto comp = [&](const char & a, const char & b) {
		return compareScores.mat_[a][b] > 0;
	};
	return std::equal(seq_.begin(), seq_.end(), otherInfo.seq_.begin(), comp);
}

std::string seqInfo::getReadId() const {
  size_t firstPeriod = name_.find(".");
  size_t firstUnder = name_.find("_", firstPeriod);
  return name_.substr(firstPeriod + 1, firstUnder - firstPeriod - 1);
}

std::string seqInfo::getStubName(bool removeChiFlag) const {
	size_t tPos = name_.rfind("_t");
	size_t fPos = name_.rfind("_f");
	std::string outString = name_;
	if (tPos == std::string::npos && fPos == std::string::npos) {
		outString = name_;
	} else if (tPos != std::string::npos && fPos != std::string::npos) {
		if (tPos > fPos) {
			outString = name_.substr(0, tPos);
		} else {
			outString = name_.substr(0, fPos);
		}
	} else if (tPos == std::string::npos) {
		outString = name_.substr(0, fPos);
	} else {
		outString = name_.substr(0, tPos);
	}

	if (removeChiFlag) {
		outString = njh::replaceString(outString, "CHI_", "");
	}
//	if (removeChiFlag) {
//		if(MetaDataInName::nameHasMetaData(outString)){
//			MetaDataInName::removeMetaDataInName(outString);
//		}
//	}
	return outString;
}

void seqInfo::setName(const std::string& newName) {
	name_ = newName + "_t" + estd::to_string(cnt_);
}

void seqInfo::addQual(const std::string & qualString) {
	addQual(stringToVector<uint32_t>(qualString));
}

void seqInfo::addQual(const std::string & qualString, uint32_t offSet) {
	qual_.clear();
	for (const auto & c : qualString) {
		qual_.emplace_back(c - offSet);
	}
	if (qual_.size() != seq_.size()) {
		std::stringstream ss;
		ss << "adding qual size does not equal seq size, qualSize: " << qual_.size()
				<< ", seqSize: " << seq_.size() << ", for " << name_;
		throw std::runtime_error { ss.str() };
	}
}

void seqInfo::addQual(const std::vector<uint32_t> & quals) {
	if (quals.size() != seq_.size()) {
		std::stringstream ss;
		ss << "adding qual size does not equal seq size, qualSize: " << quals.size()
				<< ", seqSize: " << seq_.size() << ", for " << name_ << "\n";
		printVector(quals, ", ", ss);
		ss << seq_ << "\n";
		throw std::runtime_error { ss.str() };
	}
	qual_ = quals;
}

double seqInfo::getQualCheck(uint32_t qualCutOff) const {
	uint32_t count = njh::count_if(qual_,
			[&qualCutOff](uint32_t qual) {return qual >=qualCutOff;});
	return static_cast<double>(count) / qual_.size();
}

void seqInfo::setClip(size_t upToPosNotIncluding, size_t fromPositionNotIncluding) {
	seq_ = seq_.substr(upToPosNotIncluding, fromPositionNotIncluding - upToPosNotIncluding + 1);
	qual_.erase(qual_.begin() + fromPositionNotIncluding + 1, qual_.end());
	qual_.erase(qual_.begin(), qual_.begin() + upToPosNotIncluding);
}

void seqInfo::clipOut(size_t position, size_t size){
	seq_.erase(seq_.begin() + position, seq_.begin() + position + size);
	qual_.erase(qual_.begin() + position, qual_.begin() + position + size);
}

void seqInfo::trimFront(size_t upToPosNotIncluding) {
	setClip(upToPosNotIncluding, seq_.size() - 1);
}

void seqInfo::trimBack(size_t fromPositionIncluding) {
	setClip(0, fromPositionIncluding - 1);
}

double seqInfo::getAverageQual() const {
	return static_cast<double>(getSumQual()) / qual_.size();
}

/*
 double seqInfo::getAverageErrorRate() const {
 //return 0;
 njh::randomGenerator gen;
 double sum = 0;
 for (auto q : qual_) {
 sum += std::pow(10.0, -(gen.unifRand(0.01, 4.0) / 10.0));
 //sum += std::pow(10.0, -(q / 10.0));
 }
 return sum / qual_.size();
 }*/

double seqInfo::getAverageErrorRate() const {
	double sum = 0;
	for (auto q : qual_) {
		sum += std::pow(10.0, -(q / 10.0));
	}
	return sum / qual_.size();
}

uint32_t seqInfo::getSumQual() const {
	uint32_t sum = 0;
	for (const auto& q : qual_) {
		sum += q;
	}
	return sum;
}

std::string getStubNameExternal(const std::string & name, bool removeChiFlag)  {
	size_t tPos = name.rfind("_t");
	size_t fPos = name.rfind("_f");
	std::string outString = name;
	if (tPos == std::string::npos && fPos == std::string::npos) {
		outString = name;
	} else if (tPos != std::string::npos && fPos != std::string::npos) {
		if (tPos > fPos) {
			outString = name.substr(0, tPos);
		} else {
			outString = name.substr(0, fPos);
		}
	} else if (tPos == std::string::npos) {
		outString = name.substr(0, fPos);
	} else {
		outString = name.substr(0, tPos);
	}

	if (removeChiFlag) {
		outString = njh::replaceString(outString, "CHI_", "");
	}
//	if (removeChiFlag) {
//		if(MetaDataInName::nameHasMetaData(outString)){
//			MetaDataInName::removeMetaDataInName(outString);
//		}
//	}
	return outString;
}

std::string seqInfo::getOwnSampName() const {

	if(MetaDataInName::nameHasMetaData(name_)){
		MetaDataInName meta(name_);
		if(meta.containsMeta("samp")){
			return meta.getMeta("samp");
		}else if(meta.containsMeta("sample")){
			return meta.getMeta("sample");
		}
	}
	//std::cout << name_ << std::endl;
	std::string name;
	if(MetaDataInName::nameHasMetaData(name_)){
		name = name_;
		MetaDataInName::removeMetaDataInName(name);
		name = getStubNameExternal(name, true);
	}else{
		name = getStubName(true);
	}
	//std::cout << name.substr(0,name.rfind(".")) << std::endl;
	return name.substr(0,name.rfind("."));
	/*
	std::string name = name_;
	auto firstBracket = name_.find("[");
	if(std::string::npos != firstBracket){
		name = name.substr(firstBracket);
	}
	VecStr toks = tokenizeString(name, ".");
	return njh::replaceString(toks[0], "CHI_", "");
	*/
}

bool seqInfo::nameHasMetaData() const {
	auto firstBracket = name_.find("[");
	if (std::string::npos == firstBracket) {
		return false;
	}
	auto secondBracket = name_.find("]", firstBracket);
	if (std::string::npos == secondBracket) {
		return false;
	}
	return true;
}

void seqInfo::resetMetaInName(const MetaDataInName & meta) {
	if (meta.meta_.empty()) {
		return;
	}

	if (MetaDataInName::nameHasMetaData(name_)) {
		meta.resetMetaInName(name_);
	} else {
		std::string newMeta = meta.createMetaName();
		//std::cout << name_ << std::endl;
		//std::cout << newMeta << std::endl;
		auto countPat = name_.rfind("_t");
		auto fracPat = name_.rfind("_f");
		if (std::string::npos != countPat && countPat + 2 != name_.length()
				&& isDoubleStr(name_.substr(countPat + 2))) {
			auto rest = name_.substr(countPat + 2);
			name_ = name_.substr(0, countPat) + newMeta + name_.substr(countPat);
		} else if (std::string::npos != fracPat && fracPat + 2 != name_.length()
				&& isDoubleStr(name_.substr(fracPat + 2))) {
			auto rest = name_.substr(fracPat + 2);
			name_ = name_.substr(0, fracPat) + newMeta + name_.substr(fracPat);
		} else {
			name_ = name_ + newMeta;
		}
		//std::cout << name_ << std::endl << std::endl;
	}
}

void seqInfo::processNameForMeta(std::unordered_map<std::string, std::string> & meta)const{
	auto firstBracket = name_.find("[");
	if(std::string::npos == firstBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", could not find [ in " << name_ << std::endl;
		throw std::runtime_error{ss.str()};
	}
	auto secondBracket = name_.find("]", firstBracket);
	if(std::string::npos == secondBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", could not find ] in " << name_  << " after " << firstBracket
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	auto toks = tokenizeString(name_.substr(firstBracket + 1, secondBracket - firstBracket - 1), ";");
	for(const auto & tok : toks){
		auto subToks = tokenizeString(tok, "=");
		if(2 != subToks.size()){
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__
					<< "values should be separated by one =, not " << tok
					<< std::endl;
			throw std::runtime_error{ss.str()};
		}else{
			if(meta.find(subToks[0]) != meta.end()){
				std::stringstream ss;
				ss << "Error in : " << __PRETTY_FUNCTION__ << ", key " << subToks[0]
						<< " already in meta " << std::endl;
				ss << "value is already: " << meta.find(subToks[0])->first
						<< ", attempt to add:  " << subToks[1] << std::endl;
				throw std::runtime_error { ss.str() };
			}else{
				meta.emplace(subToks[0], subToks[1]);
			}
		}
	}
}

bool seqInfo::isChimeric() const {
	return njh::beginsWith(name_, "CHI_");
	//return name_.find("CHI") != std::string::npos;
}

bool seqInfo::operator ==(const seqInfo & other) const{
	return seq_ == other.seq_ &&
			qual_ == other.qual_ &&
			on_ == other.on_ &&
			cnt_ == other.cnt_ &&
			frac_ == other.frac_;
}



void seqInfo::adjustHomopolymerRunQualities() {
  std::string condensedSeq = "";
  std::vector<uint32_t> condensedSeqQual;
  std::vector<uint32_t>condensedSeqCount;
  std::vector<std::pair<uint32_t, uint32_t>>condensedSeqQualPos;
  int currentCount = 1;
  std::vector<uint32_t> currentQuals;
  currentQuals.push_back(qual_[0]);
  std::pair<uint32_t, uint32_t> currentQualPos {0,1};
  uint32_t i = 1;
  for (; i < seq_.length(); i++) {
    if (seq_[i] == seq_[i - 1]) {
      currentQuals.push_back(qual_[i]);
      ++currentCount;
    } else {
      condensedSeq.push_back(seq_[i - 1]);
      condensedSeqQual.push_back(vectorMean(currentQuals));
      currentQualPos.second = currentQuals.size();
      condensedSeqQualPos.emplace_back(currentQualPos);
      currentQualPos.first = i;
      condensedSeqCount.push_back(currentCount);
      currentCount = 1;
      currentQuals.clear();
      currentQuals.push_back(qual_[i]);
    }
  }
  condensedSeq.push_back(seq_[i - 1]);
  condensedSeqQual.push_back(vectorMean(currentQuals));
  currentQualPos.second = currentQuals.size();
  condensedSeqQualPos.emplace_back(currentQualPos);
  condensedSeqCount.push_back(currentCount);
  qual_.clear();
  for (const auto& i : iter::range<uint64_t>(0, condensedSeq.length())) {
    addOtherVec(qual_, std::vector<uint32_t>(condensedSeqCount[i],
                                                      condensedSeqQual[i]));
  }
}



}  // namespace njhseq
