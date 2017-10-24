#include "sampleCluster.hpp"
#include <bibcpp/common.h>
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
namespace bibseq {
sampleCluster::sampleCluster(const seqInfo& initializerRead,
		const std::map<std::string, sampInfo>& infos) {
	seqBase_ = initializerRead;

	firstReadName_ = initializerRead.name_;
	firstReadCount_ = initializerRead.cnt_;

	sampName = getOwnSampName();

	for (const auto & i : infos) {
		sampInfos_[i.second.runName_] = sampInfo(i.second.runName_,
				i.second.runReadCnt_);
	}

	reads_.emplace_back(std::make_shared<readObject>(*this));
	updateInfoWithRead(*this, 0);
	remove = false;
	needToCalculateConsensus_ = false;
	setLetterCount();
}


void sampleCluster::addRead(const sampleCluster& cr) {
  // add the cluster's reads and update the info
  for (const auto& read : cr.reads_) {
    updateInfoWithRead(*read, reads_.size());
    reads_.push_back(read);
  }
  // update the fraction and totalCounts
  seqBase_.cnt_ += cr.seqBase_.cnt_;
  seqBase_.frac_ = getAveragedFrac(); //calculate average as the mean fraction between all samples
  updateName();
  //consider putting in parameters to determine if
  needToCalculateConsensus_ = true;
}

void sampleCluster::updateInfoWithRead(const readObject& read, uint32_t pos) {
  // update infos with the read
	sampleClusters_[read.sampName].push_back(pos);
  sampInfos_[read.sampName].update(read.seqBase_);
}


std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleCluster::getRepDisagreement() const {
	std::unordered_map<std::string, std::unordered_map<std::string, double>> ret;
	for (const auto & sInfo : sampInfos_) {
		for (const auto & secondSInfo : sampInfos_) {
			if (sInfo.first == secondSInfo.first) {
				continue;
			}
			if (0 == sInfo.second.numberOfClusters_
					|| 0 == secondSInfo.second.numberOfClusters_) {
				continue;
			}
			ret[sInfo.first][secondSInfo.first] = sInfo.second.fraction_
					- secondSInfo.second.fraction_;
		}
	}
	return ret;
}

void sampleCluster::resetInfos(){
	for(auto & info : sampInfos_){
		info.second.resetBasicInfo();
	}
	for(const auto & read : reads_){
		sampInfos_[read->sampName].update(read->seqBase_);
	}
	for(auto & info : sampInfos_){
		info.second.updateFraction();
	}
}

uint32_t sampleCluster::numberOfRuns() const {
	return sampleClusters_.size();
}

//
double sampleCluster::getCumulativeFrac() const {
  std::vector<double> fracs;
  for(const auto & samp : sampleClusters_){
  	double fracSum = 0;
  	for(const auto & clusPos : samp.second){
  		fracSum += reads_[clusPos]->seqBase_.frac_;
  	}
  	fracs.emplace_back(fracSum);
  }
  return vectorSum(fracs);
}

double sampleCluster::getAveragedFrac() const {
  std::vector<double> fracs;
  for(const auto & samp : sampleClusters_){
  	double fracSum = 0;
  	for(const auto & clusPos : samp.second){
  		fracSum += reads_[clusPos]->seqBase_.frac_;
  	}
  	fracs.emplace_back(fracSum);
  }
  return vectorMean(fracs);
}

double sampleCluster::getAveragedFrac(const VecStr & forClusters)const{
	return getAveragedFrac(getReadPositions(forClusters));
}

double sampleCluster::getCumulativeFrac(const VecStr & forClusters)const{
	return getCumulativeFrac(getReadPositions(forClusters));
}

double sampleCluster::getReadCnt(const VecStr & forClusters)const{
	return getReadCnt(getReadPositions(forClusters));
}

double sampleCluster::getAveragedFrac(const std::vector<uint32_t> & forClusters)const{
  std::vector<double> fracs;
  std::unordered_map<std::string, std::vector<double>> fracsBySamp;
  for(const auto & readPos : forClusters){
  	fracsBySamp[reads_[readPos]->sampName].emplace_back(reads_[readPos]->seqBase_.frac_);
  }
  for(const auto & samp : fracsBySamp){
  	fracs.emplace_back(vectorSum(samp.second));
  }
  return vectorMean(fracs);
}
double sampleCluster::getCumulativeFrac(const std::vector<uint32_t> & forClusters)const{
  std::vector<double> fracs;
  std::unordered_map<std::string, std::vector<double>> fracsBySamp;
  for(const auto & readPos : forClusters){
  	fracsBySamp[reads_[readPos]->sampName].emplace_back(reads_[readPos]->seqBase_.frac_);
  }
  for(const auto & samp : fracsBySamp){
  	fracs.emplace_back(vectorSum(samp.second));
  }
  return vectorSum(fracs);
}

double sampleCluster::getReadCnt(const std::vector<uint32_t> & forClusters)const{
	double sum = 0;
  for(const auto & readPos : forClusters){
  	sum += reads_[readPos]->seqBase_.cnt_;
  }
  return sum;
}

uint32_t sampleCluster::getSampNum(const std::vector<uint32_t> & forClusters)const{
	std::set<std::string> sampNames;
	for(const auto & readPos : forClusters){
		sampNames.emplace(reads_[readPos]->sampName);
	}
	return sampNames.size();
}

uint32_t sampleCluster::getSampNum(const VecStr & forClusters)const{
	return getSampNum(getReadPositions(forClusters));
}

double sampleCluster::getReadWeightedAveragedFrac()const{
  double sum = 0;
  //double readCnt = 0;
  for(const auto & i : sampInfos_){
  	//readCnt += i.second.readCnt_;
  	sum += i.second.runReadCnt_;
  }
  return seqBase_.cnt_/sum;
  //return readCnt/sum;
}

//so basically is for when exclusion has happened and so the total read count for sample/mid is now decreased
//this should only be called during sample collaping and never during population collapsing, that would turn things into read weighted averages
void sampleCluster::update(const std::map<std::string, sampInfo>& infos) {
  // clear and update it all
	for(const auto & info : infos){

		sampInfos_.at(info.first).updateRunReadCnt(info.second.readCnt_);
	}
  updateFractionInfo();
}
void sampleCluster::updateSampInfosFracs(){
	for(auto & sampInfo : sampInfos_){
		sampInfo.second.updateFraction();
	}
}
void sampleCluster::updateFractionInfo() {
  // clear the counts and update them
  seqBase_.frac_ = 0;
  seqBase_.cnt_ = 0;
  for (const auto& samp : sampInfos_) {
    seqBase_.cnt_ += samp.second.readCnt_;
  }

  //update the fractions for the subclusters as well
  for(auto & samp : sampleClusters_){
  	for(const auto & pos : samp.second){
  		reads_[pos]->seqBase_.frac_ = reads_[pos]->seqBase_.cnt_/sampInfos_[samp.first].runReadCnt_;
  	}
  }

  seqBase_.frac_ = getAveragedFrac();
  updateName();
}

void sampleCluster::updateName() {
  seqBase_.name_ = getStubName(false) + "_f" + estd::to_string(seqBase_.frac_);
}

void sampleCluster::setName(const std::string& newName) {
  seqBase_.name_ = newName + "_f" + estd::to_string(seqBase_.frac_);
}

std::string sampleCluster::getChimeraInfo(const std::string& delim) const {
	return vectorToString(getChimeraInfoVec(), delim);
}

VecStr sampleCluster::getChimeraInfoVec() const {
	uint32_t chiRepCnt = 0;
	uint32_t chiReadCnt = 0;
	uint32_t chiClusCnt = 0;
	for (const auto& info : sampInfos_) {
		if (info.second.chiReadCnt_ > 0) {
			++chiRepCnt;
			chiReadCnt += info.second.chiReadCnt_;
			chiClusCnt += info.second.chiNumberOfClusters_;
		}
	}
	return toVecStr(chiReadCnt, chiClusCnt, chiRepCnt);
}


VecStr sampleCluster::getClusterInfoVec() const{
	return toVecStr(seqBase_.frac_
  		, getReadWeightedAveragedFrac()
			, seqBase_.cnt_
			, numberOfRuns()
			, seqBase_.seq_
			, reads_.size()
			, getChimeraInfoVec()
			, vectorToString(readVec::getNames(reads_), ","));
}

VecStr sampleCluster::getClusterInfoHeaderVec() {
	return VecStr { "c_AveragedFrac", "c_ReadFrac", "c_ReadCnt", "c_RepCnt",
			"c_Consensus", "c_InputCnt", "c_ChiReadCnt", "c_ChiClusCnt",
			"c_ChiRepCnt", "c_InputNames" };
}

std::string sampleCluster::getClusterInfoHeader(const std::string & delim) {
	return bib::conToStr(getClusterInfoHeaderVec(), delim);
}

std::string sampleCluster::getClusterInfo(const std::string& delim ) const{
  std::stringstream outStream;
  outStream << seqBase_.frac_
  		<< delim << getReadWeightedAveragedFrac()
			<< delim << seqBase_.cnt_
			<< delim << numberOfRuns()
			<< delim << seqBase_.seq_
			<< delim << reads_.size()
			<< delim << getChimeraInfo(delim)
			<< delim << vectorToString(readVec::getNames(reads_), ",");
  return outStream.str();
}



std::string sampleCluster::getRepsInfoHeader(uint32_t maxRunCount, bool checkingExpected, const std::string & delim){
	std::stringstream ss;
  std::stringstream templateRunSring;
  templateRunSring << "RNUM_name"
  		<< delim << "RNUM_totalClusterCntExcluded"
			<< delim << "RNUM_totalCntExcluded"
			<< delim << "RNUM_totalFracExcluded"
			<< delim << "RNUM_clusterCntChiExcluded"
			<< delim << "RNUM_cntChiExcluded"
			<< delim << "RNUM_fracChiExcluded"
			<< delim << "RNUM_MapFrac"
			<< delim << "RNUM_ReadCnt"
			<< delim << "RNUM_ClusCnt"
			<< delim << "RNUM_totalReadCnt";
  for (uint32_t i = 1; i <= maxRunCount; ++i) {
  	if(i > 1){
  		ss << delim;
  	}
    ss << bib::replaceString(templateRunSring.str(), "NUM", std::to_string(i));
  }
  if(checkingExpected){
  	ss << delim << "c_bestExpected";
  }
	return ss.str();
}

std::string sampleCluster::getRepsInfo(
		const std::map<std::string, sampInfo> & input,
		const std::map<std::string, sampInfo> & excluded,
		const std::map<std::string, sampInfo> & collapsed, uint32_t maxRepCount,
		bool checkingExpected, const std::string& delim) const {
	uint32_t emptyRepAmount = 10;
	std::stringstream currentInfo;
	uint32_t infoCount = 0;
	for (const auto &info : collapsed) {
		if(infoCount > 0){
			currentInfo << delim;
		}
		++infoCount;
		auto search = sampInfos_.find(info.first);
		if (search->second.readCnt_ == 0) {
			currentInfo << repeatString(delim, emptyRepAmount);
		} else {
			currentInfo << info.first;
			if (excluded.find(info.first) == excluded.end()) {
				currentInfo << delim << 0 << delim << 0 << delim << 0 << delim << 0
						<< delim << 0 << delim << 0;
			} else {
				currentInfo << delim << excluded.at(info.first).numberOfClusters_;
				currentInfo << delim << excluded.at(info.first).readCnt_;
				currentInfo << delim
						<< excluded.at(info.first).readCnt_ / input.at(info.first).readCnt_;
				currentInfo << delim << excluded.at(info.first).chiNumberOfClusters_;
				currentInfo << delim << excluded.at(info.first).chiReadCnt_;
				currentInfo << delim
						<< excluded.at(info.first).chiReadCnt_
								/ input.at(info.first).readCnt_;
			}
			currentInfo << delim << sampInfos_.at(info.first).getReadInfo();
			currentInfo << delim << info.second.readCnt_;
		}
	}
	if (checkingExpected) {
		if(sampInfos_.size() < maxRepCount){
			for(uint32_t i = 0; i < maxRepCount - sampInfos_.size(); ++i){
				currentInfo << delim << repeatString(delim, emptyRepAmount);
			}
		}
		currentInfo << delim << expectsString;
	}
	return currentInfo.str();
}


VecStr sampleCluster::getPopHapInfoVec(double popReadCnt,
		uint32_t sampNum) const {
	return toVecStr(getCumulativeFrac() / sampNum, getCumulativeFrac(),
			seqBase_.frac_, seqBase_.cnt_ / popReadCnt, sampleClusters_.size(),
			sampleClusters_.size() / static_cast<double>(sampNum), seqBase_.cnt_,
			reads_.size(), vectorToString(readVec::getNames(reads_), ","),
			seqBase_.seq_);
}

VecStr sampleCluster::getPopHapInfoHeaderVec(){
	return VecStr { "h_PopFrac", "h_SumOfAllFracs", "h_AvgFracFoundAt",
			"h_ReadFrac", "h_SampCnt", "h_SampFrac", "h_ReadCnt", "h_ClusterCnt",
			"h_clusterNames", "h_Consesus"};
}


std::string sampleCluster::getFullPopInfoHeader(const std::string& coiHeader, const std::string & delim ){
	std::stringstream ss;
	ss 	<< "h_popUID"
			<< delim << "p_TotalInputReadCnt"
			<< delim << "p_TotalInputClusterCnt"
			<< delim << "p_TotalPopulationSampCnt"
			<< delim << "p_TotalHaplotypes"
			<< delim << coiHeader
			<< delim << "h_PopFrac"
			<< delim << "h_SumOfAllFracs"
			<< delim << "h_AvgFracFoundAt" //should change this to median
			<< delim << "h_ReadFrac"
			<< delim << "h_SampCnt"
			<< delim << "h_SampFrac"
			<< delim << "h_ReadCnt"
			<< delim << "h_ClusterCnt"
			<< delim << "h_clusterNames"
			<< delim << "h_Consesus"
			<< delim << "h_Protein";
	return ss.str();
}

std::string sampleCluster::getFullPopInfo(double popReadCnt,
		uint32_t popClusNum, uint32_t sampNum, uint32_t numOfHaps,
		const std::string& coiInfo,
		const std::string& delim) const {
  std::stringstream outStream;
  outStream << getStubName(false)
  		<< delim << popReadCnt
			<< delim << popClusNum
			<< delim << sampNum
			<< delim << numOfHaps
			<< delim << coiInfo
      << delim << getCumulativeFrac() / sampNum
			<< delim << getCumulativeFrac()
			<< delim << seqBase_.frac_
			<< delim << seqBase_.cnt_ / popReadCnt
			<< delim << sampleClusters_.size()
			<< delim << sampleClusters_.size() / static_cast<double>(sampNum)
			<< delim << seqBase_.cnt_
			<< delim << reads_.size()
			<< delim << vectorToString(readVec::getNames(reads_), ",")
			<< delim << seqBase_.seq_
			<< delim << translateRet(false, false).seqBase_.seq_;
  return outStream.str();
}



std::string sampleCluster::getPopInfoHeader(const std::string & delim){

	return bib::conToStr(getPopInfoHeaderVec(), delim);
}

VecStr sampleCluster::getPopInfoHeaderVec() {
	return VecStr { "h_popUID", "p_TotalPopulationSampCnt", "h_PopFrac",
			"h_SumOfAllFracs", "h_AvgFracFoundAt", "h_ReadFrac", "h_SampCnt",
			"h_SampFrac", "h_ReadCnt", "h_ClusterCnt", "h_clusterNames" };
}

VecStr sampleCluster::getPopInfoVec(double popReadCnt, uint32_t popClusNum,
                               uint32_t sampNum) const{
	return toVecStr(getStubName(false)
			, sampNum
      , getCumulativeFrac() / sampNum
			, getCumulativeFrac()
			, seqBase_.frac_
			, seqBase_.cnt_ / popReadCnt
			, sampleClusters_.size()
			, sampleClusters_.size() / static_cast<double>(sampNum)
			, seqBase_.cnt_
			, reads_.size()
			, vectorToString(readVec::getNames(reads_), ","));
}

std::string sampleCluster::getPopInfo(double popReadCnt, uint32_t popClusNum,
		uint32_t sampNum, const std::string& delim) const {
	return bib::conToStr(getPopInfoVec(popReadCnt, popClusNum, sampNum), delim);
}

std::vector<uint32_t> sampleCluster::getReadPositions(const VecStr & forClusters)const{
	std::vector<uint32_t> ret;
	for(const auto & readPos : iter::range(reads_.size())){
		if(bib::in(reads_[readPos]->seqBase_.name_, forClusters)){
			ret.emplace_back(readPos);
		}
	}
	return ret;
}



}  // namespace bib
