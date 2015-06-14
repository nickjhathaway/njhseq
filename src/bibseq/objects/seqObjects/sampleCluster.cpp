//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "sampleCluster.hpp"
#include <bibcpp/stdAddition.h>
namespace bibseq {
sampleCluster::sampleCluster(const readObject& initializerRead,
		const std::map<std::string, sampInfo>& infos)
		 {
	seqBase_ = initializerRead.seqBase_;
	firstReadName = initializerRead.seqBase_.name_;
	firstReadCount = initializerRead.seqBase_.cnt_;

	sampName = getOwnSampName();

	for (const auto & i : infos) {
		sampInfos_[i.second.runName_] = sampInfo(i.second.runName_,
				i.second.runReadCnt_);
	}

	reads_.emplace_back(*this);
	updateInfoWithRead(*this, 0);
	remove = false;
}


void sampleCluster::addRead(const cluster& cr) {
  // add the cluster's reads and update the info
  for (const auto& read : cr.reads_) {
    updateInfoWithRead(read, reads_.size());
    reads_.emplace_back(read);
  }
  // update the fraction and totalCounts
  seqBase_.cnt_ += cr.seqBase_.cnt_;
  seqBase_.frac_ = getAveragedFrac(); //calculate average as the mean average between all
  updateName();
  //consider putting in parameters to determine if
  needToCalculateConsensus = true;
}

void sampleCluster::updateInfoWithRead(const readObject& read, uint32_t pos) {
  // update infos with the read
	sampleClusters_[read.sampName].push_back(pos);
  sampInfos_[read.sampName].update(read);
}

void sampleCluster::resetInfos(){
	for(auto & info : sampInfos_){
		info.second.resetBasicInfo();
	}
	for(const auto & read : reads_){
		sampInfos_[read.sampName].update(read);
	}
	for(auto & info : sampInfos_){
		info.second.updateFraction();
	}
}

uint32_t sampleCluster::numberOfRuns()const {
	return sampleClusters_.size();
}

//
double sampleCluster::getCumulativeFrac()const{
  std::vector<double> fracs;
  for(const auto & samp : sampleClusters_){
  	double fracSum = 0;
  	for(const auto & clusPos : samp.second){
  		fracSum += reads_[clusPos].seqBase_.frac_;
  	}
  	fracs.emplace_back(fracSum);
  }
  return vectorSum(fracs);
}

double sampleCluster::getAveragedFrac()const{
  std::vector<double> fracs;
  for(const auto & samp : sampleClusters_){
  	double fracSum = 0;
  	for(const auto & clusPos : samp.second){
  		fracSum += reads_[clusPos].seqBase_.frac_;
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
  	fracsBySamp[reads_[readPos].sampName].emplace_back(reads_[readPos].seqBase_.frac_);
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
  	fracsBySamp[reads_[readPos].sampName].emplace_back(reads_[readPos].seqBase_.frac_);
  }
  for(const auto & samp : fracsBySamp){
  	fracs.emplace_back(vectorSum(samp.second));
  }
  return vectorSum(fracs);
}

double sampleCluster::getReadCnt(const std::vector<uint32_t> & forClusters)const{
	double sum = 0;
  for(const auto & readPos : forClusters){
  	sum += reads_[readPos].seqBase_.cnt_;
  }
  return sum;
}

uint32_t sampleCluster::getSampNum(const std::vector<uint32_t> & forClusters)const{
	std::set<std::string> sampNames;
	for(const auto & readPos : forClusters){
		sampNames.emplace(reads_[readPos].sampName);
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
  		reads_[pos].seqBase_.frac_ = reads_[pos].seqBase_.cnt_/sampInfos_[samp.first].runReadCnt_;
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
  int chiRepCnt = 0;
  int chiReadCnt = 0;
  int chiClusCnt = 0;
  for (const auto& info : sampInfos_) {
    if (info.second.chiReadCnt_ > 0) {
      ++chiRepCnt;
      chiReadCnt += info.second.chiReadCnt_;
      chiClusCnt += info.second.chiNumberOfClusters_;
    }
  }
  return vectorToString(toVecStr(chiReadCnt,chiClusCnt ,chiRepCnt), delim);
}


std::string sampleCluster::getSampleInfoHeader(const std::string & delim ){
	std::stringstream ss;
	ss 	<< "c_AveragedFrac"
			<< delim << "c_ReadFrac"
			<< delim << "c_ReadCnt"
			<< delim << "c_RepCnt"
			<< delim << "c_Consensus"
			<< delim << "c_InputCnt"
			<< delim << "c_ChiReadCnt"
			<< delim << "c_ChiClusCnt"
			<< delim << "c_ChiRepCnt"
			<< delim << "c_InputNames";
	return ss.str();
}

std::string sampleCluster::getSampleInfo(const std::string& delim ) const{
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



std::string sampleCluster::getRepsInfoHeader(uint32_t maxRunCount,bool checkingExpected, const std::string & delim){
	std::stringstream ss;
  std::stringstream templateRunSring;
  templateRunSring << "RNUM.name"
  		<< delim << "RNUM.totalClusterCntExcluded"
			<< delim << "RNUM.totalCntExcluded"
			<< delim << "RNUM.totalFracExcluded"
			<< delim << "RNUM.clusterCntChiExcluded"
			<< delim << "RNUM.cntChiExcluded"
			<< delim << "RNUM.fracChiExcluded"
			<< delim << "RNUM.MapFrac"
			<< delim << "RNUM.ReadCnt"
			<< delim << "RNUM.ClusCnt"
			<< delim << "RNUM.totalReadCnt";
  for (uint32_t i = 1; i <= maxRunCount; ++i) {
  	if(i > 1){
  		ss << delim;
  	}
    ss << replaceString(templateRunSring.str(), "NUM", std::to_string(i));
  }
  if(checkingExpected){
  	ss << delim << "bestExpected";
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
		if(collapsed.size() < maxRepCount){
			currentInfo << delim << repeatString(delim, emptyRepAmount);
		}
		currentInfo << delim << expectsString;
	}
	return currentInfo.str();
}



std::string sampleCluster::getFullPopInfoHeader(const std::string& moiHeader, const std::string & delim ){
	std::stringstream ss;
	ss 	<< "h_popUID"
			<< delim << "p_TotalInputReadCnt"
			<< delim << "p_TotalInputClusterCnt"
			<< delim << "p_TotalPopulationSampCnt"
			<< delim << "p_TotalHaplotypes"
			<< delim << moiHeader
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
		const std::string& moiInfo,
		const std::string& delim) const {
  std::stringstream outStream;
  outStream << getStubName(false)
  		<< delim << popReadCnt
			<< delim << popClusNum
			<< delim << sampNum
			<< delim << numOfHaps
			<< delim << moiInfo
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
			<< delim << getProteinFromcDNA(false);
  return outStream.str();
}

std::string sampleCluster::getFullPopInfo(double popReadCnt,
		uint32_t popClusNum, uint32_t sampNum, uint32_t numOfUniHaps,
		const std::string& moiInfo, const VecStr & forClusters,
		const std::string& delim) const {
  auto subReadPositions = getReadPositions(forClusters);
  auto readCnt = getReadCnt(subReadPositions);
  auto cumulativeFrac = getCumulativeFrac(subReadPositions);
  auto avgFrac = getAveragedFrac(subReadPositions);
  auto currentSampNum = getSampNum(subReadPositions);
  VecStr subReadNames;
  for(const auto & pos : subReadPositions){
  	subReadNames.emplace_back(reads_[pos].seqBase_.name_);
  }
  std::stringstream outStream;
  outStream << getStubName(false)
  		<< delim << popReadCnt
			<< delim << popClusNum
			<< delim << sampNum
			<< delim << numOfUniHaps
			<< delim << moiInfo
			<< delim << cumulativeFrac / sampNum
			<< delim << cumulativeFrac
			<< delim << avgFrac
			<< delim << readCnt / popReadCnt
			<< delim << currentSampNum
			<< delim << currentSampNum / static_cast<double>(sampNum)
			<< delim << readCnt
			<< delim << subReadPositions.size()
			<< delim << vectorToString(subReadNames, ",")
			<< delim << seqBase_.seq_
			<< delim << getProteinFromcDNA(false);
  return outStream.str();
}


std::string sampleCluster::getPopInfoHeader(const std::string & delim){
	std::stringstream ss;
	ss 	<< "h_popUID"
			<< delim << "p_TotalPopulationSampCnt"
			<< delim << "h_PopFrac"
			<< delim << "h_SumOfAllFracs"
			<< delim << "h_AvgFracFoundAt" //should change this to median
			<< delim << "h_ReadFrac"
			<< delim << "h_SampCnt"
			<< delim << "h_SampFrac"
			<< delim << "h_ReadCnt"
			<< delim << "h_ClusterCnt"
			<< delim << "h_clusterNames";
	return ss.str();
}

std::string sampleCluster::getPopInfo(double popReadCnt,
                                              uint32_t popClusNum,
                                              uint32_t sampNum,
                                              const std::string& delim) const {
  std::stringstream outStream;
  outStream << getStubName(false)
			<< delim << sampNum
      << delim << getCumulativeFrac() / sampNum
			<< delim << getCumulativeFrac()
			<< delim << seqBase_.frac_
			<< delim << seqBase_.cnt_ / popReadCnt
			<< delim << sampleClusters_.size()
			<< delim << sampleClusters_.size() / static_cast<double>(sampNum)
			<< delim << seqBase_.cnt_
			<< delim << reads_.size()
			<< delim << vectorToString(readVec::getNames(reads_), ",");
  return outStream.str();
}

std::vector<uint32_t> sampleCluster::getReadPositions(const VecStr & forClusters)const{
	std::vector<uint32_t> ret;
	for(const auto & readPos : iter::range(reads_.size())){
		if(bib::in(reads_[readPos].seqBase_.name_, forClusters)){
			ret.emplace_back(readPos);
		}
	}
	return ret;
}

std::string sampleCluster::getPopInfo(double popReadCnt,
                                              uint32_t popClusNum,
                                              uint32_t sampNum, const VecStr & forClusters,
                                              const std::string& delim) const {
  std::stringstream outStream;

  auto subReadPositions = getReadPositions(forClusters);
  auto readCnt = getReadCnt(subReadPositions);
  auto cumulativeFrac = getCumulativeFrac(subReadPositions);
  auto avgFrac = getAveragedFrac(subReadPositions);
  auto currentSampNum = getSampNum(subReadPositions);
  VecStr subReadNames;
  for(const auto & pos : subReadPositions){
  	subReadNames.emplace_back(reads_[pos].seqBase_.name_);
  }
  outStream << getStubName(false)
			<< delim << sampNum
      << delim << cumulativeFrac / sampNum
			<< delim << cumulativeFrac
			<< delim << avgFrac
			<< delim << readCnt / popReadCnt
			<< delim << currentSampNum
			<< delim << currentSampNum / static_cast<double>(sampNum)
			<< delim << readCnt
			<< delim << subReadPositions.size()
			<< delim << vectorToString(subReadNames, ",");
  return outStream.str();
}

}  // namespace bib
