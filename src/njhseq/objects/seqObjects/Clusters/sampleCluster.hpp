#pragma once
//
//  sampleCluster.hpp
//
//  Created by Nicholas Hathaway on 9/10/13.
//
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
#include "njhseq/objects/seqObjects/Clusters/cluster.hpp"
#include "njhseq/objects/helperObjects/sampInfo.hpp"

namespace njhseq {



class sampleCluster : public cluster {

 public:

	//template <typename T>
	sampleCluster(const seqInfo& initializerRead,
			const std::map<std::string, sampInfo>& infos);

	/** This is for reconstructing a sampleCluster from after writing it out
	 *
	 * @param initializerRead The seq base
	 * @param seqs the seqs clustered into this cluster
	 * @param infos the infos for the totals for the different replicates for this sample
	 */
	template<typename T>
	sampleCluster(const seqInfo& initializerRead,
			const std::vector<T> & seqs,
			const std::map<std::string, sampInfo>& infos){
		seqBase_ = initializerRead;

		firstReadName_ = getSeqBase(seqs.front()).name_;
		firstReadCount_ = getSeqBase(seqs.front()).cnt_;
		seqBase_.cnt_ = 0;
		sampName = getOwnSampName();

		for (const auto & i : infos) {
			sampInfos_[i.second.runName_] = sampInfo(i.second.runName_,
					i.second.runReadCnt_);
		}
		for(const auto & seq : seqs){
			sampleClusters_[getSeqBase(seq).getOwnSampName()].push_back(reads_.size());
			sampInfos_[getSeqBase(seq).getOwnSampName()].update(getSeqBase(seq));
			reads_.emplace_back(std::make_shared<readObject>(getSeqBase(seq)));
		  // update the fraction and totalCounts
		  seqBase_.cnt_ += getSeqBase(seq).cnt_;
		  seqBase_.frac_ = getAveragedFrac(); //calculate average as the mean fraction between all samples
		}

		setLetterCount();
		remove = false;
		needToCalculateConsensus_ = false;
	}

	template<typename T>
	sampleCluster(const T& initializerRead,
			const std::map<std::string, sampInfo>& infos) :
			sampleCluster(getSeqBase(initializerRead), infos) {

	}

  virtual void updateName();
  virtual void setName(const std::string& newName);



  void addRead(const sampleCluster& cr);
  void updateInfoWithRead(const readObject& read, uint32_t pos);

  void update(const std::map<std::string, sampInfo>& infos);

  void updateFractionInfo();

  void updateSampInfosFracs();
  // info
	uint32_t numberOfRuns() const;

	double getAveragedFrac() const;
	double getCumulativeFrac() const;

  double getAveragedFrac(const VecStr & forClusters)const;
  double getCumulativeFrac(const VecStr & forClusters)const;
  double getReadCnt(const VecStr & forClusters)const;
  uint32_t getSampNum(const VecStr & forClusters)const;
 private:
  double getAveragedFrac(const std::vector<uint32_t> & forClusters)const;
  double getCumulativeFrac(const std::vector<uint32_t> & forClusters)const;
  double getReadCnt(const std::vector<uint32_t> & forClusters)const;
  uint32_t getSampNum(const std::vector<uint32_t> & forClusters)const;
 public:

  std::vector<uint32_t> getReadPositions(const VecStr & forClusters)const;

  double getReadWeightedAveragedFrac()const;

  std::string getChimeraInfo(const std::string& delim = "\t") const;
  VecStr getChimeraInfoVec() const;

  std::string getClusterInfo(const std::string& delim = "\t") const;
  static std::string getClusterInfoHeader(const std::string & delim = "\t");

  VecStr getClusterInfoVec() const;
  static VecStr getClusterInfoHeaderVec();

	std::string getRepsInfo(const std::map<std::string, sampInfo> & input,
			const std::map<std::string, sampInfo> & excluded,
			const std::map<std::string, sampInfo> & collapsed,uint32_t maxRepCount, bool checkingExpected,
			const std::string& delim = "\t") const;

  static std::string getRepsInfoHeader(uint32_t maxRunCount,bool checkingExpected, const std::string & delim = "\t");

  std::string getPopInfo(double popReadCnt, uint32_t popClusNum,
                                 uint32_t sampNum,
                                 const std::string& delim = "\t") const;

  VecStr getPopInfoVec(double popReadCnt, uint32_t popClusNum,
                                 uint32_t sampNum) const;

  static std::string getPopInfoHeader(const std::string & delim = "\t");
  static VecStr getPopInfoHeaderVec();

  std::string getFullPopInfo(double popReadCnt, uint32_t popClusNum,
                                 uint32_t sampNum,uint32_t numOfHaps, const std::string& coiHeader,
                                 const std::string& delim = "\t") const;

  static std::string getFullPopInfoHeader(const std::string& coiHeader, const std::string & delim = "\t");


  VecStr getPopHapInfoVec(double popReadCnt, uint32_t sampNum) const;
	static VecStr getPopHapInfoHeaderVec();



  void resetInfos();

  std::unordered_map<std::string, std::unordered_map<std::string, double>> getRepDisagreement() const;

protected:
  std::map<std::string, std::vector<uint32_t>> sampleClusters_;
  std::map<std::string, sampInfo> sampInfos_;
public:
  const std::map<std::string, sampInfo> & sampInfos()const{return sampInfos_;}
  const std::map<std::string, std::vector<uint32_t>> & sampleClusters()const{return sampleClusters_;}

  void setSampInfos(const std::map<std::string, sampInfo> & sampInfos){
  	sampInfos_ = sampInfos;
  }

  void setSampInfosTotals(const std::map<std::string, sampInfo> & sampInfos){
  	sampInfos_ = sampInfos;
  	for(auto & info : sampInfos_){
  		info.second.resetBasicInfo();
  	}
		for (const auto & seq : reads_) {
			sampInfos_[getSeqBase(seq).getOwnSampName()].update(getSeqBase(seq));
		}
  }

  template <typename T>
  static void updateAllClusters(std::vector<T>& clusters,
                                const std::map<std::string, sampInfo>& infos) {
    for (auto& clus : clusters) {
      clus.update(infos);
    }
  }
  using size_type = baseReadObject::size_type;
};


template<>
inline sampleCluster::size_type len(const sampleCluster & read){
	return read.seqBase_.seq_.size();
}

}  // namespace njhseq


