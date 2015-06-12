#pragma once
//
//  sampleCluster.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/10/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/objects/seqObjects/cluster.hpp"
#include "bibseq/objects/helperObjects/sampInfo.hpp"

namespace bibseq {



class sampleCluster : public cluster {

 public:
  //template <typename T>
  sampleCluster(const readObject& initializerRead, const std::map<std::string, sampInfo>& infos) ;

  //sampleCluster(){}

  virtual void updateName();
  virtual void setName(const std::string& newName);



  void addRead(const cluster& cr);
  void updateInfoWithRead(const readObject& read, uint32_t pos);

  void update(const std::map<std::string, sampInfo>& infos);

  void updateFractionInfo();
  // info
  uint32_t numberOfRuns()const;

  double getAveragedFrac()const;
  double getCumulativeFrac()const;

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

  std::string getSampleInfo(const std::string& delim = "\t") const;
  static std::string getSampleInfoHeader(const std::string & delim = "\t");

	std::string getRepsInfo(const std::map<std::string, sampInfo> & input,
			const std::map<std::string, sampInfo> & excluded,
			const std::map<std::string, sampInfo> & collapsed,uint32_t maxRepCount, bool checkingExpected,
			const std::string& delim = "\t") const;

  static std::string getRepsInfoHeader(uint32_t maxRunCount,bool checkingExpected, const std::string & delim = "\t");


  std::string getPopInfo(double popReadCnt, uint32_t popClusNum,
                                 uint32_t sampNum,const VecStr & forClusters,
                                 const std::string& delim = "\t") const;

  std::string getPopInfo(double popReadCnt, uint32_t popClusNum,
                                 uint32_t sampNum,
                                 const std::string& delim = "\t") const;

  static std::string getPopInfoHeader(const std::string & delim = "\t");

  std::string getFullPopInfo(double popReadCnt, uint32_t popClusNum,
                                 uint32_t sampNum,uint32_t numOfHaps, const std::string& moiHeader,
                                 const std::string& delim = "\t") const;

	std::string getFullPopInfo(double popReadCnt, uint32_t popClusNum,
			uint32_t sampNum, uint32_t numOfUniHaps, const std::string& moiInfo,
			const VecStr & forClusters, const std::string& delim = "\t") const;

  static std::string getFullPopInfoHeader(const std::string& moiHeader, const std::string & delim = "\t");

  void resetInfos();

protected:
  std::map<std::string, std::vector<uint32_t>> sampleClusters_;
  std::map<std::string, sampInfo> sampInfos_;
public:
  const std::map<std::string, sampInfo> & sampInfos()const{return sampInfos_;}
  void setSampInfos(const std::map<std::string, sampInfo> & sampInfos){
  	sampInfos_ = sampInfos;
  }
  template <typename T>
  static void updateAllClusters(std::vector<T>& clusters,
                                const std::map<std::string, sampInfo>& infos) {
    for (auto& clus : clusters) {
      clus.update(infos);
    }
  }
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "sampleCluster.cpp"
#endif
