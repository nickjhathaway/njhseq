#pragma once
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
//
//  cluster.h
//  ampliconCluster
//
//  Created by Nicholas Hathaway on 7/17/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//

#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/seqObjects/baseCluster.hpp"
#include "bibseq/objects/seqObjects/identicalCluster.hpp"
#include "bibseq/readVectorManipulation/readVectorOperations.h"
#include "bibseq/alignment/aligner.hpp"

#include "bibseq/IO/readObjectIO.hpp"
#include "bibseq/simulation/errorProfile.hpp"

namespace bibseq {
// class simulation::errorProfile;

class cluster : public baseCluster {
 public:
  // constructors
  cluster() : baseCluster() { rejected_ = false; }

  cluster(const readObject &firstRead) : baseCluster(firstRead) {
    rejected_ = false;
    setLetterCount();
    // std::cout <<"calling the cluster first read constructor " << std::endl;
  }

  /*
  cluster(const cluster & other): baseCluster(other), rejected_(other.rejected_), mtx_(),
  		frontChiPos(other.frontChiPos), endChiPos(other.endChiPos),
  		endChimeras(other.endChimeras), allInputClusters(other.allInputClusters),
  		amountAverageMap_(other.amountAverageMap_), pValueMap_(other.pValueMap_){
  	setLetterCount();
  }*/
/*
  cluster& operator= (const cluster& other){
  	baseCluster::operator=(other);
  	rejected_ = other.rejected_;
  	frontChiPos = other.frontChiPos;
  	endChiPos = other.endChiPos;
  	endChimeras = other.endChimeras;
  	allInputClusters = other.allInputClusters;
  	amountAverageMap_ = other.amountAverageMap_;
  	pValueMap_ = other.pValueMap_;

  	return *this;
  }*/

  void addRead(const baseCluster &read);

  void removeRead(const std::string &stubName);
  void removeReads(const std::vector<readObject> &vec);
  bool rejected_;

  //std::mutex mtx_;
  // chimeras vectors
  // key is the position at which they could be chimeric with the
  // value which is the positions in the vector they are both in...
  // not ideal way to do this i guess
  std::multimap<uint32_t, uint32_t> frontChiPos;
  std::multimap<uint32_t, uint32_t> endChiPos;

  std::vector<cluster> chimeras;
  std::vector<cluster> endChimeras;
  std::vector<cluster> allInputClusters;
  // simulation stats
  std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>
      amountAverageMap_;
  std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> pValueMap_;
  // void gatherAllClusterInfo(std::ofstream & fastaFile, std::ofstream &
  // qualFile)const;
  // void gatherAllClusters(std::vector<cluster> & allClusters)const;

  void outputInfoComp(const std::string &workingDir) const;
  // void gatherInfoAboutIdenticalReadCompostion(std::map<int,int>&
  // clusterCounter)const;

  // find the seed with the largest amount of identical reads
  int getBiggestReadSize();
  //
  double getAverageSizeDifference();
  int getLargestSizeDifference();


  void writeOutInputClusters(const std::string &workingDir) const;


  // fdr calc
  double getFDRValue(uint32_t clusSize, uint32_t mismatches,
                     uint32_t observedAmount);
  // pvalue calc
  double getPValue(const seqInfo &seqBase, uint32_t numOfMismatches);
  double getPValue(const seqInfo &seqBase, aligner &alignerObj, bool local);

  //
  double getAverageFrequencyOfClusterVector() const;

  std::vector<readObject> getAllBeginingClusters(
      const std::vector<identicalCluster> &initialClusters) const;
  std::vector<std::vector<uint32_t>> alignOrigQuals(
      const std::vector<readObject> &reads, aligner &alingerObj,
      bool local) const;
  std::vector<uint32_t> getQualToConsensus(
      const readObject &read, aligner &alignerObj, bool local,
      std::map<std::string, std::vector<uint32_t>> &positions) const;
  void simOnQual(const std::vector<identicalCluster> &initialClusters,
                 aligner &alingerObj, uint32_t runTimes,
                 std::unordered_map<double, double> bestLikelihood,
                 simulation::mismatchProfile &eProfile, randomGenerator &gen);

	std::unordered_map<std::string, uint32_t> mutateConensus(
			const std::vector<std::vector<uint32_t>> &currentAlignQuals,
			std::unordered_map<double, double> bestLikelihood,
			simulation::mismatchProfile &eProfile, randomGenerator &gen);

	std::unordered_map<std::string, uint32_t> mutateConensus(double errorRate,
			simulation::mismatchProfile &eProfile, randomGenerator &gen);

	std::unordered_map<uint32_t,
			std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>>>simulate(uint32_t runTimes,
			const std::vector<std::vector<uint32_t>> &currentAlignQuals,
			std::unordered_map<double, double> bestLikelihood,
			simulation::mismatchProfile &eProfile, randomGenerator &gen);
	void postProcessSimulation(
			uint32_t runTimes,
			const std::unordered_map<
			uint32_t, std::unordered_map<
			uint32_t, std::unordered_map<uint32_t, uint32_t>>> &
			allRunCounts);
	std::unordered_map<uint32_t,
	std::unordered_map<uint32_t, std::vector<cluster>>>
	createSimReadClusters(aligner &alignerObj);
	void removeClustersOnFDR(aligner &alingerObj, double fdrCutOff,
			std::vector<cluster> &rejectedClusters);
  //////get infos
  template <typename T>
  static void getInfo(std::vector<T> &inReads, const std::string &directory,
                      const std::string &filename, int qualCheck) {
    std::map<uint32_t, std::vector<uint32_t>, std::less<uint32_t>> readsBySize;

    double totalAmountOfReads = 0;
    //std::vector<T> clusteredReads;
    for (const auto & readPos : iter::range(len(inReads))) {
      //clusteredReads.push_back(r);
      totalAmountOfReads += inReads[readPos].seqBase_.cnt_;
      readsBySize[inReads[readPos].seqBase_.cnt_].emplace_back(readPos);
    }
    readObjectIOOptions options;
    options.outFilename_ = directory + "reads";
    options.outFormat_ = "fastq";
    readObjectIO::writeVector(inReads, options);
    std::ofstream info;
    openTextFile(info, directory + filename, ".tab.txt", false, false);
    info << "ClusterSize\tnumberOfClusters\ttotalReads(%)"
            "\tShortestLength\tLongestLength\tMeanLength\tMedianLength"
         << std::endl;
    for (const auto & readsSizes : readsBySize) {
      info << readsSizes.first << "\t"
           << readsSizes.second.size() << "\t";
      double totalReads = 0;
      std::vector<uint32_t> readLengths;
      for (const auto & readsPosition : readsSizes.second) {
        totalReads += inReads[readsPosition].seqBase_.cnt_;
        readLengths.emplace_back(len(inReads[readsPosition].seqBase_.seq_));
      }
      auto shortestLength = vectorMinimum(readLengths);
      auto longestLength = vectorMaximum(readLengths);
      double medianLength = vectorMedian(readLengths);
      double meanLength = vectorMean(readLengths);
      info << getPercentageString(totalReads, totalAmountOfReads)
           << "\t";
      info << shortestLength << "\t" << longestLength << "\t" << meanLength
           << "\t" << medianLength << std::endl;
    }
    readVec::allSetQualCheck(inReads, qualCheck);
    std::map<size_t, std::pair<int, std::vector<double>>> sizeVariation;
    //get info about the reads rounded to the closet 10th
    uint64_t maxSize = 0;
    readVec::getMaxLength(inReads, maxSize);
   // size_t maxSizeRounded = maxSize / 10;
    /*for (size_t i = 0; i <= maxSizeRounded; ++i) {
      sizeVariation[i] = {0, {}};
    }*/
    for (const auto &r : inReads) {
      size_t sizeRounded = r.seqBase_.seq_.length() / 10;
      sizeVariation[sizeRounded].first += r.seqBase_.cnt_;
      sizeVariation[sizeRounded].second.emplace_back(r.fractionAboveQualCheck_);
    }
    std::ofstream sizeVariationFile;
    openTextFile(sizeVariationFile, directory + "sizeVariation.tab.txt", ".txt", false, false);
    sizeVariationFile << "sizeRange\tfrequency\tmeanFraqBase>" << qualCheck
                      << "\tmedianFraqBase>" << qualCheck << std::endl;
    for (const auto &sv : sizeVariation) {
      sizeVariationFile << sv.first * 10 << "-" << sv.first * 10 + 9 << "\t"
                        << sv.second.first << "\t"
                        << vectorMean(sv.second.second) << "\t"
                        << vectorMedian(sv.second.second) << std::endl;
    }
    std::ofstream sizeVariationGraphFile;
    openTextFile(sizeVariationGraphFile, directory + "sizeVariationGraph.tab.txt", ".txt", false, false);
    sizeVariationGraphFile << "sizeRange\tfrequency\tmeanFraqBase>" << qualCheck
                           << "\tmedianFraqBase>" << qualCheck << std::endl;
    for (const auto &sv : sizeVariation) {
      sizeVariationGraphFile << sv.first * 10 + 5 << "\t" << sv.second.first
                             << "\t" << vectorMean(sv.second.second) << "\t"
                             << vectorMedian(sv.second.second) << std::endl;
    }
  }
};


template<typename T>
cluster getConsensus(const std::vector<T> & reads, aligner & alignerObj,
		const std::string & name){
	readObject firstRead = readObject(reads.begin()->seqBase_);
	//reads.erase(reads.begin());
	cluster mainCluster(firstRead);
	std::vector<cluster> inClusters;
	for(const auto & readPos : iter::range<uint64_t>(1,reads.size())){
		inClusters.emplace_back(cluster(readObject(reads[readPos].seqBase_)));
	}
	for(const auto & clus : inClusters){
		mainCluster.addRead(clus);
	}
	mainCluster.calculateConsensus(alignerObj, true);

	mainCluster.setName(name);
	return mainCluster;
}
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "cluster.cpp"
#endif
