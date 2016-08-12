#pragma once
//
//  cluster.h
//  ampliconCluster
//
//  Created by Nicholas Hathaway on 7/17/12.
//  Copyright (c) 2012 University of Massachusetts Medical School. All rights
// reserved.
//
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
#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/seqObjects/Clusters/baseCluster.hpp"
#include "bibseq/readVectorManipulation/readVectorOperations.h"
#include "bibseq/alignment/aligner.h"

#include "bibseq/IO/SeqIO/SeqOutput.hpp"
#include "bibseq/simulation/errorProfile.hpp"

namespace bibseq {
// class simulation::errorProfile;

class cluster : public baseCluster {
 public:
  // constructors
  cluster();

  cluster(const seqInfo &firstRead);


  void addRead(const cluster &read);

  void removeRead(const std::string &stubName);
  void removeReads(const std::vector<readObject> &vec);
  bool rejected_;

  //std::mutex mtx_;


  std::vector<cluster> chimeras;
  std::vector<cluster> endChimeras;
  std::vector<cluster> allInputClusters;


  void outputInfoComp(const std::string &workingDir) const;

  // find the seed with the largest amount of identical reads
  int getBiggestReadSize();

  double getAverageSizeDifference();
  size_t getLargestSizeDifference();



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
    SeqIOOptions options = SeqIOOptions::genFastqOut(directory + "reads");
    SeqOutput writer(options);

    writer.openWrite(inReads);
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
      double medianLength = vectorMedianRef(readLengths);
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
                        << vectorMedianCopy(sv.second.second) << std::endl;
    }
    std::ofstream sizeVariationGraphFile;
    openTextFile(sizeVariationGraphFile, directory + "sizeVariationGraph.tab.txt", ".txt", false, false);
    sizeVariationGraphFile << "sizeRange\tfrequency\tmeanFraqBase>" << qualCheck
                           << "\tmedianFraqBase>" << qualCheck << std::endl;
    for (const auto &sv : sizeVariation) {
      sizeVariationGraphFile << sv.first * 10 + 5 << "\t" << sv.second.first
                             << "\t" << vectorMean(sv.second.second) << "\t"
                             << vectorMedianCopy(sv.second.second) << std::endl;
    }
  }
	using size_type = baseReadObject::size_type;
};

template<>
inline cluster::size_type len(const cluster & read) {
	return read.seqBase_.seq_.size();
}


template<typename T>
cluster getConsensus(const std::vector<T> & reads, aligner & alignerObj,
		const std::string & name){
	if(reads.empty()){
		return cluster();
	}
	cluster mainCluster(reads.begin()->seqBase_);
	if(reads.size() > 1){
		std::vector<cluster> inClusters;
		for(const auto & readPos : iter::range<uint64_t>(1,reads.size())){
			inClusters.emplace_back(cluster(reads[readPos].seqBase_));
		}
		for(const auto & clus : inClusters){
			mainCluster.addRead(clus);
		}
	}
	mainCluster.calculateConsensus(alignerObj, true);
	mainCluster.setName(name);
	return mainCluster;
}
}  // namespace bibseq


