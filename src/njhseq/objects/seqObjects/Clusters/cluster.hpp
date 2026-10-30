#pragma once
//
//  cluster.h
//
//  Created by Nicholas Hathaway on 7/17/12.

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
#include "njhseq/objects/seqObjects/readObject.hpp"
#include "njhseq/objects/seqObjects/Clusters/baseCluster.hpp"
#include "njhseq/readVectorManipulation/readVectorOperations.h"
#include "njhseq/alignment/aligner.h"

#include "njhseq/IO/SeqIO/SeqOutput.hpp"


namespace njhseq {
// class simulation::errorProfile;

class cluster : public baseCluster {
 public:
  // constructors
  cluster();

  cluster(const seqInfo &firstRead);


  void addRead(const cluster &read);




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
    for (const auto readPos : iter::range(len(inReads))) {
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

	struct snpBreakoutPars {
		uint32_t minSnps = 1;
		double snpFreqCutOff = 0.01;
		double hardSnpFreqCutOff = 0; //! range 0-1, only use snps with frequencies above this number
  	uint32_t hardCutOff = 3; //! there needs to be at least this many reads with the snp to count
  	double clusterFreqCutOff = 0.01; //! should be greater than this internal cluster frequency to breakout a new cluster
  	double strandBiasPValueCutOff = 0.01; //! if the cluster has reads from both directions, snps whose reverse complement read count is this unlikely (two sided binomial p-value) given the cluster's composition are not used
		bool useErrorModel = false; //! instead of hardCutOff and hardSnpFreqCutOff, use snps whose read count is significantly more than the cluster's background error rate would give
		double errorModelPValueCutOff = 0.01; //! upper tail binomial p-value cut off for the error model, bonferroni corrected for the number of positions and alternative bases
		double errorModelQuantile = 0.95; //! the background error rate is this quantile of the mismatch rates of every position and alternative base
		double errorModelMinErrorRate = 0.001; //! the background error rate is at least this
		uint32_t errorModelMinCount = 3; //! with the error model a snp still needs at least this many reads
		bool errorModelAsFilter = false; //! with the error model, snps also have to pass hardCutOff and hardSnpFreqCutOff, so the error model only removes snps
		QualScorePars qScorePars;
	};
  std::vector<cluster> breakoutClustersBasedOnSnps(aligner & alignerObj, const snpBreakoutPars& pars );

	/**@brief Whether a snp passes the strand composition check. If the cluster has reads from both directions, the snp's reverse complement read count is tested against the cluster's reverse complement fraction with a two sided exact binomial test and fails if it's significantly biased to one direction
	 *
	 * @param snpReadCount the read count supporting the snp
	 * @param snpCompReadCount the read count supporting the snp from reverse complemented (_Comp) reads
	 * @param compReadsInCluster the read count of reverse complemented reads in the cluster
	 * @param readsInCluster the read count of the cluster
	 * @param pValueCutOff the snp fails if the binomial p-value is below this
	 * @return true if it passes
	 */
	static bool snpPassesStrandComposition(double snpReadCount, double snpCompReadCount,
			double compReadsInCluster, double readsInCluster, double pValueCutOff);

	/**@brief Whether a read was reverse complemented during extraction. The extractors append _Comp and later steps can append a count (_t) or fraction (_f) after it, so the name has to end in _Comp or _Comp followed by one _t or _f number
	 *
	 * @param name the read name
	 * @return true if the name ends in _Comp, optionally followed by _t<number> or _f<number>
	 */
	static bool readIsRevComp(const std::string & name);

	/**@brief The background error rate for the snp error model, a quantile of the mismatch rates of every possible position and alternative base
	 *
	 * @param mismatchRates the rates of the observed mismatches, the possible mismatches not given are counted as 0
	 * @param possibleMismatches the number of possible mismatches (3 times the length of the sequence)
	 * @param quantile the quantile to use
	 * @param minErrorRate the rate returned is at least this
	 * @return the background error rate
	 */
	static double snpBackgroundErrorRate(std::vector<double> mismatchRates, uint64_t possibleMismatches,
			double quantile, double minErrorRate);

	/**@brief Whether a snp's read count is significantly more than errors alone would give, upper tail binomial test
	 *
	 * @param snpReadCount the read count with the snp
	 * @param readCount the read count of the cluster
	 * @param errorRate the background error rate
	 * @param pValueCutOff passes if the p-value is below this (already corrected for multiple testing)
	 * @return true if it passes
	 */
	static bool snpAboveBackgroundError(double snpReadCount, double readCount, double errorRate, double pValueCutOff);

	/**@brief Whether a snp passes the breakout's read count cut offs: hardCutOff and hardSnpFreqCutOff, or with useErrorModel the error model instead, or both when errorModelAsFilter is also set
	 *
	 * @param snpReadCount the read count with the snp
	 * @param readCount the read count of the cluster
	 * @param pars the breakout parameters
	 * @param backgroundErrorRate the error model's background error rate
	 * @param errorModelPValueCutOff the error model's p-value cut off (already corrected for multiple testing)
	 * @return true if it passes
	 */
	static bool snpPassesBreakoutCutOffs(double snpReadCount, double readCount, const snpBreakoutPars & pars,
			double backgroundErrorRate, double errorModelPValueCutOff);

	using size_type = baseReadObject::size_type;

	virtual Json::Value toJson() const;


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
	cluster mainCluster(getSeqBase(*reads.begin()));
	if(reads.size() > 1){
		std::vector<cluster> inClusters;
		for(const auto readPos : iter::range<uint64_t>(1,reads.size())){
			inClusters.emplace_back(cluster(getSeqBase(reads[readPos])));
		}
		for(const auto & clus : inClusters){
			mainCluster.addRead(clus);
		}
	}
	mainCluster.calculateConsensus(alignerObj, true);
	mainCluster.setName(name);
	return mainCluster;
}
}  // namespace njhseq


