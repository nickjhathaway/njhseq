#pragma once
/*
 * PopGenCalcs.hpp
 *
 *  Created on: Mar 15, 2018
 *      Author: nick
 */
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
#include "njhseq/common.h"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

#include <random>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions.hpp>
#include <utility>

namespace njhseq {


class PopGenCalculator{
public:

	struct TajimaTestRes{
		TajimaTestRes(double d, double pvalnorm, double pvalbeta) :
				d_(d), pval_normal_(pvalnorm), pval_beta_(pvalbeta) {

		}
		TajimaTestRes() = default;
		double d_{std::numeric_limits<double>::max()};
		double pval_normal_{std::numeric_limits<double>::max()};
		double pval_beta_{std::numeric_limits<double>::max()};
	};


	/**@brief Calculate tajima's test for neutrality
	 *
	 * @param nInputSeqs The number of input sequences
	 * @param nSegragtingSites The number of segregating sites
	 * @param meanPairwiseDifferences The average number of difference between input sequences
	 * @return
	 */
	static TajimaTestRes calcTajimaTest(uint32_t nInputSeqs, uint32_t nSegragtingSites, double meanPairwiseDifferences);

	struct ExpectedPloidyInfo {

		uint32_t ploidy_;

		long double expectedPolyClonal_; //!< the expected freq of polyclonal samples for given ploidy for given population frequencies
		std::unordered_map<uint32_t, long double> expectedCOIForPloidy_; //!< the expected COI or given ploidy for given population frequencies

		//currently only does ploidy up to and including 5, will throw otherwise
		static ExpectedPloidyInfo genPloidyInfo(uint32_t ploidy, const std::vector<long double> & freqs);
	};


	struct ExpectedKHeterozygosityRes {
		ExpectedKHeterozygosityRes() = default;
		ExpectedKHeterozygosityRes(uint32_t k, const std::vector<long double> & freqs);
		uint32_t k_{std::numeric_limits<uint32_t>::max()};
		long double expected_monoclonal_{std::numeric_limits<long double>::max()};
		long double k_heterozygosity_{std::numeric_limits<long double>::max()};

		static long double factorial_int(uint32_t k);
	};



	struct DiversityMeasures {

		uint32_t alleleNumber_ = 0; //!< number of unique alleles
		uint32_t doublets_ = 0; //!< number of haplotypes found twice
		uint32_t singlets_ = 0; //!< number of haplotypes found only once
		double expShannonEntropy_ = std::numeric_limits<double>::max(); //!< exp of shannon entropy base e
		double ShannonEntropyE_ = std::numeric_limits<double>::max(); //!< shannon entropy base e
		double effectiveNumOfAlleles_ = std::numeric_limits<double>::max();//!< effective number of alleles
		double heterozygostiy_  = std::numeric_limits<double>::max();//!< the expected heterozygostity (He)

		double simpsonIndex_ = std::numeric_limits<double>::max(); //!< simpson index of diversity

		std::map<uint32_t, ExpectedKHeterozygosityRes> expected_k_heterozygosities;//!< probability of choosing a specific number of unique haplotypes, key is the number of expected;

	};






	struct PopDifferentiationMeasures{

		std::unordered_map<std::string, double> hjsSample_;
		double hsSample_ = std::numeric_limits<double>::max();
		double htSample_ = std::numeric_limits<double>::max();

		double hsEst_ = std::numeric_limits<double>::max();
		double htEst_ = std::numeric_limits<double>::max();

		double gst_ = std::numeric_limits<double>::max();
		double jostD_ = std::numeric_limits<double>::max();

		double gstEst_ = std::numeric_limits<double>::max();
		double jostDEst_ = std::numeric_limits<double>::max();
		double chaoA_ = std::numeric_limits<double>::max();
		double chaoB_ = std::numeric_limits<double>::max();
		double jostDChaoEst_ = std::numeric_limits<double>::max();


		double informativenessForAssign_ = 0;

		std::unordered_map<uint32_t, double> informativenessForAssignPerHap_;
		std::unordered_map<std::string, double> informativenessForAssignPerPopulation_;

	};

	struct PopDifferentiationMeasuresPairWise{
		PopDifferentiationMeasuresPairWise(
				PopDifferentiationMeasures  genDiffMeasures,
				std::string  pop1Name,
				std::string  pop2Name): genDiffMeasures_(std::move(genDiffMeasures)),
						pop1Name_(std::move(pop1Name)),
						pop2Name_(std::move(pop2Name)){

		}
		PopDifferentiationMeasuresPairWise() = default;
		PopDifferentiationMeasures genDiffMeasures_;
		std::string pop1Name_;
		std::string pop2Name_;
		//Specific to comparing two populations
		double brayCurtisDissim_ = std::numeric_limits<double>::max();
		double brayCurtisRelativeDissim_ = std::numeric_limits<double>::max();
		//double chiSquare_ = std::numeric_limits<double>::max();
		double jaccardIndexDissim_ = std::numeric_limits<double>::max();
		double sorensenDistance_ = std::numeric_limits<double>::max();

		double matchingCoefficientDistance_ = std::numeric_limits<double>::max();
		double halfR_ = std::numeric_limits<double>::max();//!< 0.5 * 1 - correlation

		double RMSE_ = std::numeric_limits<double>::max();

		//Avalanche

		double discriminatingAvalance_ = std::numeric_limits<double>::max(); //!< not yet implemented, would take a "genetic" distance into account as well
		double plainAvalance_ = std::numeric_limits<double>::max();



		uint32_t uniqueHapsAll_ { 0 };
		uint32_t uniqueHapsShared_ { 0 };
		uint32_t uniqueHapsInPop1_ { 0 };
		uint32_t uniqueHapsInPop2_ { 0 };

		double uniqueHapsInPop1CumFreq_ {0};
		double uniqueHapsInPop2CumFreq_ {0};


	};


	struct PopHapInfo {
		PopHapInfo(const uint32_t& popUid, uint32_t count) : popUid_(popUid), unweighted_count_(count) {
		}

		PopHapInfo(const uint32_t& popUid, uint32_t count, double prob) : popUid_(popUid), unweighted_count_(count),
		                                                                  unweighted_prob_(prob) {
		}

		uint32_t popUid_;
		uint32_t unweighted_count_;
		double weighted_count_ = 0;

		double unweighted_prob_{0};
		double weighted_prob_{0};

		static double getTotalWeightedPopCount(const std::vector<PopHapInfo>& hapsForPopulation) {
			return std::accumulate(hapsForPopulation.begin(), hapsForPopulation.end(), 0,
			                       [](uint32_t weighted_count, const PopHapInfo& hap) {
				                       return weighted_count + hap.weighted_count_;
			                       });
		}

		static uint32_t getTotalPopCount(const std::vector<PopHapInfo>& hapsForPopulation) {
			return std::accumulate(hapsForPopulation.begin(), hapsForPopulation.end(), 0,
			                       [](uint32_t total, const PopHapInfo& hap) {
				                       return total + hap.unweighted_count_;
			                       });
		}

		static void setProb(std::vector<PopHapInfo>& hapsForPopulation, uint32_t total) {
			njh::for_each(hapsForPopulation, [&total](PopHapInfo& hap) {
				hap.unweighted_prob_ = hap.unweighted_count_ / static_cast<double>(total);
			});
		}

		static void setProb(std::vector<PopHapInfo>& hapsForPopulation) {
			auto total = getTotalPopCount(hapsForPopulation);
			njh::for_each(hapsForPopulation, [&total](PopHapInfo& hap) {
				hap.unweighted_prob_ = hap.unweighted_count_ / static_cast<double>(total);
			});
		}

		static void setProbWeighted(std::vector<PopHapInfo>& hapsForPopulation, double weighted_total) {
			njh::for_each(hapsForPopulation, [&weighted_total](PopHapInfo& hap) {
				hap.weighted_prob_ = hap.weighted_count_ / weighted_total;
			});
		}

		static void setProbWeighted(std::vector<PopHapInfo>& hapsForPopulation) {
			auto weighted_total = getTotalWeightedPopCount(hapsForPopulation);
			njh::for_each(hapsForPopulation, [&weighted_total](PopHapInfo& hap) {
				hap.weighted_prob_ = hap.weighted_count_ / weighted_total;
			});
		}
	};


		/**@brief Get several general measures of diversity, assumes haps are already collapsed to unique haplotypes and have frequencies set
	 *
	 * @param haps a vector of unique haplotypes
	 * @param by_unweighted_counts by default, the frequency metrics are calculated by the weighted (frac_) counts, set this to true to do freqs by count (cnt_)
	 * @return a struct with several diversity measurements
	 */
	template<typename T>
	static DiversityMeasures getGeneralMeasuresOfDiversity(const std::vector<T> & haps, bool by_unweighted_counts = false){

		std::unordered_map<std::string, uint32_t> popCounts;
		std::unordered_map<std::string, uint32_t> popCountsWeighted;
		for(const auto & seq : haps){
			popCounts[getSeqBase(seq).seq_] += getSeqBase(seq).cnt_;
		}
		for(const auto & seq : haps){
			popCountsWeighted[getSeqBase(seq).seq_] += getSeqBase(seq).frac_;
		}
		std::vector<PopGenCalculator::PopHapInfo> popHapInfos;
		uint32_t count = 0;
		for(const auto & popCount : popCounts){
			popHapInfos.emplace_back(count, popCount.second);
			popHapInfos[count].weighted_count_ = popCountsWeighted[popCount.first];
			++count;
		}
		return getGeneralMeasuresOfDiversity(popHapInfos, by_unweighted_counts);
	}

	/**@brief Get several general measures of diversity,
	 *
	 * @param haps a vector of unique haplotypes
	 * @return a struct with several diversity measurements
	 */
	template<typename T>
	static DiversityMeasures getGeneralMeasuresOfDiversityRawInput(const std::vector<T> & haps){
		std::unordered_map<std::string, uint32_t> popCounts;
		for(const auto & seq : haps){
			++popCounts[getSeqBase(seq).seq_];
		}
		std::vector<PopGenCalculator::PopHapInfo> popHapInfos;
		uint32_t count = 0;
		for(const auto & popCount : popCounts){
			popHapInfos.emplace_back(count, popCount.second);
			//set weighted_count_ as well in case used later
			popHapInfos[count].weighted_count_ = popCount.second;
			++count;
		}
		return getGeneralMeasuresOfDiversity(popHapInfos);
	}

	/**
	 * @brief get general measurements of diversity
	 * @param haps the population haplotypes with counts
	 * @param by_unweighted_counts by default, the frequency metrics are calculated by the weighted (frac_) counts, set this to true to do freqs by count (cnt_)
	 * @return a struct with several diversity metrics
	 */
	static DiversityMeasures getGeneralMeasuresOfDiversity(const std::vector<PopHapInfo> & haps, bool by_unweighted_counts = false);


	static PopDifferentiationMeasuresPairWise getPopDiffUnweighted(
			const std::string & pop1, const std::vector<PopHapInfo> & pop1Haps,
			const std::string & pop2, const std::vector<PopHapInfo> & pop2Haps,
			const std::unordered_set<uint32_t> & allPossibleHaps,
			const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDistance = std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>{});

	static PopDifferentiationMeasuresPairWise getPopDiffWeighted(
		const std::string & pop1, const std::vector<PopHapInfo> & pop1Haps,
		const std::string & pop2, const std::vector<PopHapInfo> & pop2Haps,
		const std::unordered_set<uint32_t> & allPossibleHaps,
		const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDistance = std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>{});


	template<typename T>
	static PopDifferentiationMeasures getOverallPopDiffForSeqsWeighted(
		const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>>& popSeqs) {
		if (popSeqs.size() < 2) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::unordered_map<std::string, std::vector<PopHapInfo>> hapsForPopulations;
		std::unordered_map<std::string, uint32_t> seqCounts;
		for (const auto& pop: popSeqs) {
			for (const auto& hap: *pop.second) {
				seqCounts[getSeqBase(hap).seq_] += getSeqBase(hap).cnt_;
			}
		}
		auto seqs = njh::getVecOfMapKeys(seqCounts);
		njh::sort(seqs, [&seqCounts](const std::string& seq1, const std::string& seq2) {
			return seqCounts[seq1] > seqCounts[seq2];
		});

		std::unordered_map<std::string, uint32_t> seqToPopUID;
		for (const auto pos: iter::range(seqs.size())) {
			seqToPopUID[seqs[pos]] = pos;
		}

		for (const auto& pop: popSeqs) {
			for (const auto& hap: *pop.second) {
				hapsForPopulations[pop.first].emplace_back(PopHapInfo(seqToPopUID[getSeqBase(hap).seq_], getSeqBase(hap).cnt_));
				//add weighted count as well
				hapsForPopulations[pop.first].back().weighted_count_ = getSeqBase(hap).frac_;
			}
		}
		return getOverallPopDiffWeighted(hapsForPopulations);
	}

	template<typename T>
	static PopDifferentiationMeasures getOverallPopDiffForSeqsUnweighted(
		const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>>& popSeqs) {
		if (popSeqs.size() < 2) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::unordered_map<std::string, std::vector<PopHapInfo>> hapsForPopulations;
		std::unordered_map<std::string, uint32_t> seqCounts;
		for (const auto& pop: popSeqs) {
			for (const auto& hap: *pop.second) {
				seqCounts[getSeqBase(hap).seq_] += getSeqBase(hap).cnt_;
			}
		}
		auto seqs = njh::getVecOfMapKeys(seqCounts);
		njh::sort(seqs, [&seqCounts](const std::string& seq1, const std::string& seq2) {
			return seqCounts[seq1] > seqCounts[seq2];
		});

		std::unordered_map<std::string, uint32_t> seqToPopUID;
		for (const auto pos: iter::range(seqs.size())) {
			seqToPopUID[seqs[pos]] = pos;
		}

		for (const auto& pop: popSeqs) {
			for (const auto& hap: *pop.second) {
				hapsForPopulations[pop.first].emplace_back(PopHapInfo(seqToPopUID[getSeqBase(hap).seq_], getSeqBase(hap).cnt_));
				//add weighted count as well
				hapsForPopulations[pop.first].back().weighted_count_ = getSeqBase(hap).frac_;
			}
		}
		return getOverallPopDiffUnweighted(hapsForPopulations);
	}



	static PopDifferentiationMeasures getOverallPopDiffWeighted(std::unordered_map<std::string, std::vector<PopHapInfo> > hapsForPopulations);
	static PopDifferentiationMeasures getOverallPopDiffUnweighted(std::unordered_map<std::string, std::vector<PopHapInfo> > hapsForPopulations);

	template<typename T>
	static std::unordered_map<std::string,
			std::unordered_map<std::string, PopDifferentiationMeasuresPairWise>> getPairwisePopDiffUnweighted(
			const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> & popSeqs,
			const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDistance = std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>{}) {
		if(popSeqs.size() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
			throw std::runtime_error{ss.str()};
		}

		auto keys = njh::getVecOfMapKeys(popSeqs);
		njh::sort(keys);
		std::unordered_map<std::string, uint32_t> seqCounts;
		for(const auto & pop : popSeqs){
			for(const auto & hap : *pop.second){
				seqCounts[getSeqBase(hap).seq_] += getSeqBase(hap).cnt_;
			}
		}
		auto seqs = njh::getVecOfMapKeys(seqCounts);
		njh::sort(seqs,[&seqCounts](const std::string & seq1, const std::string & seq2){
			return seqCounts[seq1] > seqCounts[seq2];
		});

		std::unordered_map<std::string, uint32_t> seqToPopUID;
		for(const auto pos : iter::range(seqs.size())){
			seqToPopUID[seqs[pos]] = pos;
		}

		std::unordered_map<std::string, std::vector<PopHapInfo> > hapsForPopulations;
		for(const auto & pop : popSeqs){
			for(const auto & hap : *pop.second){
				hapsForPopulations[pop.first].emplace_back(PopHapInfo(seqToPopUID[getSeqBase(hap).seq_], getSeqBase(hap).cnt_));
				hapsForPopulations[pop.first].back().weighted_count_ = getSeqBase(hap).frac_;
			}
		}
		//pairwiseDistance
		return getPairwisePopDiffUnweighted(hapsForPopulations, pairwiseDistance);
	}
	static std::unordered_map<std::string,
		std::unordered_map<std::string, PopDifferentiationMeasuresPairWise>> getPairwisePopDiffUnweighted(
		const std::unordered_map<std::string, std::vector<PopHapInfo>> & hapsForPopulations,
		const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDists = std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>{});

	template<typename T>
	static std::unordered_map<std::string,
		std::unordered_map<std::string, PopDifferentiationMeasuresPairWise>> getPairwisePopDiffWeighted(
		const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>>& popSeqs,
		const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>& pairwiseDistance = std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>{}) {

		if (popSeqs.size() < 2) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
			throw std::runtime_error{ss.str()};
		}

		auto keys = njh::getVecOfMapKeys(popSeqs);
		njh::sort(keys);
		std::unordered_map<std::string, uint32_t> seqCounts;
		for (const auto& pop: popSeqs) {
			for (const auto& hap: *pop.second) {
				seqCounts[getSeqBase(hap).seq_] += getSeqBase(hap).cnt_;
			}
		}
		auto seqs = njh::getVecOfMapKeys(seqCounts);
		njh::sort(seqs, [&seqCounts](const std::string& seq1, const std::string& seq2) {
			return seqCounts[seq1] > seqCounts[seq2];
		});

		std::unordered_map<std::string, uint32_t> seqToPopUID;
		for (const auto pos: iter::range(seqs.size())) {
			seqToPopUID[seqs[pos]] = pos;
		}

		std::unordered_map<std::string, std::vector<PopHapInfo>> hapsForPopulations;
		for (const auto& pop: popSeqs) {
			for (const auto& hap: *pop.second) {
				hapsForPopulations[pop.first].emplace_back(PopHapInfo(seqToPopUID[getSeqBase(hap).seq_], getSeqBase(hap).cnt_));
				hapsForPopulations[pop.first].back().weighted_count_ = getSeqBase(hap).frac_;
			}
		}
		//pairwiseDistance
		return getPairwisePopDiffWeighted(hapsForPopulations, pairwiseDistance);
	}

	static std::unordered_map<std::string,
		std::unordered_map<std::string, PopDifferentiationMeasuresPairWise>> getPairwisePopDiffWeighted(
		const std::unordered_map<std::string, std::vector<PopHapInfo>>& hapsForPopulations,
		const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>& pairwiseDists = std::unordered_map<
			uint32_t, std::unordered_map<uint32_t, double>>{});


	class FisherExactFor2x2{
	public:

		struct FisherExactFor2x2Result{

			double oddsRatio_ {std::numeric_limits<double>::max()};
			double pValue_ {std::numeric_limits<double>::max()};
			double lowerConfInterval_ {std::numeric_limits<double>::max()};
			double upperConfInterval_ {std::numeric_limits<double>::max()};

		};

		struct FisherExactFor2x2Input{
			enum class Tail{
				TWOTAILED,
				GREATER,
				LESSER
			};
			double confInterval{0.95};/**< should be in (0,1) */

			//contingency table
			//     [,1] [,2]
			// [1,] TP   FN
			// [2,] FP   TN

			//     [,1] [,2]
			// [1,] a    b
			// [2,] c    d

			uint32_t TP{0};
			uint32_t FN{0};
			uint32_t FP{0};
			uint32_t TN{0};

			[[nodiscard]] uint32_t total() const {
				return  TP + FP + FN + TN;
			}



			Tail pvalTail{Tail::TWOTAILED};

		};

		static FisherExactFor2x2Result runFisherExactOn2x2(const FisherExactFor2x2Input & inPars);

		static std::vector<int32_t> genSupportVec(int32_t hi, int32_t lo){
			std::vector<int32_t> support(hi + 1 - lo);
			njh::iota(support, lo);
			return support;
		}
	};

};


}  // namespace njhseq





