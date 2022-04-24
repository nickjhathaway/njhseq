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

namespace njhseq {


class PopGenCalculator{
public:

	struct TajimaTestRes{
		TajimaTestRes(double d, double pvalnorm, double pvalbeta) :
				d_(d), pval_normal_(pvalnorm), pval_beta_(pvalbeta) {

		}
		TajimaTestRes() {

		}
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

		double expectedPolyClonal_; //!< the expected freq of polyclonal samples for given ploidy for given population frequencies
		std::unordered_map<uint32_t, double> expectedCOIForPloidy_; //!< the expected COI or given ploidy for given population frequencies


		double getMaxExpPloidy() const;

		//currently only does ploidy up to and including 5, will throw otherwise
		static ExpectedPloidyInfo genPloidyInfo(uint32_t ploidy, const std::vector<long double> & freqs);
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

		ExpectedPloidyInfo ploidy2_;//!< info when sampling 2 haplotypes
		ExpectedPloidyInfo ploidy3_;//!< info when sampling 3 haplotypes
		ExpectedPloidyInfo ploidy4_;//!< info when sampling 4 haplotypes
		ExpectedPloidyInfo ploidy5_;//!< info when sampling 5 haplotypes



	};


	/**@brief Get several general measures of diversity, assumes haps are already collapsed to unique haplotypes and have frequencies set
	 *
	 * @param haps a vector of unique haplotypes
	 * @return a struct with several diversity measurements
	 */
	template<typename T>
	static DiversityMeasures getGeneralMeasuresOfDiversity(const std::vector<T> & haps){

		std::unordered_map<std::string, uint32_t> popCounts;
		for(const auto & seq : haps){
			popCounts[getSeqBase(seq).seq_] += getSeqBase(seq).cnt_;
		}
		std::vector<PopGenCalculator::PopHapInfo> popHapInfos;
		uint32_t count = 0;
		for(const auto & popCount : popCounts){
			popHapInfos.emplace_back(count, popCount.second);
			++count;
		}
		return getGeneralMeasuresOfDiversity(popHapInfos);
//
//
//		DiversityMeasures res;
//
//		res.alleleNumber_ = haps.size();
//		double sumOfSquares = 0;
//		double sumOfLogFreqTimesFreq = 0;
//
//
//		for (const auto & hap : haps) {
//			const seqInfo & seqRef = getSeqBase(hap);
//			if (1 == seqRef.cnt_) {
//				++res.singlets_;
//			} else if (2 == seqRef.cnt_) {
//				++res.doublets_;
//			}
//			sumOfSquares += std::pow(seqRef.frac_, 2.0);
//			sumOfLogFreqTimesFreq += seqRef.frac_ * std::log(seqRef.frac_);
//		}
//
//		res.heterozygostiy_ = 1 - sumOfSquares;
//		res.effectiveNumOfAlleles_ = std::pow(sumOfSquares, -1);
//		res.ShannonEntropyE_ = -sumOfLogFreqTimesFreq;
//		res.expShannonEntropy_ = std::exp(-sumOfLogFreqTimesFreq);
//
//		return res;
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
			++count;
		}
		return getGeneralMeasuresOfDiversity(popHapInfos);
	}




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
				const PopDifferentiationMeasures & genDiffMeasures,
				const std::string & pop1Name,
				const std::string & pop2Name): genDiffMeasures_(genDiffMeasures),
						pop1Name_(pop1Name),
						pop2Name_(pop2Name){

		}
		PopDifferentiationMeasuresPairWise(){
			//here for convenience
		}
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


	//
	//template<typename T>
	//PopDifferentiationMeasures getOverallPopDiff(const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> & popSeqs){
	//	if(popSeqs.size() < 2){
	//		std::stringstream ss;
	//		ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
	//		throw std::runtime_error{ss.str()};
	//	}
	//	PopDifferentiationMeasures ret;
	//
	//	std::unordered_map<std::string, uint32_t> subPopSizes;
	//	std::unordered_set<std::string> allHapSeqs;
	//	std::unordered_map<std::string, std::unordered_map<std::string, double>> freqsForPopForHapSeq;
	//	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsForPopForHapSeq;
	//
	//	for(const auto & subPop : popSeqs){
	//		subPopSizes[subPop.first] =
	//				std::accumulate(subPop.second->begin(), subPop.second->end(), 0,
	//						[](uint32_t total, const T & seq) {
	//							return total + std::round(getSeqBase(seq).cnt_);
	//						});
	//		double sumOfSquares = 0;
	//		for(const auto & seq : *subPop.second){
	//			allHapSeqs.emplace(getSeqBase(seq).seq_);
	//			freqsForPopForHapSeq[subPop.first][getSeqBase(seq).seq_] = getSeqBase(seq).frac_;
	//			countsForPopForHapSeq[subPop.first][getSeqBase(seq).seq_] = std::round(getSeqBase(seq).cnt_);
	//			sumOfSquares += std::pow(getSeqBase(seq).frac_, 2.0);
	//		}
	//		ret.hjsSample_[subPop.first] = 1 - sumOfSquares;
	//	}
	//
	//	double sumOfHjsSample = 0;
	//	for(const auto & subPop : ret.hjsSample_){
	//		sumOfHjsSample += subPop.second;
	//	}
	//	ret.hsSample_ = sumOfHjsSample/ret.hjsSample_.size();
	//	double jtSample = 0;
	//	for(const auto & hapSeq : allHapSeqs){
	//		double freqSum = 0;
	//		for( auto & subPop : freqsForPopForHapSeq){
	//			freqSum += subPop.second[hapSeq];
	//		}
	//		jtSample += std::pow(freqSum/popSeqs.size(), 2.0);
	//	}
	//	ret.htSample_ = 1 - jtSample;
	//
	//	double harmonicMean = 0;
	//	double sumOfInverses = 0;
	//	uint32_t totalHaps = 0;
	//	for(const auto & popSize : subPopSizes){
	//		totalHaps += popSize.second;
	//		sumOfInverses += 1.0/popSize.second;
	//	}
	//	harmonicMean = subPopSizes.size()/sumOfInverses;
	//
	//	ret.hsEst_ = (harmonicMean/(harmonicMean - 1)) * ret.hsSample_;
	//	ret.htEst_ = ret.htSample_ + (ret.hsEst_)/(harmonicMean * popSeqs.size());
	//
	//	ret.gst_ = (ret.htSample_ - ret.hsSample_)/ret.htSample_;
	//	ret.jostD_ = ((ret.htSample_ - ret.hsSample_)/(1 - ret.hsSample_)) * (popSeqs.size()/(popSeqs.size() - 1));
	//
	//	ret.gstEst_ =  (ret.htEst_ - ret.hsEst_)/ret.htEst_;
	//	ret.jostDEst_ = ((ret.htEst_ - ret.hsEst_)/(1 - ret.hsEst_)) * (popSeqs.size()/(popSeqs.size() - 1));
	//
	//
	//	double a = 0;
	//	for(const auto & hapSeq : allHapSeqs){
	//		double sumOfFreqs = 0;
	//		double sumOfSqaureFreqs = 0;
	//		for( auto & subPop : freqsForPopForHapSeq){
	//			sumOfFreqs += subPop.second[hapSeq];
	//			sumOfSqaureFreqs += std::pow(subPop.second[hapSeq], 2.0);
	//		}
	//		a += (std::pow(sumOfFreqs, 2.0) - sumOfSqaureFreqs)/(subPopSizes.size() - 1);
	//	}
	//	ret.chaoA_ = a;
	//	double b = 0;
	//	for(const auto & hapSeq : allHapSeqs){
	//		for(auto & subPop : countsForPopForHapSeq){
	//			if(subPop.second[hapSeq] > 0){
	//				b += (subPop.second[hapSeq] *(subPop.second[hapSeq] - 1) )/static_cast<double>(subPopSizes[subPop.first] * (subPopSizes[subPop.first] - 1));
	//			}
	//		}
	//	}
	//	ret.chaoB_ = b;
	//	ret.jostDChaoEst_ = 1 - (a/b);
	//	return ret;
	//}
	//


	struct PopHapInfo {
		PopHapInfo(const uint32_t & popUid, uint32_t count): popUid_(popUid), count_(count){

		}
		PopHapInfo(const uint32_t & popUid, uint32_t count, double prob): popUid_(popUid), count_(count), prob_(prob){

		}

		uint32_t popUid_;
		uint32_t count_;

		double prob_{0};


		static uint32_t getTotalPopCount(const std::vector<PopHapInfo> & hapsForPopulation){
			return std::accumulate(hapsForPopulation.begin(), hapsForPopulation.end(), 0, [](uint32_t total, const PopHapInfo & hap){
				return total + hap.count_;
			});
		}

		static void setProb(std::vector<PopHapInfo> & hapsForPopulation, uint32_t total){
			njh::for_each(hapsForPopulation, [&total](PopHapInfo & hap){
				hap.prob_ = hap.count_/static_cast<double>(total);
			});
		}
		static void setProb(std::vector<PopHapInfo> & hapsForPopulation){
			auto total = getTotalPopCount(hapsForPopulation);
			njh::for_each(hapsForPopulation, [&total](PopHapInfo & hap){
				hap.prob_ = hap.count_/static_cast<double>(total);
			});
		}
	};


	static DiversityMeasures getGeneralMeasuresOfDiversity(const std::vector<PopHapInfo> & haps);


	static PopDifferentiationMeasures getOverallPopDiff(std::unordered_map<std::string, std::vector<PopHapInfo> > hapsForPopulations);


	static PopDifferentiationMeasuresPairWise getPopDiff(
			const std::string & pop1, const std::vector<PopHapInfo> & pop1Haps,
			const std::string & pop2, const std::vector<PopHapInfo> & pop2Haps,
			const std::unordered_set<uint32_t> & allPossibleHaps,
			const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDistacne = std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>{});

	template<typename T>
	static PopDifferentiationMeasures getOverallPopDiffForSeqs(const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> & popSeqs){
		if(popSeqs.size() < 2){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " error, popSeqs should at least be size 2 not " << popSeqs.size() << "\n";
			throw std::runtime_error{ss.str()};
		}
		std::unordered_map<std::string, std::vector<PopHapInfo> > hapsForPopulations;
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

		for(const auto & pop : popSeqs){
			for(const auto & hap : *pop.second){
				hapsForPopulations[pop.first].emplace_back(PopHapInfo(seqToPopUID[getSeqBase(hap).seq_], getSeqBase(hap).cnt_));
			}
		}
		return getOverallPopDiff(hapsForPopulations);
	}


	static std::unordered_map<std::string,
			std::unordered_map<std::string, PopDifferentiationMeasuresPairWise>> getPairwisePopDiff(
			const std::unordered_map<std::string, std::vector<PopHapInfo>> & hapsForPopulations,
			const std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> & pairwiseDists = std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>>{});


	template<typename T>
	static std::unordered_map<std::string,
			std::unordered_map<std::string, PopDifferentiationMeasuresPairWise>> getPairwisePopDiff(
			const std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> & popSeqs,
			const std::vector<std::vector<double>> & pairwiseDistacne = std::vector<std::vector<double>>{}) {
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
			}
		}
		return getPairwisePopDiff(hapsForPopulations);
//		for(const auto keyPos : iter::range(keys.size())){
//			for(const auto secondKeyPos : iter::range(keyPos)){
//	//			std::unordered_map<std::string, std::shared_ptr<std::vector<T>>> currentPair;
//	//			currentPair[keys[keyPos]] = popSeqs.at(keys[keyPos]);
//	//			currentPair[keys[secondKeyPos]] = popSeqs.at(keys[secondKeyPos]);
//	//			auto popMeasures = getOverallPopDiff(currentPair);
//				auto popMeasures = getPopDiff(keys[keyPos], hapsForPopulations[keys[keyPos]],
//						keys[secondKeyPos], hapsForPopulations[keys[secondKeyPos]]);
//				ret[keys[keyPos]][keys[secondKeyPos]] = popMeasures;
//			}
//		}
//		return ret;
	}

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

		static FisherExactFor2x2Result runFisherExactOn2x2(const FisherExactFor2x2Input & inPars){
			//implementation based on R implement of 2x2 fisher exact
			uint32_t k = inPars.TP + inPars.FN;
			uint32_t m = inPars.TP + inPars.FP;
			uint32_t n = inPars.FN + inPars.TN;
			uint32_t total = inPars.total();
			uint32_t lo = n > k ? 0U : k - n; //std::max(0, k - n);
			uint32_t hi = std::min(k, m);

			std::vector<uint32_t> support(hi + 1 - lo);
			njh::iota(support, lo);
			//m, k, total
			boost::math::hypergeometric_distribution<double> hg_dist(m, k, total);
			std::vector<double> logdc;
			logdc.reserve(support.size())	;
			for(const auto val : support){
				logdc.emplace_back(log(boost::math::pdf(hg_dist, val)));
			}
			auto dnhyper = [&logdc,&support](double ncp){
				std::vector<double> ret(logdc.size());
				for(const auto idx : iter::range(logdc.size())){
					ret[idx] = logdc[idx] + log(ncp) * support[idx];
				}
				auto maxret = vectorMaximum(ret);
				for(const auto idx : iter::range(logdc.size())){
					ret[idx] = exp(ret[idx] - maxret);
				}
				auto sumret = vectorSum(ret);
				for(const auto idx : iter::range(logdc.size())){
					ret[idx] = ret[idx]/sumret;
				}
				return ret;
			};
//		auto test_dnhyper = dnhyper(1);
//		std::cout << "test_dnhyper: "  << std::endl;
//		for(const auto idx : iter::range(test_dnhyper.size())){
//			std::cout << test_dnhyper[idx] << std::endl;
//		}

			auto mnhyper = [&lo,&hi,&dnhyper,&support](double ncp){
				if (ncp == 0){
					return static_cast<double>(lo);
				}
				if (ncp == std::numeric_limits<double>::infinity()){
					return static_cast<double>(hi);
				}
				double sum = 0;
				auto dnhyper_res = dnhyper(ncp);
				for(const auto idx : iter::range(support.size())){
					sum += dnhyper_res[idx] * support[idx];
				}
				return sum;
			};


			auto pnhyper = [&lo,&hi,&dnhyper,&support,&hg_dist,&inPars](double q, double ncp =1, bool upperTail = false)  {
				if (ncp == 1) {
					if(upperTail){
						//1 - phyper(..., lower.tail = TRUE) is the same as phyper(..., lower.tail = FALSE)
						return 1 - boost::math::cdf(hg_dist, inPars.TP - 1);//will fail if TP == 0
					} else {
						return boost::math::cdf(hg_dist, inPars.TP);
					}
				}
				if (ncp == 0) {
					if(upperTail){
						return q <= lo ? 1. : 0.;
					}else{
						return q >= lo ? 1. : 0.;
					}
				}
				if (ncp == std::numeric_limits<double>::infinity()){
					if(upperTail){
						return q <= hi ? 1. : 0.;
					}else{
						return q >= hi ? 1. : 0.;
					}
				}
				double sum = 0;
				auto dnhyper_res = dnhyper(ncp);
				if(upperTail){
					for(const auto idx : iter::range(support.size())){
						if(support[idx] >= q){
							sum += dnhyper_res[idx];
						}
					}
				}else{
					for(const auto idx : iter::range(support.size())){
						if(support[idx] <= q){
							sum += dnhyper_res[idx];
						}
					}
				}
				return sum;
			};

			auto pnhyper_twotailed = [&lo,&hi,&dnhyper,&support,&inPars](const double test_odds)  {
				if (test_odds == 0){
					return inPars.TP == lo ? 1. : 0.;
				} else if (test_odds == std::numeric_limits<double>::infinity()){
					return inPars.TP == hi ? 1. : 0.;
				} else {
					double relErr = 1 + 1e-07;
					auto dnhyper_res = dnhyper(test_odds);
					double sum = 0;
					double testAgainst = dnhyper_res[inPars.TP - lo] * relErr;
					for(const auto idx : iter::range(support.size())){
						if(dnhyper_res[idx] <= testAgainst){
							sum += dnhyper_res[idx];
						}
					}
					return sum;
				}
			};
			auto mle = [&lo,&hi,&mnhyper](const double x) {
				if (x == lo) {
					return 0.0;
				}
				if (x == hi){
					return std::numeric_limits<double>::infinity();
				}
				boost::uintmax_t max_iter=500; // Set max iterations.

				boost::math::tools::eps_tolerance<double> tol(30); //Set the eps tolerance.

//



				double mu = mnhyper(1);
				if (mu > x) {
					auto funcToRoot = [&mnhyper,&x](double input){
						return mnhyper(input) - x;
					};
					auto r1 = boost::math::tools::toms748_solve(funcToRoot, 0.0, 1.0, tol, max_iter); // use the toms solve algorithm.
					return r1.first;
					//uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
				} else if (mu < x) {
					auto funcToRoot = [&mnhyper,&x](double input){
						return mnhyper(1/input) - x;
					};
					auto r1 = boost::math::tools::toms748_solve(funcToRoot, std::numeric_limits<double>::epsilon(), 1.0, tol, max_iter); // use the toms solve algorithm.
					return 1/r1.first;
					//1 / uniroot(function(t) mnhyper(1 / t) - x, c(.Machine$double.eps, 1))$root
				}
				return 1.0;
			};

			auto upperConf = [&pnhyper,&hi](const double x, const double alpha){
				if (x == hi) {
					return std::numeric_limits<double>::infinity();
				}
				boost::uintmax_t max_iter=500; // Set max iterations.

				boost::math::tools::eps_tolerance<double> tol(30); //Set the eps tolerance.

				double p = pnhyper(x, 1);
				if (p < alpha) {
					auto funcToRoot = [&pnhyper,&x,&alpha](double t){
						return pnhyper(x, t) - alpha;
					};
					auto r1 = boost::math::tools::toms748_solve(funcToRoot, 0.0, 1.0, tol, max_iter); // use the toms solve algorithm.
					return r1.first;
					//uniroot(function(t) pnhyper(x, t) - alpha, c(0, 1))$root
				} else if (p > alpha) {
					auto funcToRoot = [&pnhyper,&x,&alpha](double t){
						return pnhyper(x, 1/t) - alpha;
					};
					auto r1 = boost::math::tools::toms748_solve(funcToRoot, std::numeric_limits<double>::epsilon(), 1.0, tol, max_iter); // use the toms solve algorithm.
					return 1/r1.first;
					//1 / uniroot(function(t) pnhyper(x, 1 / t) - alpha, c(.Machine$double.eps, 1))$root
				} else{
					return 1.0;
				}
			};

			auto lowerConf = [&pnhyper,&lo](const double x, const double alpha){
				if (x == lo) {
					return std::numeric_limits<double>::infinity();
				}
				boost::uintmax_t max_iter=500; // Set max iterations.

				boost::math::tools::eps_tolerance<double> tol(30); //Set the eps tolerance.

				double p = pnhyper(x, 1, true);
				if (p > alpha) {
					auto funcToRoot = [&pnhyper,&x,&alpha](double t){
						return pnhyper(x, t, true) - alpha;
					};
					auto r1 = boost::math::tools::toms748_solve(funcToRoot, 0.0, 1.0, tol, max_iter); // use the toms solve algorithm.
					return r1.first;
					//uniroot(function(t) pnhyper(x, t) - alpha, c(0, 1))$root
				} else if (p < alpha) {
					auto funcToRoot = [&pnhyper,&x,&alpha](double t){
						return pnhyper(x, 1/t, true) - alpha;
					};
					auto r1 = boost::math::tools::toms748_solve(funcToRoot, std::numeric_limits<double>::epsilon(), 1.0, tol, max_iter); // use the toms solve algorithm.
					return 1/r1.first;
					//1 / uniroot(function(t) pnhyper(x, 1 / t) - alpha, c(.Machine$double.eps, 1))$root
				} else{
					return 1.0;
				}
			};

			FisherExactFor2x2Result ret;
			//odds
			ret.oddsRatio_ = mle(inPars.TP);

			//p-value and confidence interval
			switch (inPars.pvalTail) {
				case FisherExactFor2x2Input::Tail::TWOTAILED  :
					ret.lowerConfInterval_ = lowerConf(inPars.TP, (1 - inPars.confInterval)/2);
					ret.upperConfInterval_ = upperConf(inPars.TP, (1 - inPars.confInterval)/2);
					ret.pValue_ = pnhyper_twotailed(1);
					break; //optional
				case FisherExactFor2x2Input::Tail::LESSER  :
					ret.lowerConfInterval_ = 0;
					ret.upperConfInterval_ = upperConf(inPars.TP, 1 - inPars.confInterval);
					ret.pValue_ = pnhyper(inPars.TP, 1, false);
					break; //optional
				case FisherExactFor2x2Input::Tail::GREATER  :
					ret.lowerConfInterval_ = lowerConf(inPars.TP, 1 - inPars.confInterval);
					ret.upperConfInterval_ = 0;
					ret.pValue_ = pnhyper(inPars.TP, 1, true);
					break; //optional
			}

			return ret;
		}

		static std::vector<int32_t> genSupportVec(int32_t hi, int32_t lo){
			std::vector<int32_t> support(hi + 1 - lo);
			njh::iota(support, lo);
			return support;
		}
	};

};


}  // namespace njhseq





