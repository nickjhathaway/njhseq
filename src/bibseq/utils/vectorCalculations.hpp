#pragma once
//
//  vectorCalculations.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 1/3/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
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
//arma
#include "bibseq/utils/utils.hpp"
#include <vector>
#include <algorithm>
//#include <armadillo>

/// various functions to calculate stats on vectors of any number

namespace bibseq {
template <typename T>
double vectorMedianCopy(std::vector<T> scores) {
  double median = 0.0;
  if (scores.size() != 0) {
    std::size_t size = scores.size();
    std::sort(scores.begin(), scores.end());
    if (size % 2 == 0) {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2.0;
    } else {
      median = scores[size / 2];
    }
  }
  return median;
}

template <typename T>
double vectorMedianRef(std::vector<T>& scores) {
  double median = 0.0;
  if (scores.size() != 0) {
    std::size_t size = scores.size();
    std::sort(scores.begin(), scores.end());
    if (size % 2 == 0) {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2.0;
    } else {
      median = scores[size / 2];
    }
  }
  return median;
}

/*
template <typename T>
double medianTrue(const std::vector<T>& ov){
   std::vector<double> nv(ov.begin(), ov.end());
   return arma::median(arma::vec(nv.data(), nv.size(), false));
}

template<typename T>
T vectorMedian(const std::vector<T>& v) {
	if (v.size() != 0) {
		return arma::median(arma::Col<T>(v));
	}
	return 0;
}*/
/*
template <typename T>
double vectorMedian(std::vector<T> scores) {
  double median = 0.0;
  if (scores.size() != 0) {
    std::size_t size = scores.size();
    std::sort(scores.begin(), scores.end());
    if (size % 2 == 0) {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2.0;
    } else {
      median = scores[size / 2];
    }
  }
  return median;
}*/



template <typename T>
T vectorMinimum(const std::vector<T> & scores) {
  if (scores.size() == 0) {
    return 0;
  }
  auto minEl = std::min_element(scores.begin(), scores.end())	;
  return *minEl;
}

template <typename T>
T vectorMaximum(const std::vector<T> &  scores) {
  if (scores.size() == 0) {
    return 0;
  }
  auto maxEl = std::max_element(scores.begin(), scores.end())	;
  return *maxEl;
}

template <typename T>
T vectorSum(const std::vector<T>& scores) {
	T sum = std::accumulate(scores.begin(), scores.end(), 0.0);
  return sum;
}

template <typename T>
double vectorMean(const std::vector<T>& scores) {
  if (scores.size() != 0) {
  	auto sum = vectorSum(scores);
    return sum / static_cast<double>(scores.size());
  }
  return 0;
}

template <typename T>
double vectorVarianceSamp(const std::vector<T>& scores) {
  if (scores.size() != 0) {
    double meanScore = vectorMean(scores);
    double sumOfSquares = 0.0;
    for (const auto& value : scores) {
      sumOfSquares += std::pow((value - meanScore), 2.0);
    }
    return sumOfSquares / (scores.size() - 1);
  }
  return 0;
}

template <typename T>
double vectorStandardDeviationSamp(const std::vector<T>& scores) {
  if (scores.size() != 0) {
    double vars = vectorVarianceSamp(scores);
    return std::pow(vars, 0.5);
  }
  return 0;
}

template <typename T>
double vectorSEMSamp(const std::vector<T>& scores) {
  if (scores.size() != 0) {
    double vars = vectorVarianceSamp(scores);
    return std::pow(vars, 0.5) / std::sqrt(scores.size());
  }
  return 0;
}

template <typename T>
double vectorVariancePop(const std::vector<T>& scores) {
  if (scores.size() != 0) {
    double meanScore = vectorMean(scores);
    double sumOfSquares = 0.0;
    for (const auto& value : scores) {
      sumOfSquares += std::pow((value - meanScore), 2.0);
    }
    return sumOfSquares / scores.size();
  }
  return 0;
}

template <typename T>
double vectorStandardDeviationPop(const std::vector<T>& scores) {
  if (scores.size() != 0) {
    double vars = vectorVariancePop(scores);
    return std::pow(vars, 0.5);
  }
  return 0;
}

template <typename T>
double vectorSEMPop(const std::vector<T>& scores) {
  if (scores.size() != 0) {
    double vars = vectorVariancePop(scores);
    return std::pow(vars, 0.5) / std::sqrt(scores.size());
  }
  return 0;
}

template <typename T>
std::vector<double> vectorOfZScores(const std::vector<T>& vec, double givenMean,
                                    double givenStd) {
  std::vector<double> ans;
  if (givenMean == 0) {
  	std::stringstream ss;
    ss << "mean can't be zero" << std::endl;
    throw std::runtime_error{ss.str()};
  }
  for (const auto& value : vec) {
    ans.emplace_back((value - givenStd) / givenMean);
  }
  return ans;
}

template <typename T>
std::vector<double> vectorOfZScoresPop(const std::vector<T>& vec) {
  double meanOfScores = vectorMean(vec);
  double stdOfScores = vectorStandardDeviationPop(vec);
  return vectorOfZScores(vec, meanOfScores, stdOfScores);
}
template <typename T>
std::vector<double> vectorOfZScoresSamp(const std::vector<T>& vec) {
  double meanOfScores = vectorMean(vec);
  double stdOfScores = vectorStandardDeviationSamp(vec);
  return vectorOfZScores(vec, meanOfScores, stdOfScores);
}

template<typename T>
double getPearsonCoefficientZScores(const std::vector<T> & zScores1,
		const std::vector<T> & zScores2){
	if(zScores1.size() != zScores2.size()){
		std::stringstream ss;
		ss << "getPearsonCoefficient" << std::endl;
		ss << "scores size must equal" << std::endl;
		throw std::runtime_error{ss.str()};
	}
	return (1/(zScores1.size() -1 )) * std::inner_product(zScores1.begin(), zScores1.end(),
			zScores2.begin(), 0);
}
template<typename T>
double getPearsonCoefficient(const std::vector<T> & scores1,
		const std::vector<T> & scores2){
	return getPearsonCoefficientZScores(vectorOfZScoresSamp(scores1),
			vectorOfZScoresSamp(scores2));
}

template <typename T>
void outputStringVectorMap(const std::map<std::string, std::vector<T>> info,
                           std::ostream& out) {
  // get the mean, range, and median of several vectors
  for (const auto& kv : info) {
    double meanProf = vectorMean(kv.second);
    double maximumProf = vectorMaximum(kv.second);
    double minimumProf = vectorMinimum(kv.second);
    double medianProf = vectorMedian(kv.second);
    out << kv.first << " average: " << meanProf << " median: " << medianProf
        << " range " << minimumProf << ":" << maximumProf << std::endl;
  }
}

template <typename T>
void outputMeanMedainRangeStd(const std::vector<T>& vec, std::ostream& out) {
  double meanProf = mean(vec);
  double maximumProf = maximum(vec);
  double minimumProf = minimum(vec);
  double medianProf = median(vec);
  double stdProf = standardDeviation(vec);
  out << " average: " << meanProf << " median: " << medianProf << " range "
      << minimumProf << ":" << maximumProf << " std: " << stdProf << std::endl;
}



template <typename T>
std::map<std::string, double> getStatsOnVec(const std::vector<T>& vec) {
  return {{"mean", vectorMean(vec)},
          {"median", vectorMedianCopy(vec)},
          {"max", vectorMaximum(vec)},
          {"min", vectorMinimum(vec)},
          {"std", vectorStandardDeviationSamp(vec)},
          {"sum", vectorSum(vec)}};
}
template <typename T>
std::map<std::string, double> getStatsOnVecMore(const std::vector<T>& vec) {
  return {{"mean", vectorMean(vec)},
          {"median", vectorMedianCopy(vec)},
          {"max", vectorMaximum(vec)},
          {"min", vectorMinimum(vec)},
          {"std", vectorStandardDeviationSamp(vec)},
          {"sum", vectorSum(vec)},
          {"sem", vectorSEMSamp(vec)}};
}


template<typename T>
T getSumFromVecStr(const VecStr & strNums){
	auto converted = bib::lexical_cast_con<VecStr, std::vector<T>>(strNums);
	return vectorSum(converted);
}


}  // namespace bib
