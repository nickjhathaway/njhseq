#pragma once
//
//  vectorCalculations.hpp
//
//  Created by Nick Hathaway on 1/3/13.
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
#include "njhseq/utils/utils.hpp"
#include <vector>
#include <algorithm>
#include <boost/math/statistics/univariate_statistics.hpp>

//#include <armadillo>

/// various functions to calculate stats on vectors of any number

namespace njhseq {
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
double vectorSum(const std::vector<T>& scores) {
	double sum = std::accumulate(scores.begin(), scores.end(), 0.0);
  return sum;
}

template <typename T>
T vectorSumSameType(const std::vector<T>& scores) {
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
		ss << __PRETTY_FUNCTION__ << std::endl;
		ss << "scores size must equal" << std::endl;
		throw std::runtime_error{ss.str()};
	}
	return 1/(zScores1.size() -1 ) * std::inner_product(zScores1.begin(), zScores1.end(),
                                                      zScores2.begin(), 0);
}

template<typename T>
double getPearsonCoefficient(const std::vector<T> &scores1,
                             const std::vector<T> &scores2) {
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
	auto converted = njh::lexical_cast_con<VecStr, std::vector<T>>(strNums);
	return vectorSum(converted);
}


template<typename T>
double lins_concordance_correlation(
    const std::vector<T>& x,
    const std::vector<T>& y){
  // If all elements are exactly equal
  bool all_equal = true;
  for (size_t i = 0; i < x.size(); ++i) {
    if (x[i] != y[i]) { all_equal = false; break; }
  }
  if (all_equal) {
    return 1.0;
  }

  if (x.size() != y.size() || x.empty()) {
    throw std::invalid_argument("Vectors must have same non-zero length");
  }
  // Compute means
  const double mean_x = boost::math::statistics::mean(x);
  const double mean_y = boost::math::statistics::mean(y);

  // Compute variances (population)
  const double var_x = boost::math::statistics::variance(x);
  const double var_y = boost::math::statistics::variance(y);

  // Compute covariance
  double cov_xy = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    cov_xy += (x[i] - mean_x) * (y[i] - mean_y);
  }
  cov_xy /= static_cast<double>(x.size());

  // Pearson correlation coefficient
  const double r = cov_xy / std::sqrt(var_x * var_y);

  // Lin’s CCC formula:
  // CCC = (2 * r * σx * σy) / (σx² + σy² + (μx - μy)²)
  const double sd_x = std::sqrt(var_x);
  const double sd_y = std::sqrt(var_y);

  const double numerator = 2.0 * r * sd_x * sd_y;
  const double denominator = var_x + var_y + std::pow(mean_x - mean_y, 2.0);

  return numerator / denominator;
}


class ConcordanceCalculator {
public:
  enum class CccCiMethod { ZTransform, BCa };

  struct CccCiResult {
    double ccc = std::numeric_limits<double>::quiet_NaN();
    double conf_level = 0.95;
    CccCiMethod method = CccCiMethod::ZTransform;

    double lower = std::numeric_limits<double>::quiet_NaN();
    double upper = std::numeric_limits<double>::quiet_NaN();

    // Only meaningful for BCa:
    std::size_t n_boot = 0;
    std::size_t n_used = 0;
};

  // ---- Normal CDF / inverse CDF (Acklam) ----
  static double normal_cdf(const double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2.0));
  }

  static double normal_inv_cdf(const double p) {
    if (!(p > 0.0 && p < 1.0)) {
      throw std::invalid_argument("normal_inv_cdf: p must be in (0,1)");
    }
    constexpr double
        a1 = -3.969683028665376e+01, a2 = 2.209460984245205e+02,
        a3 = -2.759285104469687e+02, a4 = 1.383577518672690e+02,
        a5 = -3.066479806614716e+01, a6 = 2.506628277459239e+00;
    constexpr double
        b1 = -5.447609879822406e+01, b2 = 1.615858368580409e+02,
        b3 = -1.556989798598866e+02, b4 = 6.680131188771972e+01,
        b5 = -1.328068155288572e+01;
    constexpr double
        c1 = -7.784894002430293e-03, c2 = -3.223964580411365e-01,
        c3 = -2.400758277161838e+00, c4 = -2.549732539343734e+00,
        c5 = 4.374664141464968e+00, c6 = 2.938163982698783e+00;
    constexpr double
        d1 = 7.784695709041462e-03, d2 = 3.224671290700398e-01,
        d3 = 2.445134137142996e+00, d4 = 3.754408661907416e+00;

    const double plow = 0.02425, phigh = 1.0 - plow;
    double q;

    if (p < plow) {
      q = std::sqrt(-2.0 * std::log(p));
      return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
             ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    if (p > phigh) {
      q = std::sqrt(-2.0 * std::log(1.0 - p));
      return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
             ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }

    q = p - 0.5;
    const double r = q * q;
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
           (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
  }

  static double quantile_sorted_linear(const std::vector<double> &sorted, const double p) {
    if (sorted.empty()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (p <= 0.0) {
      return sorted.front();
    }
    if (p >= 1.0) {
      return sorted.back();
    }

    const double n = static_cast<double>(sorted.size());
    const double h = 1.0 + (n - 1.0) * p; // R-like type=7 style
    const std::size_t i = static_cast<std::size_t>(std::floor(h));
    const double frac = h - static_cast<double>(i);

    const std::size_t idx0 = (i <= 1) ? 0 : (i - 1);
    const std::size_t idx1 = std::min(idx0 + 1, sorted.size() - 1);

    return sorted[idx0] + frac * (sorted[idx1] - sorted[idx0]);
  }

  // ---- CCC point estimate (population var/cov) ----
  template<typename T>
  static double lins_ccc_point(const std::vector<T> &x, const std::vector<T> &y) {
    static_assert(std::is_arithmetic_v<T>, "T must be numeric");
    if (x.size() != y.size() || x.empty()) {
      throw std::invalid_argument("Vectors must have same non-zero length");
    }
    bool all_equal = true;
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (x[i] != y[i]) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return 1.0;

    double mean_x = 0.0, mean_y = 0.0;
    double M2x = 0.0, M2y = 0.0, C = 0.0;
    double n = 0.0;

    for (std::size_t i = 0; i < x.size(); ++i) {
      const double xi = static_cast<double>(x[i]);
      const double yi = static_cast<double>(y[i]);

      n += 1.0;
      const double dx = xi - mean_x;
      const double dy = yi - mean_y;

      mean_x += dx / n;
      mean_y += dy / n;

      M2x += dx * (xi - mean_x);
      M2y += dy * (yi - mean_y);
      C += dx * (yi - mean_y);
    }

    const double N = static_cast<double>(x.size());
    const double var_x = M2x / N;
    const double var_y = M2y / N;
    const double cov_xy = C / N;

    const double denom = var_x + var_y + std::pow(mean_x - mean_y, 2.0);
    if (denom == 0.0) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return 2.0 * cov_xy / denom;
  }

  // ---- DescTools-matching "z-transform" CI ----
  // Mirrors the CCC() implementation in DescTools (Lin 2000 variance + atanh transform). :contentReference[oaicite:1]{index=1}
  template<typename T>
  static void ccc_ci_ztransform_desctools(
    const std::vector<T> &x,
    const std::vector<T> &y,
    const double conf_level,
    double &est,
    double &lwr,
    double &upr
  ) {
    const std::size_t k = x.size();
    if (k < 3) {
      throw std::invalid_argument("z-transform CI requires n >= 3 (uses k-2).");
    }

    // Means
    double xb = 0.0, yb = 0.0;
    for (std::size_t i = 0; i < k; ++i) {
      xb += static_cast<double>(x[i]);
      yb += static_cast<double>(y[i]);
    }
    xb /= static_cast<double>(k);
    yb /= static_cast<double>(k);

    // Sample variances (ddof=1) and sample covariance -> correlation
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < k; ++i) {
      const double dx = static_cast<double>(x[i]) - xb;
      const double dy = static_cast<double>(y[i]) - yb;
      sxx += dx * dx;
      syy += dy * dy;
      sxy += dx * dy;
    }
    const double varx_sample = sxx / static_cast<double>(k - 1);
    const double vary_sample = syy / static_cast<double>(k - 1);
    const double cov_sample = sxy / static_cast<double>(k - 1);

    const double sd2 = std::sqrt(varx_sample); // sd(dat$x)
    const double sd1 = std::sqrt(vary_sample); // sd(dat$y)

    if (sd1 == 0.0 || sd2 == 0.0) {
      est = lins_ccc_point(x, y);
      lwr = upr = std::numeric_limits<double>::quiet_NaN();
      return;
    }

    const double r = cov_sample / (sd1 * sd2); // cor(dat$x, dat$y)

    // Population variances (DescTools: var()* (k-1)/k)
    const double sx2 = varx_sample * static_cast<double>(k - 1) / static_cast<double>(k);
    const double sy2 = vary_sample * static_cast<double>(k - 1) / static_cast<double>(k);

    const double sxy_pop = r * std::sqrt(sx2 * sy2);
    const double p = 2.0 * sxy_pop / (sx2 + sy2 + std::pow(yb - xb, 2.0)); // CCC
    est = p;

    const double u = (yb - xb) / std::pow(sx2 * sy2, 0.25);

    const double alpha = 1.0 - conf_level;
    const double zv = normal_inv_cdf(1.0 - alpha / 2.0);

    // sep per DescTools (Lin 2000), then transform via atanh (inverse hyperbolic tangent). :contentReference[oaicite:2]{index=2}
    const double sep = std::sqrt(
      (
        (1.0 - r * r) * (p * p) * (1.0 - p * p) / (r * r)
        + (2.0 * std::pow(p, 3.0) * (1.0 - p) * (u * u) / r)
        - (0.5 * std::pow(p, 4.0) * std::pow(u, 4.0) / (r * r))
      ) / static_cast<double>(k - 2)
    );

    const double t = 0.5 * std::log((1.0 + p) / (1.0 - p)); // atanh(p)
    const double set = sep / (1.0 - p * p);

    const double llt = t - zv * set;
    const double ult = t + zv * set;

    // back-transform tanh
    auto tanh_from_atanh = [](const double z) {
      const double e2 = std::exp(2.0 * z);
      return (e2 - 1.0) / (e2 + 1.0);
    };

    lwr = tanh_from_atanh(llt);
    upr = tanh_from_atanh(ult);
  }

  // ---- BCa bits (paired bootstrap + jackknife acceleration) ----
  template<typename T, typename URNG>
  static std::vector<double> bootstrap_ccc_replicates(
    const std::vector<T> &x,
    const std::vector<T> &y,
    const std::size_t n_boot,
    URNG &rng
  ) {
    std::uniform_int_distribution<std::size_t> pick(0, x.size() - 1);
    std::vector<T> xb(x.size()), yb(y.size());
    std::vector<double> boots;
    boots.reserve(n_boot);

    for (std::size_t b = 0; b < n_boot; ++b) {
      for (std::size_t i = 0; i < x.size(); ++i) {
        const std::size_t j = pick(rng);
        xb[i] = x[j];
        yb[i] = y[j];
      }
      if (const double c = lins_ccc_point(xb, yb); std::isfinite(c)) {
        boots.push_back(c);
      }
    }
    return boots;
  }

  template<typename T>
  static double bca_acceleration_jackknife(const std::vector<T> &x, const std::vector<T> &y) {
    const std::size_t n = x.size();
    if (n < 3) {
      return 0.0;
    }
    std::vector<double> jack;
    jack.reserve(n);

    std::vector<T> xj;
    xj.reserve(n - 1);
    std::vector<T> yj;
    yj.reserve(n - 1);

    for (std::size_t leave = 0; leave < n; ++leave) {
      xj.clear();
      yj.clear();
      for (std::size_t i = 0; i < n; ++i) {
        if (i == leave) continue;
        xj.push_back(x[i]);
        yj.push_back(y[i]);
      }
      if (const double t = lins_ccc_point(xj, yj); std::isfinite(t)) jack.push_back(t);
    }
    if (jack.size() < 3) {
      return 0.0;
    }

    const double mean = std::accumulate(jack.begin(), jack.end(), 0.0) / static_cast<double>(jack.size());
    double num = 0.0, den = 0.0;
    for (const double t: jack) {
      const double d = mean - t;
      num += d * d * d;
      den += d * d;
    }
    if (den == 0.0) {
      return 0.0;
    }
    return num / (6.0 * std::pow(den, 1.5));
  }

  // ---- Main CCC function ----
  template<typename T>
  static CccCiResult lins_ccc_with_ci(
    const std::vector<T> &x,
    const std::vector<T> &y,
    double conf_level = 0.95,
    const std::string &ci = "z-transform", // "z-transform" | "bca"
    std::size_t n_boot = 2000, // only for bca
    const std::uint64_t seed = 0 // only for bca
  ) {
    static_assert(std::is_arithmetic_v<T>, "T must be numeric");
    if (x.size() != y.size() || x.empty()) {
      throw std::invalid_argument("Vectors must have same non-zero length");
    }
    if (!(conf_level > 0.0 && conf_level < 1.0)) {
      throw std::invalid_argument("conf_level must be in (0,1)");
    }
    CccCiResult out;
    out.conf_level = conf_level;

    if (ci == "z-transform") {
      out.method = CccCiMethod::ZTransform;
      ccc_ci_ztransform_desctools(x, y, conf_level, out.ccc, out.lower, out.upper);
      return out;
    }

    if (ci == "bca") {
      out.method = CccCiMethod::BCa;
      out.ccc = lins_ccc_point(x, y);
      out.n_boot = n_boot;

      if (n_boot < 100) {
        throw std::invalid_argument("n_boot should be >= 100 (prefer 1000-10000+).");
      }
      std::mt19937_64 rng(seed ? seed : std::random_device{}());
      auto boots = bootstrap_ccc_replicates(x, y, n_boot, rng);
      out.n_used = boots.size();

      if (boots.size() < std::max<std::size_t>(50, n_boot / 10)) {
        throw std::runtime_error("Too many invalid bootstrap replicates (often zero-variance resamples).");
      }
      std::sort(boots.begin(), boots.end());

      const double alpha = 1.0 - conf_level;

      std::size_t count_less = 0;
      for (const double b: boots) {
        if (b < out.ccc) {
          ++count_less;
        }
      }

      double p_less = static_cast<double>(count_less) / static_cast<double>(boots.size());
      const double eps = 1.0 / (2.0 * static_cast<double>(boots.size()));
      p_less = std::min(1.0 - eps, std::max(eps, p_less));

      const double z0 = normal_inv_cdf(p_less);
      const double a = bca_acceleration_jackknife(x, y);

      const double z_lo = normal_inv_cdf(alpha / 2.0);
      const double z_hi = normal_inv_cdf(1.0 - alpha / 2.0);

      auto adj_p = [&z0,&a](const double z) {
        const double num = z0 + z;
        const double den = 1.0 - a * num;
        const double z_adj = z0 + num / den;
        return normal_cdf(z_adj);
      };

      double p1 = std::clamp(adj_p(z_lo), 0.0, 1.0);
      double p2 = std::clamp(adj_p(z_hi), 0.0, 1.0);

      out.lower = quantile_sorted_linear(boots, p1);
      out.upper = quantile_sorted_linear(boots, p2);
      return out;
    }

    throw std::invalid_argument("ci must be \"z-transform\" or \"bca\"");
  }
};


}  // namespace njh
