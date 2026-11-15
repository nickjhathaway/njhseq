/*
 * distinctColors.cpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "distinctColors.hpp"
#include "kMeans.hpp"

namespace njhseq {

namespace {

double labFInv(double t) {
	const double cubed = t * t * t;
	return cubed > 0.008856 ? cubed : (t - 16.0 / 116.0) / 7.787;
}

double srgbGamma(double c) {
	return c <= 0.0031308 ? 12.92 * c : 1.055 * std::pow(c, 1.0 / 2.4) - 0.055;
}

//sRGB 0-1 from Lab (L, a, b) with the D65 white point, false when it's outside sRGB
bool labToRgb(const KMeans::Point & lab, double & red, double & green, double & blue) {
	const double fy = (lab[0] + 16.0) / 116.0;
	const double fx = fy + lab[1] / 500.0;
	const double fz = fy - lab[2] / 200.0;
	const double x = 0.95047 * labFInv(fx);
	const double y = 1.0 * labFInv(fy);
	const double z = 1.08883 * labFInv(fz);
	red = srgbGamma(3.2404542 * x - 1.5371385 * y - 0.4985314 * z);
	green = srgbGamma(-0.9692660 * x + 1.8760108 * y + 0.0415560 * z);
	blue = srgbGamma(0.0556434 * x - 0.2040259 * y + 1.0572252 * z);
	return red >= 0 && red <= 1 && green >= 0 && green <= 1 && blue >= 0 && blue <= 1;
}

std::string toHex(const KMeans::Point & lab) {
	double red = 0;
	double green = 0;
	double blue = 0;
	labToRgb(lab, red, green, blue);
	static const char digits[] = "0123456789abcdef";
	std::string ret = "#";
	for (const auto channel : {red, green, blue}) {
		const auto val = static_cast<uint32_t>(std::round(channel * 255));
		ret.push_back(digits[val / 16]);
		ret.push_back(digits[val % 16]);
	}
	return ret;
}

}  // namespace

VecStr genDistinctColorHexStrs(uint32_t count, const DistinctColorsPars & pars) {
	if (0 == count) {
		return {};
	}
	//candidates, a grid in Lab within sRGB and the lightness and chroma ranges
	std::vector<KMeans::Point> candidates;
	for (uint32_t lPos = 0; pars.lightnessMin_ + lPos * pars.gridStep_ <= pars.lightnessMax_; ++lPos) {
		for (uint32_t aPos = 0; -100 + aPos * pars.gridStep_ <= 100; ++aPos) {
			for (uint32_t bPos = 0; -100 + bPos * pars.gridStep_ <= 100; ++bPos) {
				const KMeans::Point lab{pars.lightnessMin_ + lPos * pars.gridStep_, -100 + aPos * pars.gridStep_, -100 + bPos * pars.gridStep_};
				const double chroma = std::sqrt(lab[1] * lab[1] + lab[2] * lab[2]);
				double red = 0;
				double green = 0;
				double blue = 0;
				if (chroma >= pars.chromaMin_ && chroma <= pars.chromaMax_ && labToRgb(lab, red, green, blue)) {
					candidates.emplace_back(lab);
				}
			}
		}
	}
	if (candidates.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no colors are in sRGB within lightness " << pars.lightnessMin_ << "-" << pars.lightnessMax_
				<< " and chroma " << pars.chromaMin_ << "-" << pars.chromaMax_ << "\n";
		throw std::runtime_error{ss.str()};
	}
	const auto candidatesMean = KMeans::mean(candidates);
	//spread out starting colors, refined by k-means with each center kept on a candidate
	KMeans::Pars kPars;
	kPars.maxIterations_ = pars.maxIterations_;
	kPars.snapToPoints_ = true;
	const auto starts = KMeans::farthestPoints(candidates, std::min<uint32_t>(count, candidates.size()), candidatesMean);
	const auto clusters = KMeans::run(candidates, starts, kPars);
	//ordered so each color is the farthest from the ones before it, starting with the one farthest from the mean
	const auto order = KMeans::farthestPoints(clusters.centers_, clusters.centers_.size(), candidatesMean);
	VecStr ret;
	ret.reserve(count);
	for (uint32_t pos = 0; pos < count; ++pos) {
		ret.emplace_back(toHex(clusters.centers_[order[pos % order.size()]]));
	}
	return ret;
}

std::vector<njh::color> genDistinctColors(uint32_t count, const DistinctColorsPars & pars) {
	std::vector<njh::color> ret;
	for (const auto & hex : genDistinctColorHexStrs(count, pars)) {
		ret.emplace_back(hex);
	}
	return ret;
}

}  // namespace njhseq
