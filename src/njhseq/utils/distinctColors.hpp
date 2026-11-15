#pragma once

/*
 * distinctColors.hpp
 *
 *  Created on: Nov 14, 2026
 *      Author: nick
 */

#include "njhseq/common.h"
#include <njhcpp/graphics/color.hpp>

namespace njhseq {

/**@brief Settings for genDistinctColors(), the region of CIE Lab colors to pick from, the defaults are iwanthue's default preset
 *
 */
struct DistinctColorsPars {
	double lightnessMin_ = 35;
	double lightnessMax_ = 80;
	double chromaMin_ = 30;
	double chromaMax_ = 80;
	double gridStep_ = 5; /**< spacing of the candidate colors in L, a and b*/
	uint32_t maxIterations_ = 50; /**< the most k-means rounds*/
};

/**@brief Colors as far apart as possible, in the style of iwanthue (https://medialab.github.io/iwanthue/)
 *
 * Candidates are a grid of CIE Lab colors (D65) that are in sRGB and within the lightness and chroma ranges,
 * the starting colors are picked one at a time as the candidate farthest from those already picked (KMeans::farthestPoints),
 * then refined with k-means over the candidates (KMeans::run, centers snapped to unused candidates), and ordered so each color is the farthest from the ones before it.
 * Every tie goes to the lowest index so the colors are always the same, seqServer's utils.js distinctColorHexStrs() has the same algorithm so the viewers and njhseq agree.
 * With more colors asked for than there are candidates, the candidates repeat.
 *
 * @param count the number of colors
 * @param pars the region of colors to use
 * @return hex strings, e.g. #a1b2c3
 */
[[nodiscard]] VecStr genDistinctColorHexStrs(uint32_t count, const DistinctColorsPars & pars = DistinctColorsPars{});

/**@brief genDistinctColorHexStrs() as colors
 *
 */
[[nodiscard]] std::vector<njh::color> genDistinctColors(uint32_t count, const DistinctColorsPars & pars = DistinctColorsPars{});

}  // namespace njhseq
