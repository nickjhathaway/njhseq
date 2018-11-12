#pragma once
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
//
//  alignerUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/23/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//


#include <njhcpp/jsonUtils.h>
#include "njhseq/utils.h"
#include "njhseq/alignment/alignerUtils/mismatch.hpp"
#include "njhseq/alignment/alignerUtils/gaps.hpp"
#include "njhseq/alignment/alignerUtils/QualScorePars.hpp"

namespace njhseq {

/**@brief Information for one strain in a distance comparison between another one
 *
 */
class DistanceMet {
public:

	uint32_t identities_ = 0; /**< number of bases matched out of total*/
	double identity_ = 0;     /**< percentage of bases matched out of total*/
	uint32_t covered_ = 0;    /**< number of bases covered including gaps, matches, mismatches*/
	double coverage_ = 0;     /**< percentage of bases overlapped */
	/**@brief reset all numbers to zero
	 *
	 */
	void reset();
  /**@brief convert to json representation
   *
   * @return Json::Value object
   */
	Json::Value toJson()const;
};

/**@brief The distance comparison between two strains/variants
 *
 */
class DistanceComp {
public:

	//the following three percents are of the overlap of the alignment only
	uint32_t basesInAln_ = 0; /**< number of bases in the overlap alignment including gaps */
	double percentMismatch_ = 0; /**< percent of the overlap alignment that is mismatches */
	double percentMatch_ = 0; /**< percent of the overlap alignment that is matches */
	double percentGaps_ = 0; /**< percent of the overlap alignment that is gap (not counting gaps but gap bases) */
	uint32_t overLappingEvents_ = 0; /**< The total number of events that could happen */
	double eventBasedIdentity_ = 0; /**< matches over number of overlapping events */
	uint32_t overLappingEventsHq_ = 0; /**< The total number of events that could happen */
	double eventBasedIdentityHq_ = 0; /**< matches over number of overlapping events */
	DistanceMet ref_; //alignObj A
	DistanceMet query_; //alignObj B

	std::map<uint32_t, mismatch> mismatches_; /**< The mismatches in the alignment, key is alignment position*/
	std::map<uint32_t, mismatch> lowKmerMismatches_;/**< The mismatches in the alignment that were classified as low kmer frequency, key is alignment position*/
	std::map<uint32_t, gap> alignmentGaps_;/**< The gaps in the alignment, key is alignment position*/

	/**@brief reset all numbers to zero
	 *
	 */
	void reset();

	/**@brief Get number of events, the sum of gaps and mismatches
	 *
	 * @param countLowKmer Whether to include the low kmer mismatches in the event count
	 * @return The number of gaps and mismatches
	 */
	uint32_t getNumOfEvents(bool countLowKmer) const;



  /**@brief convert to json representation
   *
   * @return Json::Value object
   */
	Json::Value toJson() const;
};


/**@brief The comparison between two different variants
 *
 */
class comparison {
public:

  double oneBaseIndel_ = 0; /**< Number of one base indels, can be weight for indels in homopolymers*/
  double twoBaseIndel_ = 0; /**< Number of two base indels, can be weight for indels in homopolymers*/
  double largeBaseIndel_ = 0; /**< Number of larger than two base indels, can be weight for indels in homopolymers*/
  uint32_t hqMismatches_ = 0; /**< Number of high quality mismatches*/
  uint32_t lqMismatches_ = 0; /**< Number of low quality mismatches*/
  uint32_t lowKmerMismatches_ = 0; /**< Number of low kmer frequency mismatches*/

  uint32_t highQualityMatches_ = 0; /**< Number of high quality matches*/
  uint32_t lowQualityMatches_ = 0; /**< Number of low quality matches*/

  double alnScore_ = 0; /**< alnignment score*/

  DistanceComp distances_; /**< The distance metrics in the comparison of the two variants*/

  std::string refName_; /**< The name of the variant that goes first the alignment (alignObjectA) which is normally a reference/major variant*/
  std::string queryName_; /**< The name of the variant that goes goes the alignment (alignObjectB) which is normally a query/minor variant*/

  void recalcMismatchQuality(const QualScorePars & pars);

	void setEventBaseIdentity();
	void setEventBaseIdentityHq();

	/**@brief reset all numbers to zero
	 *
	 */
  void resetCounts();

  /**@brief Whether the generated errors pass the threshold held in this comparison
   *
   * @param generatedError The generated error profile in an alignment comparison
   * @return Whether the generatedErrors fall below the errors in this comparison object
   */
	bool passErrorProfile(const comparison& generatedError) const;

	/**@brief Whether the generated errors pass the event based identity threshold held in this comparison
	 *
	 * @param generatedError The generated error profile in an alignment comparison
	 * @return Whether the generatedErrors falls above event based identity threshold in this comparison object
	 */
	bool passIdThreshold(const comparison& generatedError) const;

	/**@brief Whether the generated errors pass the event based identity threshold (calculated by using only hq errors) held in this comparison
	 *
	 * @param generatedError The generated error profile in an alignment comparison
	 * @return Whether the generatedErrors falls above event based identity threshold in this comparison object
	 */
	bool passIdThresholdHq(const comparison& generatedError) const;

	/**@brief Whether the generated errors pass the event based identity threshold and the errors held in this comparison
	 *
	 * @param generatedError The generated error profile in an alignment comparison
	 * @return Whether the generatedErrors falls above event based identity threshold and below the errors in this comparison object
	 */
	bool passIdAndErrorThreshold(const comparison& generatedError) const;

  /**@brief convert to json representation
   *
   * @return Json::Value object
   */
  Json::Value toJson() const;

};

std::ostream & operator <<(std::ostream & out, const comparison & comp);



}  // namespace njhseq


