#pragma once
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
/*
 * ExtractionStator.hpp
 *
 *  Created on: Sep 29, 2015
 *      Author: nick
 */


#include "bibseq/utils.h"

namespace bibseq {

class ExtractionStator {
public:
	struct extractCounts {
		//good, bad, and contamination should add up to total
		uint32_t good_ = 0;
		//badReverse_, containsNs, minLenBad_, maxLenBad_, qualityFailed_ should add up to bad
		uint32_t bad_ = 0;

		uint32_t badReverse_ = 0;

		uint32_t containsNs_ = 0;
		uint32_t minLenBad_ = 0;
		uint32_t maxLenBad_ = 0;
		uint32_t qualityFailed_ = 0;

		uint32_t contamination_ = 0;

		uint32_t getTotal() const {
			return good_ + bad_ + contamination_;
		}
	};

  enum class extractCase{
  	GOOD,BADREVERSE,CONTAINSNS,MINLENBAD,MAXLENBAD,QUALITYFAILED,CONTAMINATION
  };

  ExtractionStator();
  /**@brief Construct with initial counts
   *
   * @param totalReadCount The total count of the input reads for extraction
   * @param readsUnrecBarcode The number of reads that have unrecognized barcodes
   * @param smallFrags The number of reads that are small fragments
   */
	ExtractionStator(uint32_t totalReadCount, uint32_t readsUnrecBarcode,
			uint32_t smallFrags);

private:
  std::map<std::string,std::map<bool,extractCounts>> counts_; /**< The counts per MID and primer */
  std::map<std::string,std::map<bool,uint32_t>> failedForward_; /**< The counts of failing forward primer for MID */
  uint32_t totalReadCount_ = 0;
  uint32_t readsUnrecBarcode_ = 0;
  uint32_t smallFrags_ = 0;
public:
  /**@brief Increase the count failing the forward primer
   *
   * @param midName The name of the MID
   * @param seqName The name of the sequence to determine direction it was found in
   */
	void increaseFailedForward(const std::string & midName,
			const std::string & seqName);
	/**@brief Increase the count for the failure of a given case
	 *
	 * @param midName The name of the MID and the primer
	 * @param seqName The name of the sequence to determine direction it was found in
	 * @param eCase The case that the read failed
	 */
	void increaseCounts(const std::string & midName, const std::string & seqName,
			extractCase eCase);

	//output

	/**@brief Output stats for the primer plus MID
	 *
	 * @param out the output file
	 * @param delim the delimiter to use for the columns
	 */
	void outStatsPerName(std::ostream & out, const std::string & delim);
	/**@brief Output stats of failing forward primer for a MID
	 *
	 * @param out the output file
	 * @param delim the delimiter to sue for the columns
	 */
	void outFailedForwardStats(std::ostream & out, const std::string & delim);
	/**@brief
	 *
	 * @param out the output file
	 * @param delim the delimiter to sue for the columns
	 */
	void outTotalStats(std::ostream & out, const std::string & delim);

};

}  // namespace bibseq
