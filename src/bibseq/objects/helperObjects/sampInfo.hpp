#pragma once
//
//  sampInfo.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 05/03/15.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
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
#include "bibseq/common/allSystemIncludes.h"
#include "bibseq/objects/seqObjects/readObject.hpp"

namespace bibseq {
class sampInfo {
public:

  sampInfo();
  sampInfo(const seqInfo& cr);
  sampInfo(const std::string & runName, double totalRunCnt);
  // updates

  /**@brief Update info with this read
   *
   * @param cr The read to update info with
   */
  void update(const seqInfo& cr);

  /**@brief Set runReadCnt to a new value
   *
   * @param runReadCnt The new total count for this sample
   */
  void updateRunReadCnt(double runReadCnt);

  /**@brief set fraction by doing readCnt_/runReadCnt_
   *
   */
  void updateFraction();

  /**@brief reset all numbers but runReadCnt_ to zero
   *
   */
  void resetBasicInfo();

  // samp info
  std::string runName_; /**< The sample name associated with this info*/
  double runReadCnt_; /**< total number of reads associated with runName_(MID rep) in all clusters */
  // sub cluster info
  double readCnt_; /**< amount of reads associated with runName_ (MID rep) */
  double fraction_; /**< fraction of rep reads in this cluster, */
  uint32_t numberOfClusters_;  /**< total number of sub clusters in this cluster */
  // chimeric info
  double chiReadCnt_; /**< The amount of chimeric reads for this cluster*/
  uint32_t chiNumberOfClusters_; /**< The number of chimeric clusters for this cluster*/

  /**@brief Get fraction, total read cnt, and cluster number
   *
   * @param delim The delimiter
   * @return A sting with fraction, totalCnt, and cluster number delimited
   */
  std::string getReadInfo(const std::string& delim = "\t") const;
  /**@brief Get fraction re-calculated with cnt, total read cnt, and cluster number
   *
   * @param cnt The new sample total count
   * @param delim The delimiter
   * @return A sting with fraction, totalCnt, and cluster number delimited
   */
  std::string getReadInfo(uint32_t cnt, const std::string& delim = "\t") const;
  /**@brief Get chimeric fraction, total chimera read cnt, and chimera cluster number
   *
   * @param delim The delimiter
   * @return A sting with chimeric fraction, total chimera read cnt, and chimera cluster number delimited
   */
  std::string getChimeraInfo(const std::string& delim = "\t") const;
  /**@brief Get chimeric fraction re-calculated with cnt, total chimera read cnt, and chimera cluster number
   *
   * @param cnt The new sample total count
   * @param delim The delimiter
   * @return A sting with chimeric fraction, total chimera read cnt, and chimera cluster number delimited
   */
  std::string getChimeraInfo(uint32_t cnt, const std::string& delim = "\t") const;
};


}  // namespace bibseq


