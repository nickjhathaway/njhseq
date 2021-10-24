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
//  mismatch.hpp
//
//  Created by Nicholas Hathaway on 9/26/13.
//

#include <njhcpp/jsonUtils.h>

#include "njhseq/utils/stringUtils.hpp"
#include "njhseq/alignment/alignerUtils/QualScorePars.hpp"


namespace njhseq {

class mismatch {

 public:
  mismatch(const char rBase, uint8_t rQual,
           const std::vector<uint8_t>& rLeadQual,
           const std::vector<uint8_t>& rTrailQual, uint32_t rBasePos,
           const char sBase, uint8_t sQual,
           const std::vector<uint8_t>& sLeadQual,
           const std::vector<uint8_t>& sTrailQual, uint32_t sBasePos,
           int kFreqByPos, int kFreq)
      : refBase(rBase),
        refQual(rQual),
				refBasePos(rBasePos),
        transition(false),
        refLeadingQual(rLeadQual),
        refTrailingQual(rTrailQual),
        seqBase(sBase),
        seqQual(sQual),
        seqLeadingQual(sLeadQual),
        seqTrailingQual(sTrailQual),
        seqBasePos(sBasePos),
        kMerFreqByPos(kFreqByPos),
        kMerFreq(kFreq){}
  //mismatch info shared between all seqs sharing this mismatch
  char refBase;
  uint8_t refQual;
  uint32_t refBasePos;
  bool transition;
  uint32_t freq = 1;
  double frac_ = 0.0;
  std::vector<uint8_t> refLeadingQual;
  std::vector<uint8_t> refTrailingQual;

  //specific to a single sequence
  char seqBase;
  uint8_t seqQual;
  std::vector<uint8_t> seqLeadingQual;
  std::vector<uint8_t> seqTrailingQual;
  uint32_t seqBasePos;
  uint32_t kMerFreqByPos;
  uint32_t kMerFreq;

  bool highQuality(const QualScorePars & qScorePars) const;
  bool highQualityJustRef(const QualScorePars & qScorePars) const;
  bool highQualityJustSeq(const QualScorePars & qScorePars) const;

  std::string outputInfoString() const ;
  Json::Value toJson()const;

  void setTransitionTransverstion();

  static bool isMismatchTransition(const char& baseA, const char& baseB);

};
}  // namespace njhseq


