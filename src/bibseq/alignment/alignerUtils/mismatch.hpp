#pragma once
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
//
//  mismatch.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/26/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include <bibcpp/jsonUtils.h>

#include "bibseq/utils/stringUtils.hpp"
#include "bibseq/alignment/alignerUtils/QualScorePars.hpp"


namespace bibseq {

class mismatch {

 public:
  mismatch(const char rBase, uint32_t rQual,
           const std::vector<uint32_t>& rLeadQual,
           const std::vector<uint32_t>& rTrailQual, uint32_t rBasePos,
           const char sBase, uint32_t sQual,
           const std::vector<uint32_t>& sLeadQual,
           const std::vector<uint32_t>& sTrailQual, uint32_t sBasePos,
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
  uint32_t refQual;
  uint32_t refBasePos;
  bool transition;
  uint32_t freq = 1;
  double frac_ = 0.0;
  std::vector<uint32_t> refLeadingQual;
  std::vector<uint32_t> refTrailingQual;

  //specific to a single sequence
  char seqBase;
  uint32_t seqQual;
  std::vector<uint32_t> seqLeadingQual;
  std::vector<uint32_t> seqTrailingQual;
  uint32_t seqBasePos;
  uint32_t kMerFreqByPos;
  uint32_t kMerFreq;

  bool highQuality(const QualScorePars & qScorePars) const;

  std::string outputInfoString() const ;
  Json::Value toJson()const;

  void setTransitionTransverstion();

  static bool isMismatchTransition(const char& baseA, const char& baseB);

};
}  // namespace bibseq


