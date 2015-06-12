#pragma once
//
//  mismatch.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/26/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils/stringUtils.hpp"
#include <bibcpp/jsonUtils.h>

namespace bibseq {

class mismatch {

 public:
  mismatch(const char& rBase, uint32_t rQual,
           const std::vector<uint32_t>& rLeadQual,
           const std::vector<uint32_t>& rTrailQual, uint32_t rBasePos,
           const char& sBase, uint32_t sQual,
           const std::vector<uint32_t>& sLeadQual,
           const std::vector<uint32_t>& sTrailQual, uint32_t sBasePos,
           int kFreqByPos, int kFreq)
      : refBase(rBase),
        refQual(rQual),
				refBasePos(rBasePos),
        transition(isMismatchTransition(rBase, sBase)),
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


  std::string outputInfoString() const ;
  Json::Value outputJson()const;

  static bool isMismatchTransition(const char& baseA, const char& baseB);

};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "mismatch.cpp"
#endif
