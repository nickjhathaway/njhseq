#pragma once
//
//  mismatch.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/26/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils/stringUtils.hpp"

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
        refLeadingQual(rLeadQual),
        refTrailingQual(rTrailQual),
        refBasePos(rBasePos),
        seqBase(sBase),
        seqQual(sQual),
        seqLeadingQual(sLeadQual),
        seqTrailingQual(sTrailQual),
        seqBasePos(sBasePos),
        kMerFreqByPos(kFreqByPos),
        kMerFreq(kFreq),
        transition(isMismatchTransition(refBase, seqBase)) {}
  char refBase;
  int refQual;
  std::vector<uint32_t> refLeadingQual;
  std::vector<uint32_t> refTrailingQual;
  uint32_t refBasePos;
  char seqBase;
  int seqQual;
  std::vector<uint32_t> seqLeadingQual;
  std::vector<uint32_t> seqTrailingQual;
  uint32_t seqBasePos;
  int kMerFreqByPos;
  int kMerFreq;
  bool transition;
  /*
  void tempOutputInfo(std::ostream & out){
      out<<refBasePos<<"\t";
      if (transition) {
          out<<"transition\t";
      }else{
          out<<"transversion\t";
      }
      out<<refBase<<"\t"<<refQual<<"\t"<<seqBase<<"\t"<<seqQual;
  }
  std::string tempOutputInfoString(){
      std::stringstream out;
      out<<refBasePos<<"\t";
      if (transition) {
          out<<"transition\t";
      }else{
          out<<"transversion\t";
      }
      out<<refBase<<"\t"<<refQual<<"\t"<<seqBase<<"\t"<<seqQual;
      return out.str();
  }*/
  static bool isMismatchTransition(const char& baseA, const char& baseB);
  std::string outputInfoString() const {
    std::stringstream out;

    if (transition) {
      out << "transition\t";
    } else {
      out << "transversion\t";
    }
    out << refBasePos << "\t" << refBase << "\t" << refQual << "\t"
        << vectorToString(refLeadingQual, ",") << "\t"
        << vectorToString(refTrailingQual, ",") << "\t" << seqBasePos << "\t"
        << seqBase << "\t" << seqQual << "\t"
        << vectorToString(seqLeadingQual, ",") << "\t"
        << vectorToString(seqTrailingQual, ",") << "\t" << kMerFreqByPos << "\t"
        << kMerFreq;
    return out.str();
  }
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "mismatch.cpp"
#endif
