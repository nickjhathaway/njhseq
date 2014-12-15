#pragma once
//
//  seqInfo.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 03/03/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/alignment/substituteMatrix.hpp"
#include <bibcpp/jsonUtils.h>

namespace bibseq {

struct seqInfo {
  // constructors
  seqInfo() : name_(""), seq_(""), cnt_(1), frac_(0) {}
  seqInfo(const std::string& name, const std::string& seq,
          const std::vector<uint32_t>& qual)
      : name_(name), seq_(seq), qual_(qual), cnt_(1), frac_(0) {}
  seqInfo(const std::string& name, const std::string& seq,
          const std::vector<uint32_t>& qual, double cnt)
      : name_(name), seq_(seq), qual_(qual), cnt_(cnt), frac_(0) {}
  seqInfo(const std::string& name, const std::string& seq)
      : name_(name),
        seq_(seq),
        qual_(std::vector<uint32_t>(seq.size(), 40)),
        cnt_(1),
        frac_(0) {}

  seqInfo(const std::string& name, const std::string& seq,
          const std::string& stringQual)
      : name_(name),
        seq_(seq),
        qual_(stringToVector<uint32_t>(stringQual)),
        cnt_(1),
        frac_(0) {}

  seqInfo(const std::string& name, const std::string& seq,
          const std::string& stringQual, int off_set)
      : name_(name),
        seq_(seq),
        qual_(std::vector<uint32_t>(0)),
        cnt_(1),
        frac_(0) {
    for (size_t i = 0; i < stringQual.size(); ++i) {
      qual_.push_back(stringQual[i] - off_set);
    }
  }
  seqInfo(const std::string& name, const std::string& seq,
          const std::vector<uint32_t>& qual, double cnt, double frac)
      : name_(name), seq_(seq), qual_(qual), cnt_(cnt), frac_(frac) {}
  // Members
  std::string name_;
  std::string seq_;
  std::vector<uint32_t> qual_;
  double cnt_;
  double frac_;

  // functions
  // get a sub-portion of the read
  seqInfo getSubRead(uint32_t pos, uint32_t size) const;
  seqInfo getSubRead(uint32_t pos) const;
  // changing the seq and qual of the read
  void prepend(const std::string& seq, const std::vector<uint32_t>& qual);
  void append(const std::string& seq, const std::vector<uint32_t>& qual);
  void prepend(const std::string& seq, uint32_t defaultQuality = 40);
  void append(const std::string& seq, uint32_t defaultQuality = 40);
  bool checkLeadQual(int pos, uint32_t secondayQual, int out = 5) const;
  bool checkTrailQual(int pos, uint32_t secondayQual, int out = 5) const;
  bool checkPrimaryQual(int pos, uint32_t primaryQual) const;
  bool checkQual(int pos, uint32_t primaryQual, uint32_t secondayQual,
                 int out = 5) const;
  int findLowestNeighborhoodQual(int posA, int out = 5) const;
  const std::vector<uint32_t> getLeadQual(int posA, int out = 5) const;
  const std::vector<uint32_t> getTrailQual(int posA, int out = 5) const;
  // get a quality string from vector
  std::string getQualString() const;
  std::string getFastqQualString(int offset) const;
  static std::string getFastqString(const std::vector<uint32_t>& quals,
                                    uint32_t offset);
  // remove base at position
  void removeBase(size_t pos);
  void removeLowQualityBases(uint32_t qualCutOff);
  // handle gaps
  void removeGaps();
  // change name
  void markAsChimeric();
  //reverse complement read
  void reverseComplementRead(bool mark = false);

  //comparison
  bool degenCompare(const seqInfo & otherInfo,
  		const substituteMatrix & compareScores)const;
  // output
  void outPutFastq(std::ostream& fastqFile) const;
  void outPutSeq(std::ostream& fastaFile) const;
  void outPutSeqAnsi(std::ostream& fastaFile) const;
  void outPutQual(std::ostream& qualFile) const;
  // description
  void printDescription(std::ostream& out, bool deep) const;
  Json::Value toJson()const;
  const static std::unordered_map<char, uint32_t> ansiBaseColor;

};

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "seqInfo.cpp"
#endif
