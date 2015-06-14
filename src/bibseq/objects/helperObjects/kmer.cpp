//
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
#include "kmer.hpp"

namespace bibseq {

void kmer::addPosition(uint32_t pos) {
  positions_.push_back(pos);
  ++count_;
}

void kmer::addPosition(uint32_t pos, const std::string& name,
                       uint32_t numReads) {
  positions_.push_back(pos);
  count_ += numReads;
  auto it = names_.find(name);
  // check to see if it occured more than once on this read
  if (it == names_.end()) {
    readCnt_ += numReads;
    names_[name] = 1;
  } else {
    ++(it->second);
  }
}
void kmer::infoLine(std::ostream& out, bool header) const{
	if(header){
		out << "k\tcount\tpositions" << std::endl;
	}else{
		out << k_ << "\t" << count_ << "\t"
		      << vectorToString(positions_, ",") << std::endl;
	}
}
void kmer::multipleInfoLine(const std::vector<kmer> & kmers,
		std::ostream& out) {
	kmers.front().infoLine(out, true);
	for_each(kmers, [&](const kmer & currentK){ currentK.infoLine(out);});
}
void kmer::multipleInfoLine(const std::unordered_map<std::string, kmer> & kmers,
		std::ostream& out) {
	kmers.begin()->second.infoLine(out, true);
	for_each(kmers, [&](decltype(*kmers.begin()) & currentK){ currentK.second.infoLine(out);});
}

void kmer::outputInfo(std::ostream& out) const {
  out << "k: " << k_ << " count: " << count_ << " positions_"
      << vectorToString(positions_, ",") << std::endl;
}

// comparisons
bool kmer::operator<(const kmer& otherKmer) const {
  if (positions_[0] < otherKmer.positions_[0]) {
    return true;
  } else {
    return false;
  }
}

bool kmer::operator>(const kmer& otherKmer) const {
  if (positions_[0] > otherKmer.positions_[0]) {
    return true;
  } else {
    return false;
  }
}

bool kmer::operator==(const kmer& otherKmer) const {
  if (positions_[0] == otherKmer.positions_[0]) {
    return true;
  } else {
    return false;
  }
}
bool kmer::operator>=(const kmer& otherKmer) const {
  if (positions_[0] >= otherKmer.positions_[0]) {
    return true;
  } else {
    return false;
  }
}

bool kmer::operator<=(const kmer& otherKmer) const {
  if (positions_[0] <= otherKmer.positions_[0]) {
    return true;
  } else {
    return false;
  }
}

bool kmer::operator==(const std::string& kString) const {
  return (k_ == kString);
}

bool kmerMaps::isKmerLowFrequency(const std::string& kmer, uint32_t position,
                                  bool kmersByPositions, uint32_t cutOff) {
  if (kmersByPositions) {
    return (kmersByPos_[position].find(kmer) == kmersByPos_[position].end() ||
            kmersByPos_[position][kmer].readCnt_ <= cutOff);
  } else {
    return (kmersAnyWhere_.find(kmer) == kmersAnyWhere_.end() ||
            kmersAnyWhere_[kmer].readCnt_ <= cutOff);
  }
}

void kmerMaps::outputKmerInfo(kmerMaps kMaps, std::ostream& out) {
  out << "pos\tkmer\treadCount\treadCountAnywhere" << std::endl;
  for (const auto& kIter : kMaps.kmersByPos_) {
    for (const auto& kPos : kIter.second) {
      out << kIter.first << "\t" << kPos.second.k_ << "\t"
          << kPos.second.readCnt_ << "\t"
          << kMaps.kmersAnyWhere_[kPos.second.k_].readCnt_ << std::endl;
    }
  }
}

}  // namespace bib
