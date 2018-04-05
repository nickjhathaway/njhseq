#pragma once
//
//  kmerCalculator.hpp
//
//  Created by Nick Hathaway on 1/4/13.
//
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bibseq/objects/kmer/kmer.hpp"
#include "bibseq/objects/seqObjects/Paired/PairedCluster.hpp"

namespace bibseq {

class kmerCalculator {

 public:
  static std::unordered_map<std::string, kmer> indexKmer(const std::string& str,
                                                         uint32_t kmerSize);

  static void indexKmer(
      const std::string& str, uint32_t kmerSize,
      std::unordered_map<std::string, kmer>& kTuppleOccurences);

  static void indexKmer(
      const std::string& str, uint32_t kmerSize, const std::string& name,
      uint32_t readCount,
      std::unordered_map<std::string, kmer>& kTuppleOccurences);

  template <typename T>
  static std::unordered_map<std::string, kmer> indexKmerForStringMap(
      const std::vector<T>& reads, uint32_t kmerSize);
  template <typename T>
  static std::unordered_map<std::string, kmer> indexKmerForStringMap(
      const std::vector<T>& reads,const std::vector<uint64_t> & positons,
      uint32_t kmerSize);
/*
  template <>
  static std::unordered_map<std::string, kmer> indexKmerForStringMap(
      const std::vector<PairedCluster>& reads, uint32_t kmerSize);
  template <>
  static std::unordered_map<std::string, kmer> indexKmerForStringMap(
      const std::vector<PairedCluster>& reads,const std::vector<uint64_t> & positons,
      uint32_t kmerSize);*/

  template <typename T>
  static std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      indexKmerForPosMap(const std::vector<T>& reads, uint32_t kmerSize,
      		bool expandPos = false, uint32_t expandSize = 5);
  template <typename T>
  static std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      indexKmerForPosMap(const std::vector<T>& reads,const std::vector<uint64_t> & positions,
      		uint32_t kmerSize,
      		bool expandPos = false, uint32_t expandSize = 5);
  /*
  template <>
  static std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      indexKmerForPosMap(const std::vector<PairedCluster>& reads, uint32_t kmerSize,
      		bool expandPos = false, uint32_t expandSize = 5);
  template <>
  static std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      indexKmerForPosMap(const std::vector<PairedCluster>& reads,const std::vector<uint64_t> & positions,
      		uint32_t kmerSize,
      		bool expandPos = false, uint32_t expandSize = 5);*/


  static void addKmerForPos(
  		std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & outputKmers,
  		uint64_t pos, const std::string currentKmer, const seqInfo & sInfo) ;

  static std::pair<uint32_t, uint32_t> determineExpandPos(uint32_t expandSize, uint32_t pos, const seqInfo & sInfo){
    uint32_t startPos = 0;
    uint32_t stopPos = 0;
    if(expandSize >= pos){
    	startPos = 0;
    }else{
    	startPos = pos - expandSize;
    }
    if((expandSize + pos + 1 )> sInfo.seq_.size()){
    	stopPos = sInfo.seq_.size();
    }else{
    	stopPos = pos + expandSize + 1;
    }
    return {startPos, stopPos};
  }

  static void addStrKmersPos(
  		std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & outputKmers,
  		const seqInfo & sInfo, uint32_t kmerLen);

  static void addStrKmersPosExpand(
  		std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & outputKmers,
  		const seqInfo & sInfo, uint32_t kmerLen, uint32_t expandSize);



  // simply get all kmers, even repeats
  static VecStr getAllKmers(const std::string& seq, uint32_t kLength);
};

template <typename T>
std::unordered_map<std::string, kmer> kmerCalculator::indexKmerForStringMap(
    const std::vector<T>& reads, uint32_t kmerSize) {
	std::vector<uint64_t> positions(reads.size());
	bib::iota<uint64_t>(positions, 0);
	return indexKmerForStringMap(reads, positions, kmerSize);
}

template <typename T>
std::unordered_map<std::string, kmer> kmerCalculator::indexKmerForStringMap(
    const std::vector<T>& reads,const std::vector<uint64_t> & positons,
    uint32_t kmerSize) {
  std::unordered_map<std::string, kmer> kTuppleOccurences;
  for(const auto & readPos : positons){
  	const auto & read = reads[readPos];
  	indexKmer(read.seqBase_.seq_, kmerSize, read.seqBase_.name_,
  	              read.seqBase_.cnt_, kTuppleOccurences);
  }
  return kTuppleOccurences;
}

template <typename T>
std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
kmerCalculator::indexKmerForPosMap(const std::vector<T>& reads,
                                   uint32_t kmerSize, bool expandPos,
                                   uint32_t expandSize) {
	std::vector<uint64_t> positions(reads.size());
	bib::iota<uint64_t>(positions, 0);
	return indexKmerForPosMap(reads, positions, kmerSize,expandPos, expandSize);
}




template <typename T>
std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
kmerCalculator::indexKmerForPosMap(const std::vector<T>& reads,
																	 const std::vector<uint64_t> & positions,
                                   uint32_t kmerSize, bool expandPos,
                                   uint32_t expandSize) {
  std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      outputKmers;
  if(expandPos){
    for (const auto& readPos : positions) {
    	const auto & read = reads[readPos];
    	addStrKmersPosExpand(outputKmers, read.seqBase_, kmerSize, expandSize);
    }
  }else{
    for (const auto& readPos : positions) {
    	const auto & read = reads[readPos];
    	addStrKmersPos(outputKmers, read.seqBase_, kmerSize);
    }
  }
  return outputKmers;
}


template <>
inline std::unordered_map<std::string, kmer> kmerCalculator::indexKmerForStringMap(
    const std::vector<PairedCluster>& reads,const std::vector<uint64_t> & positons,
    uint32_t kmerSize){
  std::unordered_map<std::string, kmer> kTuppleOccurences;
  for(const auto & readPos : positons){
  	const auto & read = reads[readPos];
  	indexKmer(read.seqBase_.seq_, kmerSize, read.seqBase_.name_,
  	              read.seqBase_.cnt_, kTuppleOccurences);
  	indexKmer(read.mateSeqBase_.seq_, kmerSize, read.seqBase_.name_,
  	  	              read.seqBase_.cnt_, kTuppleOccurences);
  }
  return kTuppleOccurences;
}

template <>
inline std::unordered_map<std::string, kmer> kmerCalculator::indexKmerForStringMap(
    const std::vector<PairedCluster>& reads, uint32_t kmerSize){
	std::vector<uint64_t> positions(reads.size());
	bib::iota<uint64_t>(positions, 0);
	return indexKmerForStringMap(reads, positions, kmerSize);
}


template <>
inline std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
kmerCalculator::indexKmerForPosMap(const std::vector<PairedCluster>& reads,const std::vector<uint64_t> & positions,
    		uint32_t kmerSize,
    		bool expandPos , uint32_t expandSize ){
  std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      outputKmers;
  if(expandPos){
    for (const auto& readPos : positions) {
    	const auto & read = reads[readPos];
    	addStrKmersPosExpand(outputKmers, read.seqBase_, kmerSize, expandSize);
    	addStrKmersPosExpand(outputKmers, read.mateSeqBase_, kmerSize, expandSize);
    }
  }else{
    for (const auto& readPos : positions) {
    	const auto & read = reads[readPos];
    	addStrKmersPos(outputKmers, read.seqBase_, kmerSize);
    	addStrKmersPos(outputKmers, read.mateSeqBase_, kmerSize);
    }
  }
  return outputKmers;
}

template <>
inline std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
kmerCalculator::indexKmerForPosMap(const std::vector<PairedCluster>& reads, uint32_t kmerSize,
    		bool expandPos, uint32_t expandSize ){
	std::vector<uint64_t> positions(reads.size());
	bib::iota<uint64_t>(positions, 0);
	return indexKmerForPosMap(reads, positions, kmerSize,expandPos, expandSize);
}


}  // namespace bibseq


