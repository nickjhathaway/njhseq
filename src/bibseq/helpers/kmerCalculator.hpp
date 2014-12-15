#pragma once
//
//  kmerCalculator.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 1/4/13.
//  Copyright (c) 2013 Nick Hathaway. All rights reserved.
//

#include "bibseq/objects/helperObjects/kmer.hpp"

namespace bibseq {

class kmerCalculator {

 public:
  kmerCalculator() {}
  static std::unordered_map<std::string, kmer> indexKmer(const std::string& str,
                                                         uint32_t kmerSize);

  static void indexKmer(
      const std::string& str, uint32_t kmerSize,
      std::unordered_map<std::string, kmer>& kTuppleOccurences);

  static void indexKmer(
      const std::string& str, uint32_t kmerSize, const std::string& name,
      uint32_t readCount,
      std::unordered_map<std::string, kmer>& kTuppleOccurences);

  template <class T>
  static std::unordered_map<std::string, kmer> indexKmerForStringMap(
      const std::vector<T>& reads, uint32_t kmerSize);
  template <class T>
  static std::unordered_map<std::string, kmer> indexKmerForStringMap(
      const std::vector<T>& reads,const std::vector<uint64_t> & positons,
      uint32_t kmerSize);

  template <class T>
  static std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      indexKmerForPosMap(const std::vector<T>& reads, uint32_t kmerSize,
      		bool expandPos = false, uint32_t expandSize = 5);
  template <class T>
  static std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      indexKmerForPosMap(const std::vector<T>& reads,const std::vector<uint64_t> & positions,
      		uint32_t kmerSize,
      		bool expandPos = false, uint32_t expandSize = 5);

  // get both kmer type maps
  template <class T>
  static kmerMaps indexKmerMpas(const std::vector<T>& reads,
                                uint32_t kmerLength, uint32_t runCutoff,
                                uint32_t qualRunCutOff,
                            		bool expandPos = false, uint32_t expandSize = 5) {
    double total = 0;
    for (const auto& read : reads) {
      total += read.seqBase_.cnt_;
    }
    return kmerMaps(indexKmerForPosMap(reads, kmerLength, expandPos, expandSize),
                    indexKmerForStringMap(reads, kmerLength), kmerLength,
                    runCutoff, qualRunCutOff, total);
  }

  template <class T>
  static kmerMaps indexKmerMpas(const std::vector<T>& reads,
  															const std::vector<uint64_t> & positions,
                                uint32_t kmerLength, uint32_t runCutoff,
                                uint32_t qualRunCutOff,
                            		bool expandPos = false,
                            		uint32_t expandSize = 5) {
    double total = 0;
    for (const auto& readPos : positions) {
      total += reads[readPos].seqBase_.cnt_;
    }
    return kmerMaps(indexKmerForPosMap(reads,positions, kmerLength, expandPos, expandSize),
                    indexKmerForStringMap(reads, positions, kmerLength), kmerLength,
                    runCutoff, qualRunCutOff, total);
  }
  // simply get all kmers, even repeats
  static VecStr getAllKmers(const std::string& seq, uint32_t kLength);
};

template <class T>
std::unordered_map<std::string, kmer> kmerCalculator::indexKmerForStringMap(
    const std::vector<T>& reads, uint32_t kmerSize) {
  std::unordered_map<std::string, kmer> kTuppleOccurences;
  for_each(reads, [&](const T& read) {
    indexKmer(read.seqBase_.seq_, kmerSize, read.seqBase_.name_,
              read.seqBase_.cnt_, kTuppleOccurences);
  });
  return kTuppleOccurences;
}

template <class T>
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

template <class T>
std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
kmerCalculator::indexKmerForPosMap(const std::vector<T>& reads,
                                   uint32_t kmerSize, bool expandPos,
                                   uint32_t expandSize) {
  std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>>
      outputKmers;
  if(expandPos){
    for (const auto& read : reads) {
      uint32_t cursor = 0;
      while (cursor + kmerSize <= read.seqBase_.seq_.length()) {
        std::string currentKmer = read.seqBase_.seq_.substr(cursor, kmerSize);
        uint32_t startPos = 0;
        uint32_t stopPos = 0;
        if(expandSize >= cursor){
        	startPos = 0;
        }else{
        	startPos = cursor - expandSize;
        }
        if((expandSize + cursor + 1 )> len(read.seqBase_.seq_)){
        	stopPos = len(read.seqBase_.seq_);
        }else{
        	stopPos = cursor + expandSize + 1;
        }
        for(const auto & nCursor : iter::range(startPos, stopPos)){
          if (outputKmers.find(nCursor) == outputKmers.end()) {
            outputKmers.emplace(nCursor, std::unordered_map<std::string, kmer>{
                {currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                   read.seqBase_.cnt_)}});
          } else {
            if (outputKmers[nCursor].find(currentKmer) ==
                outputKmers[nCursor].end()) {
              outputKmers[nCursor].emplace(
                  currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                    read.seqBase_.cnt_));
            } else {
              outputKmers[nCursor][currentKmer]
                  .addPosition(cursor, read.seqBase_.name_, read.seqBase_.cnt_);
            }
          }
        }
        ++cursor;
      }
    }
  }else{
    for (const auto& read : reads) {
      uint32_t cursor = 0;
      while (cursor + kmerSize <= read.seqBase_.seq_.length()) {
        std::string currentKmer = read.seqBase_.seq_.substr(cursor, kmerSize);
        if (outputKmers.find(cursor) == outputKmers.end()) {
          std::unordered_map<std::string, kmer> tempKmerSubMap = {
              {currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                 read.seqBase_.cnt_)}};
          outputKmers.emplace(cursor, tempKmerSubMap);
        } else {
          if (outputKmers[cursor].find(currentKmer) ==
              outputKmers[cursor].end()) {
            outputKmers[cursor].emplace(
                currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                  read.seqBase_.cnt_));
          } else {
            outputKmers[cursor][currentKmer]
                .addPosition(cursor, read.seqBase_.name_, read.seqBase_.cnt_);
          }
        }
        ++cursor;
      }
    }
  }

  return outputKmers;
}


template <class T>
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
      uint32_t cursor = 0;
      while (cursor + kmerSize <= read.seqBase_.seq_.length()) {
        std::string currentKmer = read.seqBase_.seq_.substr(cursor, kmerSize);
        uint32_t startPos = 0;
        uint32_t stopPos = 0;
        if(expandSize >= cursor){
        	startPos = 0;
        }else{
        	startPos = cursor - expandSize;
        }
        if((expandSize + cursor + 1 )> len(read.seqBase_.seq_)){
        	stopPos = len(read.seqBase_.seq_);
        }else{
        	stopPos = cursor + expandSize + 1;
        }
        for(const auto & nCursor : iter::range(startPos, stopPos)){
          if (outputKmers.find(nCursor) == outputKmers.end()) {
            outputKmers.emplace(nCursor, std::unordered_map<std::string, kmer>{
                {currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                   read.seqBase_.cnt_)}});
          } else {
            if (outputKmers[nCursor].find(currentKmer) ==
                outputKmers[nCursor].end()) {
              outputKmers[nCursor].emplace(
                  currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                    read.seqBase_.cnt_));
            } else {
              outputKmers[nCursor][currentKmer]
                  .addPosition(cursor, read.seqBase_.name_, read.seqBase_.cnt_);
            }
          }
        }
        ++cursor;
      }
    }
  }else{
    for (const auto& readPos : positions) {
    	const auto & read = reads[readPos];
      uint32_t cursor = 0;
      while (cursor + kmerSize <= read.seqBase_.seq_.length()) {
        std::string currentKmer = read.seqBase_.seq_.substr(cursor, kmerSize);
        if (outputKmers.find(cursor) == outputKmers.end()) {
          std::unordered_map<std::string, kmer> tempKmerSubMap = {
              {currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                 read.seqBase_.cnt_)}};
          outputKmers.emplace(cursor, tempKmerSubMap);
        } else {
          if (outputKmers[cursor].find(currentKmer) ==
              outputKmers[cursor].end()) {
            outputKmers[cursor].emplace(
                currentKmer, kmer(currentKmer, cursor, read.seqBase_.name_,
                                  read.seqBase_.cnt_));
          } else {
            outputKmers[cursor][currentKmer]
                .addPosition(cursor, read.seqBase_.name_, read.seqBase_.cnt_);
          }
        }
        ++cursor;
      }
    }
  }
  return outputKmers;
}

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "kmerCalculator.cpp"
#endif
