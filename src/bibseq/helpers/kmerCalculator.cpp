#include "kmerCalculator.hpp"

namespace bibseq {
void kmerCalculator::addStrKmersPos(
		std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & outputKmers,
		const seqInfo & sInfo, uint32_t kmerLen){
  uint32_t cursor = 0;
  while (cursor + kmerLen <= sInfo.seq_.length()) {
    std::string currentKmer = sInfo.seq_.substr(cursor, kmerLen);
    addKmerForPos(outputKmers,cursor, currentKmer, sInfo);
    ++cursor;
  }
}

void kmerCalculator::addStrKmersPosExpand(
		std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & outputKmers,
		const seqInfo & sInfo, uint32_t kmerLen, uint32_t expandSize){
  uint32_t cursor = 0;
  while (cursor + kmerLen <= sInfo.seq_.length()) {
    std::string currentKmer = sInfo.seq_.substr(cursor, kmerLen);
    auto kPositions = determineExpandPos(expandSize, cursor, sInfo);
    for(const auto & nCursor : iter::range(kPositions.first, kPositions.second)){
    	addKmerForPos(outputKmers,nCursor, currentKmer, sInfo);
    }
    ++cursor;
  }
}
void kmerCalculator::addKmerForPos(
		std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & outputKmers,
		uint64_t pos, const std::string currentKmer, const seqInfo & sInfo) {
  auto search = outputKmers.find(pos);
  if (search == outputKmers.end()) {
    outputKmers.emplace(pos, std::unordered_map<std::string, kmer>{
        {currentKmer, kmer(currentKmer, pos, sInfo.name_,
                           sInfo.cnt_)}});
  } else {
  	auto subSearch = search->second.find(currentKmer);
    if (subSearch ==  search->second.end()) {
    	search->second.emplace(
          currentKmer, kmer(currentKmer, pos, sInfo.name_,
                            sInfo.cnt_));
    } else {
    	subSearch->second.addPosition(pos, sInfo.name_, sInfo.cnt_);
    }
  }
}

std::unordered_map<std::string, kmer> kmerCalculator::indexKmer(
    const std::string& str, uint32_t kmerSize) {
  std::unordered_map<std::string, kmer> kTuppleOccurences;
  indexKmer(str, kmerSize, kTuppleOccurences);
  return kTuppleOccurences;
}

void kmerCalculator::indexKmer(
    const std::string& str, uint32_t kmerSize,
    std::unordered_map<std::string, kmer>& kTuppleOccurences) {
  uint32_t cursor = 0;
  while (cursor + kmerSize <= static_cast<uint32_t>(str.size())) {
  	auto currentKmer = str.substr(cursor, kmerSize);
  	auto search = kTuppleOccurences.find(currentKmer);
    if ( search == kTuppleOccurences.end()) {
      kTuppleOccurences.emplace(currentKmer,
                                kmer(currentKmer, cursor));
    } else {
    	search->second.addPosition(cursor);
    }
    ++cursor;
  }
}
void kmerCalculator::indexKmer(
    const std::string& str, uint32_t kmerSize, const std::string& name,
    uint32_t readCount,
    std::unordered_map<std::string, kmer>& kTuppleOccurences) {
  uint32_t cursor = 0;

  while (cursor + kmerSize <= static_cast<uint32_t>(str.size())) {
  	auto currentKmer = str.substr(cursor, kmerSize);
  	auto search = kTuppleOccurences.find(currentKmer);
    if (search == kTuppleOccurences.end()) {
      kTuppleOccurences.emplace(
      		currentKmer,
          kmer(currentKmer, cursor, name, readCount));
    } else {
    	search->second.addPosition(cursor, name, readCount);
    }
    ++cursor;
  }
}

VecStr kmerCalculator::getAllKmers(const std::string& seq, uint32_t kLength) {
  uint32_t cursor = 0;
  VecStr ans;
  while (cursor + kLength <= seq.size()) {
    ans.emplace_back(seq.substr(cursor, kLength));
    ++cursor;
  }
  return ans;
}
}  // namespace bib
