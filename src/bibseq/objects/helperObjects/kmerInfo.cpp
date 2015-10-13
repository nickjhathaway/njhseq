#include "kmerInfo.hpp"
#include "bibseq/helpers/seqUtil.hpp"

namespace bibseq {
kmerInfo::kmerInfo(const std::string & seq, uint32_t kLength, bool setReverse): kLen_(kLength), seqLen_(seq.size()), infoSet_(true){
	setKmers(seq,kLength, setReverse);
}


void kmerInfo::setKmers(const std::string & seq, uint32_t kLength, bool setReverse){
	infoSet_ = true;
	kLen_ = kLength;
	seqLen_ = seq.size();
	kmers_.clear();
	for(const auto & pos : iter::range(seq.size() + 1 - kLen_)){
		auto currentK = seq.substr(pos, kLen_);
		auto k = kmers_.find(currentK);
		if(k!= kmers_.end()){
			k->second.addPosition(pos);
		}else{
			kmers_[currentK] = kmer(currentK, pos);
		}
	}
	if(setReverse){
		kmersRevComp_.clear();
		std::string reverseComplement = seqUtil::reverseComplement(seq, "DNA");
		for(const auto & pos : iter::range(reverseComplement.size() + 1 - kLen_)){
			auto currentK = reverseComplement.substr(pos, kLen_);
			auto k = kmersRevComp_.find(currentK);
			if(k!= kmersRevComp_.end()){
				k->second.addPosition(pos);
			}else{
				kmersRevComp_[currentK] = kmer(currentK, pos);
			}
		}
	}
}

std::pair<uint32_t, double> kmerInfo::compareKmers(const kmerInfo & info) const{
	uint32_t kShared = 0;
	//uint32_t kLen = kmers_.begin()->first.size();
	for(const auto & k : kmers_){
		auto otherK = info.kmers_.find(k.first);
		if(otherK!= info.kmers_.end()){
			kShared += std::min(otherK->second.count_, k.second.count_);
		}
	}
	//std::cout << info.seqLen_ << std::endl;
	//std::cout << seqLen_ << std::endl;
	return {kShared,
		kShared/static_cast<double>(std::min(info.seqLen_, seqLen_) + 1 - kLen_)};
}

std::pair<uint32_t, double> kmerInfo::compareKmersRevComp(const kmerInfo & info) const{
	uint32_t kShared = 0;
	//uint32_t kLen = kmers_.begin()->first.size();
	for(const auto & k : kmers_){
		auto otherK = info.kmersRevComp_.find(k.first);
		if(otherK!= info.kmersRevComp_.end()){
			kShared += std::min(otherK->second.count_, k.second.count_);
		}
	}
	return {kShared,
		kShared/static_cast<double>(std::min(info.seqLen_,
				seqLen_) + 1 - kLen_)};
}

std::pair<uint32_t, double> kmerInfo::compareKmers(const kmerInfo & info,
		uint32_t startPos, uint32_t windowSize) const{
	uint32_t kmersShared = 0;
	double maxKMers = windowSize - kLen_ + 1;
	for(const auto & k : kmers_){
		uint32_t kmersInWindow = 0;
		for(const auto & pos : k.second.positions_){
			if(pos >= startPos && pos <= (startPos + windowSize - kLen_)){
				++kmersInWindow;
			}
			if (pos > (startPos + windowSize - kLen_)){
				break;
			}
		}
		uint32_t kmersInWindowOther = 0;
		auto otherK = info.kmers_.find(k.first);
		if(otherK != info.kmers_.end()){
			for(const auto & pos : otherK->second.positions_){
				if(pos >= startPos && pos <= (startPos + windowSize - kLen_)){
					++kmersInWindowOther;
				}
				if (pos > (startPos + windowSize - kLen_)){
					break;
				}
			}
		}
		kmersShared += std::min(kmersInWindow, kmersInWindowOther);
	}
	return {kmersShared, kmersShared/maxKMers};
}

std::pair<uint32_t, double> kmerInfo::compareSubKmersToFull(const kmerInfo & info,
		uint32_t startPos, uint32_t windowSize) const{
	uint32_t kmersShared = 0;
	double maxKMers = windowSize - kLen_ + 1;
	for(const auto & k : kmers_){
		uint32_t kmersInWindow = 0;
		for(const auto & pos : k.second.positions_){
			if(pos >= startPos && pos <= (startPos + windowSize - kLen_)){
				++kmersInWindow;
			}
			if (pos > (startPos + windowSize - kLen_)){
				break;
			}
		}
		uint32_t kmersInWindowOther = 0;
		auto otherK = info.kmers_.find(k.first);
		if(otherK != info.kmers_.end()){
			kmersInWindowOther = otherK->second.count_;
		}
		kmersShared += std::min(kmersInWindow, kmersInWindowOther);
	}
	return {kmersShared, kmersShared/maxKMers};
}

std::vector<std::pair<uint32_t, double>> kmerInfo::slideCompareKmers(const kmerInfo & info,
		uint32_t windowSize, uint32_t windowStepSize
		) const {
	std::vector<std::pair<uint32_t, double>> ret;
	uint64_t minLen = std::min(seqLen_, info.seqLen_);
	for(const auto & pos : iter::range<uint32_t>(0, minLen - windowSize + 1, windowStepSize)){
		ret.emplace_back(compareKmers(info, pos, windowSize));
	}
	return ret;
}

std::vector<std::pair<uint32_t, double>> kmerInfo::slideCompareSubKmersToFull(const kmerInfo & info,
		uint32_t windowSize, uint32_t windowStepSize
		) const {
	std::vector<std::pair<uint32_t, double>> ret;
	for(const auto & pos : iter::range<uint32_t>(0,seqLen_ - windowSize + 1, windowStepSize)){
		ret.emplace_back(compareSubKmersToFull(info, pos, windowSize));
	}
	return ret;
}



}  // namespace bib
