#include "kmerInfo.hpp"
#include "bibseq/helpers/seqUtil.hpp"

namespace bibseq {
kmerInfo::kmerInfo(const std::string & seq, uint32_t kLength): kLen_(kLength), seqLen_(seq.size()), infoSet_(true){
	setKmers(seq,kLength);
}


void kmerInfo::setKmers(const std::string & seq,uint32_t kLength){
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
}

std::pair<uint32_t, double> kmerInfo::compareKmers(const kmerInfo & info) const{
	uint32_t kShared = 0;
	for(const auto & k : kmers_){
		auto otherK = info.kmers_.find(k.first);
		if(otherK!= info.kmers_.end()){
			kShared += std::min(otherK->second.count_, k.second.count_);
		}
	}
	return {kShared,
		kShared/static_cast<double>(std::min(info.seqLen_, seqLen_) + 1 - kLen_)};
}



}  // namespace bib
