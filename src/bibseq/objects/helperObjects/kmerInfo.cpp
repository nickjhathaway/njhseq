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
