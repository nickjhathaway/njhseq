/*
 * BamAlnsCache.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: nick
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "BamAlnsCache.hpp"


namespace bibseq {

void BamAlnsCache::add(const BamTools::BamAlignment & aln){
	cache_[aln.Name] = std::make_shared<BamTools::BamAlignment>(aln);
}

void BamAlnsCache::remove(const std::string & name){
	cache_.erase(name);
}


VecStr BamAlnsCache::getNames() const{
	return bib::getVecOfMapKeys(cache_);
}


std::shared_ptr<BamTools::BamAlignment> BamAlnsCache::get(
		const std::string & name) {
	if (has(name)) {
		return cache_[name];
	} else {
		return nullptr;
	}
}

bool BamAlnsCache::has(const std::string & name) const {
	return bib::has(cache_, name);
}


std::unordered_map<int32_t, int32_t> BamAlnsCache::getMinPos() const {
	std::unordered_map<int32_t, int32_t> ret;
	for (const auto & aln : cache_) {
		if (bib::has(ret, aln.second->RefID)) {
			if (aln.second->Position < ret[aln.second->RefID]) {
				ret[aln.second->RefID] = aln.second->Position;
			}
		} else {
			ret[aln.second->RefID] = aln.second->Position;
		}
	}
	return ret;
}

std::unordered_map<int32_t, int32_t> BamAlnsCache::getMaxEndPos() const {
	std::unordered_map<int32_t, int32_t> ret;
	for (const auto & aln : cache_) {
		if (bib::has(ret, aln.second->RefID)) {
			if (aln.second->GetEndPosition() > ret[aln.second->RefID]) {
				ret[aln.second->RefID] = aln.second->GetEndPosition();
			}
		} else {
			ret[aln.second->RefID] = aln.second->GetEndPosition();
		}
	}
	return ret;
}

BamAlnsCache::~BamAlnsCache(){

}


}  // namespace bibseq
