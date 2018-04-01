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
#include "BamAlnsCacheWithRegion.hpp"


namespace bibseq {

void BamAlnsCacheWithRegion::addWithRegion(const BamTools::BamAlignment & aln, const GenomicRegion & region){
	cache_[aln.Name] = std::make_shared<BamTools::BamAlignment>(aln);
	regionCache_[aln.Name] = std::make_shared<GenomicRegion>(region);
}

void BamAlnsCacheWithRegion::addWithRegion(const BamTools::BamAlignment & aln, const std::shared_ptr<GenomicRegion> & region){
	cache_[aln.Name] = std::make_shared<BamTools::BamAlignment>(aln);
	regionCache_[aln.Name] = region;
}

void BamAlnsCacheWithRegion::remove(const std::string & name){
	cache_.erase(name);
	regionCache_.erase(name);
}

std::shared_ptr<GenomicRegion> BamAlnsCacheWithRegion::getRegion(
		const std::string & name) {
	if (has(name)) {
		return regionCache_[name];
	} else {
		return nullptr;
	}
}

BamAlnsCacheWithRegion::~BamAlnsCacheWithRegion(){

}



}  // namespace bibseq
