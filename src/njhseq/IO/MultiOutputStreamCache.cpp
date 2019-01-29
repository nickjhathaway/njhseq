//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
/*
 * MultiOutputStream.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: nick
 */

#include "MultiOutputStreamCache.hpp"

namespace njhseq {


void MultiOutputStreamCache::addOutputStream(const std::string & uid, const OutOptions & opts) {
	writers_.addStream(uid, opts);
}


void MultiOutputStreamCache::add(const std::string & uid, const std::string & line) {
	writers_.containsStreamThrow(uid);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	cache_[uid].push_back(line);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	++cacheSize_;
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if (cacheSize_ == cacheLimit_) {
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		writeCache();
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
}

void MultiOutputStreamCache::add(const std::string & uid, const std::vector<std::string> & lines) {
	writers_.containsStreamThrow(uid);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	for (const auto & line : lines) {
		cache_[uid].push_back(line);
		++cacheSize_;
		if (cacheSize_ == cacheLimit_) {
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			writeCache();
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
}

void MultiOutputStreamCache::writeCache() {
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	for (const auto & lines : cache_) {
		if(lines.second.empty()){
			continue;
		}

		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		writers_.openWrite(lines.first, lines.second);
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	for (auto & lines : cache_) {
		lines.second.clear();
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	cacheSize_ = 0;
}


void MultiOutputStreamCache::closeOutAll() {
	writeCache();
	writers_.closeOutAll();
}


void MultiOutputStreamCache::closeOutForReopeningAll() {
	writeCache();
	writers_.closeOutForReopeningAll();
}


uint32_t MultiOutputStreamCache::getCacheLimit() const {
	return cacheLimit_;
}


void MultiOutputStreamCache::setCacheLimit(uint32_t cacheLimit) {
	cacheLimit_ = cacheLimit;
}


void MultiOutputStreamCache::setOpenLimit(uint32_t fileOpenLimit){
	writers_.setOpenLimit(fileOpenLimit);
}


MultiOutputStreamCache::~MultiOutputStreamCache(){
	closeOutForReopeningAll();
}

void MultiOutputStreamCache::containsReaderThrow(const std::string & uid) const {
	writers_.containsStreamThrow(uid);
}

bool MultiOutputStreamCache::containsReader(const std::string & uid) const {
	return writers_.containsStream(uid);
}


}  // namespace njhseq
