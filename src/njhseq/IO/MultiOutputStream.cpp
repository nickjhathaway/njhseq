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

#include <unordered_map>

#include "MultiOutputStream.hpp"

namespace njhseq {

bool MultiOutputStream::containsStream(const std::string & uid) const {
	return outStreams_.find(uid) != outStreams_.end();
}

void MultiOutputStream::addStream(const std::string & uid,
		const OutOptions & opts) {
	if (containsStream(uid)) {
		std::stringstream ss;
		ss << njh::bashCT::bold << njh::bashCT::red << __PRETTY_FUNCTION__
				<< njh::bashCT::resetAdd(njh::bashCT::red)
				<< ": trying to add a reader that already exists: " << uid
				<< njh::bashCT::reset;
		throw std::runtime_error { ss.str() };
	} else {
		outStreams_.emplace(uid, std::make_unique<OutputStreamWrap>(opts));
	}
}



void MultiOutputStream::write(const std::string & uid, const std::string & line, bool flush) {
	containsStreamThrow(uid);
	std::lock_guard<std::mutex> lock(outStreams_.at(uid)->mut_);
	outStreams_.at(uid)->write(line, flush);
}



void MultiOutputStream::openWrite(const std::string & uid, const std::string & line, bool flush) {
	openOut(uid);
	std::lock_guard<std::mutex> lock(outStreams_.at(uid)->mut_);
	outStreams_.at(uid)->writeNoCheck(line, flush);
}

void MultiOutputStream::openWrite(const std::string & uid, const std::vector<std::string> & lines, bool flush) {
	openOut(uid);
	std::lock_guard<std::mutex> lock(outStreams_.at(uid)->mut_);
	for(const auto & line : lines){
		outStreams_.at(uid)->writeNoCheck(line, flush);
	}
}






void MultiOutputStream::closeOutAll() {
	for (auto & outStream : outStreams_) {
		std::lock_guard<std::mutex> lock(outStream.second->mut_);
		if (outStream.second->outOpen()) {
			outStream.second->closeOut();
		}
	}
}

void MultiOutputStream::closeOutForReopeningAll() {
	for (auto & outStream : outStreams_) {
		std::lock_guard<std::mutex> lock(outStream.second->mut_);
		if (outStream.second->outOpen()) {
			outStream.second->closeOutForReopening();
		}
	}
}

void MultiOutputStream::closeNext() {
	auto nextUp = outsOpen_.front();
	outsOpen_.pop_front();
	while (openPriorityCounts_[nextUp] > 1) {
		--openPriorityCounts_[nextUp];
		nextUp = outsOpen_.front();
		outsOpen_.pop_front();
	}
	auto outStream = outStreams_.find(nextUp);
	std::lock_guard<std::mutex> lock(outStream->second->mut_);
	outStream->second->closeOutForReopening();
	openPriorityCounts_[nextUp] = 0;
}

uint32_t MultiOutputStream::getOpenLimit() const {
	return outOpenLimit_;
}

void MultiOutputStream::setOpenLimit(uint32_t limit) {
	if(limit == 0){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", can't set limit to 0" << std::endl;
		throw std::runtime_error{ss.str()};
	}
	outOpenLimit_ = limit;
}

void MultiOutputStream::containsStreamThrow(const std::string & uid) const {
	if (!containsStream(uid)) {
		std::stringstream ss;
		ss << njh::bashCT::bold << njh::bashCT::red << __PRETTY_FUNCTION__
				<< njh::bashCT::resetAdd(njh::bashCT::red)
				<< ": trying to use reader that hasn't been added: " << uid
				<< njh::bashCT::reset;
		throw std::runtime_error { ss.str() };
	}
}


void MultiOutputStream::openOut(const std::string & uid){
	//std::cout << __PRETTY_FUNCTION__ << 1 << std::endl;
	//std::cout << uid << std::endl;
	containsStreamThrow(uid);
	//std::cout << __PRETTY_FUNCTION__ << 2 << std::endl;
	mut_.lock();
	std::lock_guard<std::mutex> lock(outStreams_.at(uid)->mut_);
	//std::cout << __PRETTY_FUNCTION__ << 3 << std::endl;
	auto outStream = outStreams_.find(uid);
	//std::cout << __PRETTY_FUNCTION__ << 4 << std::endl;
	if (!outStream->second->outOpen()) {
		//std::cout << __PRETTY_FUNCTION__ << 5 << std::endl;
		if (outCurrentlyOpen_ >= outOpenLimit_) {
			//std::cout << __PRETTY_FUNCTION__ << 6 << std::endl;
			closeNext();
			//std::cout << __PRETTY_FUNCTION__ << 7 << std::endl;
		} else {
			//std::cout << __PRETTY_FUNCTION__ << 8 << std::endl;
			++outCurrentlyOpen_;
			//std::cout << __PRETTY_FUNCTION__ << 9 << std::endl;
		}
		//std::cout << __PRETTY_FUNCTION__ << 10 << std::endl;
		outStream->second->openOut();
		//std::cout << __PRETTY_FUNCTION__ << 11 << std::endl;
	}
	//std::cout << __PRETTY_FUNCTION__ << 12 << std::endl;
	outsOpen_.push_back(uid);
	//std::cout << __PRETTY_FUNCTION__ << 13 << std::endl;
	++openPriorityCounts_[uid];
	//std::cout << __PRETTY_FUNCTION__ << 14 << std::endl;
	mut_.unlock();
	//std::cout << __PRETTY_FUNCTION__ << 15 << std::endl;
}


}  // namespace njhseq


