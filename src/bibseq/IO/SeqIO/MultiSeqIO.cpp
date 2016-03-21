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
/*
 * MultiSeqIO.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */


#include "MultiSeqIO.hpp"

namespace bibseq {

bool MultiSeqIO::containsReader(const std::string & uid) const {
	return readIos_.find(uid) != readIos_.end();
}

void MultiSeqIO::addReader(const std::string & uid,
		const SeqIOOptions & opts) {
	if (containsReader(uid)) {
		std::stringstream ss;
		ss << bib::bashCT::bold << bib::bashCT::red << __PRETTY_FUNCTION__
				<< bib::bashCT::resetAdd(bib::bashCT::red)
				<< ": trying to add a reader that already exists: " << uid
				<< bib::bashCT::reset;
		throw std::runtime_error { ss.str() };
	} else {
		readIos_.emplace(uid, std::make_unique<SeqIO>(opts));
	}
}


void MultiSeqIO::openWrite(const std::string & uid,
		const std::string & line) {
	openOut(uid);
	std::lock_guard<std::mutex> lock(readIos_.at(uid)->mut_);
		auto readIO = readIos_.find(uid);
	readIO->second->out_.writeNoCheck(line);
}

void MultiSeqIO::openWriteFlow(const std::string & uid,
		const sffObject & read) {
	openOut(uid);
	std::lock_guard<std::mutex> lock(readIos_.at(uid)->mut_);
	auto readIO = readIos_.find(uid);
	readIO->second->out_.writeNoCheckFlow(read);
}

void MultiSeqIO::closeInAll() {
	for (auto & readIo : readIos_) {
		std::lock_guard<std::mutex> lock(readIo.second->mut_);
		if (readIo.second->in_.inOpen()) {
			readIo.second->in_.closeIn();
		}
	}
}

void MultiSeqIO::closeOutAll() {
	for (auto & readIo : readIos_) {
		std::lock_guard<std::mutex> lock(readIo.second->mut_);
		if (readIo.second->out_.outOpen()) {
			readIo.second->out_.closeOut();
		}
	}
}

void MultiSeqIO::closeOutForReopeningAll() {
	for (auto & readIo : readIos_) {
		std::lock_guard<std::mutex> lock(readIo.second->mut_);
		if (readIo.second->out_.outOpen()) {
			readIo.second->out_.closeOutForReopening();
		}
	}
}

void MultiSeqIO::closeNext() {
	auto nextUp = outsOpen_.front();
	outsOpen_.pop_front();
	while (openPriorityCounts_[nextUp] > 1) {
		--openPriorityCounts_[nextUp];
		nextUp = outsOpen_.front();
		outsOpen_.pop_front();
	}
	auto readIO = readIos_.find(nextUp);
	std::lock_guard<std::mutex> lock(readIO->second->mut_);
	readIO->second->out_.closeOutForReopening();
	openPriorityCounts_[nextUp] = 0;
}

uint32_t MultiSeqIO::getOpenLimit() const {
	return outOpenLimit_;
}

void MultiSeqIO::setOpenLimit(uint32_t limit) {
	if(limit == 0){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", can't set limit to 0" << std::endl;
		throw std::runtime_error{ss.str()};
	}
	outOpenLimit_ = limit;
}

void MultiSeqIO::containsReaderThrow(const std::string & uid) const {
	if (!containsReader(uid)) {
		std::stringstream ss;
		ss << bib::bashCT::bold << bib::bashCT::red << __PRETTY_FUNCTION__
				<< bib::bashCT::resetAdd(bib::bashCT::red)
				<< ": trying to use reader that hasn't been added: " << uid
				<< bib::bashCT::reset;
		throw std::runtime_error { ss.str() };
	}
}


void MultiSeqIO::openOut(const std::string & uid){
	containsReaderThrow(uid);
	mut_.lock();
	std::lock_guard<std::mutex> lock(readIos_.at(uid)->mut_);
	auto readIO = readIos_.find(uid);
	if (!readIO->second->out_.outOpen()) {
		if (outCurrentlyOpen_ >= outOpenLimit_) {
			closeNext();
		} else {
			++outCurrentlyOpen_;
		}
		readIO->second->out_.openOut();
	}
	outsOpen_.push_back(uid);
	++openPriorityCounts_[uid];
	mut_.unlock();
}


}  // namespace bibseq


