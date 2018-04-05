#pragma once
/*
 * BioDataFileIO.hpp
 *
 *  Created on: Aug 1, 2016
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "bibseq/common.h"
#include "bibseq/IO.h"




namespace bibseq {


template<typename DATA>
class BioDataFileIO {
public:

	BioDataFileIO(const IoOptions & ioOpts):ioOpts_(ioOpts){}
	BioDataFileIO(const BioDataFileIO& other): ioOpts_(other.ioOpts_){

	}
	BioDataFileIO& operator=(const BioDataFileIO& other) = delete;

	const IoOptions ioOpts_;
	std::unique_ptr<InputStream> inFile_;
	std::unique_ptr<OutputStream> out_;
	std::string currentLine_{""};

	using value_type = DATA;

	void openIn(){
		if(isOutOpen() && ioOpts_.in_.inFilename_ == ioOpts_.out_.outName()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, can't open file " << ioOpts_.in_.inFilename_ << " for reading while it is open for writing" << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(bib::strToLowerRet(ioOpts_.in_.inFilename_.string()) != "stdin" &&  !ioOpts_.in_.inExists()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, file " << ioOpts_.in_.inFilename_ << " doesn't exist" << "\n";
			throw std::runtime_error{ss.str()};
		}
		inFile_ = std::make_unique<InputStream>(ioOpts_.in_);
		if(!inFile_->good()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error in opening " << ioOpts_.in_.inFilename_ << "\n";
			throw std::runtime_error{ss.str()};
		}
	}

	void openOut(){
		if(isInOpen() && ioOpts_.in_.inFilename_ == ioOpts_.out_.outName()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, can't open file " << ioOpts_.in_.inFilename_ << " for writing while it is open for reading" << "\n";
			throw std::runtime_error{ss.str()};
		}
		out_ = std::make_unique<OutputStream>(ioOpts_.out_);
	}

	void closeIn(){
		if(isInOpen()){
			inFile_ = nullptr;
		}
	}
	void closeOut(){
		if(isOutOpen()){
			out_ = nullptr;
		}
	}

	bool isInOpen() const{
		return nullptr != inFile_;
	}

	bool isOutOpen() const{
		return  nullptr != out_;
	}

	bool readNextRecord(DATA & record) {
		if(!isInOpen()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, infile_ not open " << std::endl;
			throw std::runtime_error{ss.str()};
		}
		std::string line = "";
		if (bib::files::crossPlatGetline(*inFile_, line)) {
			while ((line.front() == '#' || line.empty()) && !inFile_->eof()) {
				bib::files::crossPlatGetline(*inFile_, line);
			}
			if (line.front() == '#' || line.empty()) {
				return false;
			} else {
				record = DATA(line);
				currentLine_ = line;
				return true;
			}
		} else {
			return false;
		}
	}

	std::shared_ptr<DATA> readNextRecord() {
		if(!isInOpen()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, infile_ not open " << std::endl;
			throw std::runtime_error{ss.str()};
		}
		std::string line = "";
		if (bib::files::crossPlatGetline(*inFile_, line)) {
			while ((line.front() == '#' || line.empty()) && !inFile_->eof()) {
				bib::files::crossPlatGetline(*inFile_, line);
			}
			if (line.front() == '#' || line.empty()) {
				return nullptr;
			} else {
				currentLine_ = line;
				return std::make_shared<DATA>(line);
			}
		} else {
			return nullptr;
		}
	}

	void write(const DATA & record, std::function<void(const DATA & d,std::ostream & out)> func){
		if(!isOutOpen()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, outFile_ not open " << std::endl;
			throw std::runtime_error{ss.str()};
		}
		writeNoCheck(record, func);
	}

	void openWrite(const DATA & record, std::function<void(const DATA & d,std::ostream & out)> func){
		if(!isOutOpen()){
			openOut();
		}
		writeNoCheck(record, func);
	}

	void writeNoCheck(const DATA & record, std::function<void(const DATA & d,std::ostream & out)> func){
		func(record, *out_);
	}

	void write(const std::vector<DATA> & records, std::function<void(const DATA & d,std::ostream & out)> func){
		if(!isOutOpen()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, outFile_ not open " << std::endl;
			throw std::runtime_error{ss.str()};
		}
		for(const auto & record : records){
			writeNoCheck(record, func);
		}
	}

	void openWrite(const std::vector<DATA> & records, std::function<void(const DATA & d,std::ostream & out)> func){
		if(!isOutOpen()){
			openOut();
		}
		for(const auto & record : records){
			writeNoCheck(record, func);
		}
	}


	virtual ~BioDataFileIO(){
		inFile_ = nullptr;
		out_ = nullptr;
	}




};





}  // namespace bibseq
