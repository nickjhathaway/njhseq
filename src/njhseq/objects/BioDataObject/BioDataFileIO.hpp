#pragma once
/*
 * BioDataFileIO.hpp
 *
 *  Created on: Aug 1, 2016
 *      Author: nick
 */
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
#include "njhseq/common.h"
#include "njhseq/IO.h"




namespace njhseq {


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

	std::string possibleHeader_{""};

	using value_type = DATA;

	void openIn(){
		if(isOutOpen() && ioOpts_.in_.inFilename_ == ioOpts_.out_.outName()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, can't open file " << ioOpts_.in_.inFilename_ << " for reading while it is open for writing" << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(njh::strToLowerRet(ioOpts_.in_.inFilename_.string()) != "stdin" &&  !ioOpts_.in_.inExists()){
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
		{
			//add input file header if it begins with a #
			if(isInOpen() && njh::strToLowerRet(ioOpts_.in_.inFilename_.string()) != "stdin"){
				InputStream tempIn(ioOpts_.in_);
				std::string line;
				njh::files::crossPlatGetline(tempIn, line);
				while(njh::beginsWith(line, "#")){
					(*out_) << line << std::endl;
					njh::files::crossPlatGetline(tempIn, line);
				}
			}
		}
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

	bool isOutOpen() const {
		return nullptr != out_;
	}

	std::vector<DATA> readAll() {
		std::vector<DATA> ret;
		DATA record;
		openIn();
		while (readNextRecord(record)) {
			ret.emplace_back(record);
		}
		closeIn();
		return ret;
	}

	std::vector<std::shared_ptr<DATA>> readAllPtrs() {
		std::vector<std::shared_ptr<DATA>> ret;
		std::shared_ptr<DATA> record = readNextRecord();
		openIn();
		while (nullptr != record) {
			ret.emplace_back(std::make_shared<DATA>(*record));
			record = readNextRecord();
		}
		closeIn();
		return ret;
	}


	bool readNextRecord(DATA & record) {
		if(!isInOpen()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, infile_ not open " << std::endl;
			throw std::runtime_error{ss.str()};
		}
		std::string line = "";
		if (njh::files::crossPlatGetline(*inFile_, line)) {
			if(line.front() == '#' && "" == possibleHeader_){
				possibleHeader_ = line;
			}
			while ((line.front() == '#' || line.empty()) && !inFile_->eof()) {
				njh::files::crossPlatGetline(*inFile_, line);
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
		if (njh::files::crossPlatGetline(*inFile_, line)) {
			if(line.front() == '#' && "" == possibleHeader_){
				possibleHeader_ = line;
			}
			while ((line.front() == '#' || line.empty()) && !inFile_->eof()) {
				njh::files::crossPlatGetline(*inFile_, line);
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





}  // namespace njhseq
