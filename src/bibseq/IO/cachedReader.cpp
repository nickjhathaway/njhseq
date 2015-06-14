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
#include "cachedReader.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {

void cachedReader::seek(uint64_t filePosition){
	is_.seekg(filePosition);
	lineBuffer_ = VecStr(bufferMax_);
	refillBuffer();
}
cachedReader::cachedReader(std::istream & is): is_(is){
	lineBuffer_ = VecStr(bufferMax_);
	refillBuffer();
}

bool cachedReader::refillBuffer(){
	//std::cout << bib::bashCT::cyan << "refillBuffer" << bib::bashCT::reset  << std::endl;
	lineNum_ = 0;
	std::string tempStr;
	while(lineNum_ < bufferMax_ && std::getline(is_, tempStr)){
		lineBuffer_[lineNum_] = tempStr;
		//std::cout << "\tlineNum_: " << lineNum_ << std::endl;
		//std::cout << "\t\t" << lineBuffer_[lineNum_] << std::endl;
		//std::cout << bib::colorBool(is_.eof()) <<std::endl;
		++lineNum_;
	}

	bufferPos_ = 0;
	if(lineNum_ == 0){
		//didn't read in anything, probably at end of the file
		if(!is_.eof()){
			throw std::runtime_error{"error in reading during cachedReader::refillBuffer()"};
		}
		doneReadering_ = true;
		//std::cout << "refillBuffer success : " << bib::colorBool(false) << std::endl;
		return false;
	}else{
		//std::cout << "refillBuffer success : " << bib::colorBool(true) << std::endl;
		// did read in some in the buffer
		return true;
	}
}

const std::string & cachedReader::currentLine(){
	//std::cout << bib::bashCT::bold << "Reading currentLine: " << bib::bashCT::reset  << std::endl;
	if(bufferPos_ >= bufferMax_){
		refillBuffer();
	}
	//std::cout << bib::bashCT::blue << "currentLine: " << bufferPos_ << bib::bashCT::reset << std::endl;
	return lineBuffer_[bufferPos_];
}

bool cachedReader::setNextLine(){
	//std::cout << bib::bashCT::bold << "setting next line: " << bib::bashCT::reset  << std::endl;
	++bufferPos_;
	//std::cout << bufferPos_ << ":" << lineNum_ << " of " << lineBuffer_.size()<< std::endl;
	if(bufferPos_ < lineNum_){
		if(bufferPos_ == lineNum_ && is_.eof()){
			doneReadering_ = true;
		}
		return true;
	}else{
		return refillBuffer();
	}
}

bool cachedReader::done(){
	return doneReadering_;
}



}  // namespace bibseq
