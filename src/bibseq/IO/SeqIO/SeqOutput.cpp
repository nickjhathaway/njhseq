//
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
/*
 * SeqOutput.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */


#include "SeqOutput.hpp"
#include <bibcpp/bashUtils.h>
#include "bibseq/IO/fileUtils.hpp"
namespace bibseq {
SeqOutput::SeqOutput(const SeqIOOptions & options) :
				ioOptions_(options) {}


bool SeqOutput::outOpen() const {
	return outOpen_;
}

void SeqOutput::openOut() {
	if(!outOpen_){
		std::stringstream ss;
		auto prim_outOpts = ioOptions_.out_;
		auto sec_outOpts = ioOptions_.out_;
		switch (ioOptions_.outFormat_) {
			case SeqIOOptions::outFormats::FASTA:
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				break;
			case SeqIOOptions::outFormats::FASTAGZ:
				if(!bib::endsWith(prim_outOpts.outExtention_, ".gz")){
					prim_outOpts.outExtention_.append(".gz");
				}
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				break;
			case SeqIOOptions::outFormats::FASTQ:
				if("" == prim_outOpts.outExtention_){
					prim_outOpts.outExtention_ = ".fastq";
				}
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				break;
			case SeqIOOptions::outFormats::FASTQGZ:
				if("" == prim_outOpts.outExtention_){
					prim_outOpts.outExtention_ = ".fastq.gz";
				}
				if(!bib::endsWith(prim_outOpts.outExtention_, ".gz")){
					prim_outOpts.outExtention_.append(".gz");
				}
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				break;
			case SeqIOOptions::outFormats::FASTAQUAL:
				prim_outOpts.outExtention_ = ".fasta";
				sec_outOpts.outExtention_ = ".fasta.qual";
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				secondaryOut_ = std::make_unique<OutputStream>(sec_outOpts);
				break;
			case SeqIOOptions::outFormats::FASTQPAIRED:
				prim_outOpts.outExtention_ = "_R1.fastq";
				sec_outOpts.outExtention_ = "_R2.fastq";
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				secondaryOut_ = std::make_unique<OutputStream>(sec_outOpts);
				break;
			case SeqIOOptions::outFormats::FASTQPAIREDGZ:
				prim_outOpts.outExtention_ = "_R1.fastq.gz";
				sec_outOpts.outExtention_ = "_R2.fastq.gz";
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				secondaryOut_ = std::make_unique<OutputStream>(sec_outOpts);
				break;
			case SeqIOOptions::outFormats::FLOW:
				prim_outOpts.outExtention_ = ".dat";
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				break;
			case SeqIOOptions::outFormats::FLOWMOTHUR:
				prim_outOpts.outExtention_ = ".flow";
				primaryOut_ = std::make_unique<OutputStream>(prim_outOpts);
				break;
			default:
				ss << "Unrecognized out file type : in" << __PRETTY_FUNCTION__
						<< std::endl;
				ss << "Acceptable types are fasta, fastq, and fastaQual" << std::endl;
				throw std::runtime_error { ss.str() };
				break;
		}
		outOpen_ = true;
	}
}

void SeqOutput::closeOut() {
	primaryOut_ = nullptr;
	secondaryOut_ = nullptr;
	outOpen_ = false;
}

void SeqOutput::closeOutForReopening() {
	closeOut();
	ioOptions_.out_.append_ = true;
}



void SeqOutput::write(const seqInfo & read) {
	if (!outOpen_) {
		throw std::runtime_error {
				"Error in readObjectIOOpt, attempted to write when out files aren't open, out file: " + ioOptions_.out_.outName().string() };
	}
	writeNoCheck(read);
}

void SeqOutput::writeNoCheck(const seqInfo & read) {
	switch (ioOptions_.outFormat_) {
	case SeqIOOptions::outFormats::FASTA:
	case SeqIOOptions::outFormats::FASTAGZ:
		read.outPutSeq(*primaryOut_);
		break;
	case SeqIOOptions::outFormats::FASTQ:
	case SeqIOOptions::outFormats::FASTQGZ:
		read.outPutFastq(*primaryOut_);
		break;
	case SeqIOOptions::outFormats::FASTAQUAL:
		read.outPutSeq(*primaryOut_);
		read.outPutQual(*secondaryOut_);
		break;
	default:
		throw std::runtime_error { bib::bashCT::boldRed(
				"in " + std::string(__PRETTY_FUNCTION__)
						+ " : unrecognized out format case") };
		break;
	}
}

void SeqOutput::openWrite(const seqInfo & read) {
	if (!outOpen_) {
		openOut();
	}
	writeNoCheck(read);
}

void SeqOutput::writeNoCheck(const PairedRead & seq) {
	if (SeqIOOptions::outFormats::FASTQPAIRED == ioOptions_.outFormat_ || SeqIOOptions::outFormats::FASTQPAIREDGZ == ioOptions_.outFormat_) {
		seq.seqBase_.outPutFastq(*primaryOut_);
		seqInfo mateInfo = seq.mateSeqBase_;
		if(seq.mateRComplemented_){
			mateInfo.reverseComplementRead(false, true);
		}
		mateInfo.outPutFastq(*secondaryOut_);
	} else {
		throw std::runtime_error {
				"Error in " + std::string(__PRETTY_FUNCTION__)
						+ ", attempted to write paired read when output option is not fastqPaired"};
	}
}

void SeqOutput::write(const PairedRead & read) {
	if (!outOpen_) {
		throw std::runtime_error { "Error in " + std::string(__PRETTY_FUNCTION__)
				+ ", attempted to write when out files aren't open" };
	}
	writeNoCheck(read);
}

void SeqOutput::openWrite(const PairedRead & read) {
	if (!outOpen_) {
		openOut();
	}
	writeNoCheck(read);
}

void SeqOutput::openWriteFlow(const sffObject & read) {
	if (!outOpen_) {
		openOut();
	}
	if (SeqIOOptions::outFormats::FLOW == ioOptions_.outFormat_
			|| SeqIOOptions::outFormats::FLOWMOTHUR == ioOptions_.outFormat_) {
		read.outPutFlowValuesPryo(*primaryOut_);
	} else {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< ", attempting to write flow data when the output format is not FLOW or FLOWMOTHUR"<< std::endl;
	}
}

void SeqOutput::writeNoCheckFlow(const sffObject & read) {
	read.outPutFlowValuesPryo(*primaryOut_);
}

void SeqOutput::openWrite(const std::string & line){
	if (!outOpen_) {
		openOut();
	}
	(*primaryOut_) << line << std::endl;
}

void SeqOutput::writeNoCheck(const std::string & line){
	(*primaryOut_) << line << std::endl;
}

size_t SeqOutput::tellpPri(){
	return primaryOut_->tellp();
}

void SeqOutput::seekpPri(size_t pos){
	primaryOut_->seekp(pos);
}

size_t SeqOutput::tellpSec(){
	return secondaryOut_->tellp();
}
//Error in readObjectIOOpt, attempted to write when out files aren't open
void SeqOutput::seekpSec(size_t pos){
	secondaryOut_->seekp(pos);
}

bfs::path SeqOutput::getPrimaryOutFnp() const{
	if(nullptr == primaryOut_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error primaryOut_ is not currently set" << "\n";
		throw std::runtime_error{ss.str()};
	}
	return primaryOut_->outOpts_.outName();
}

bfs::path SeqOutput::getSecondaryOutFnp() const{
	if(nullptr == secondaryOut_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error primaryOut_ is not currently set" << "\n";
		throw std::runtime_error{ss.str()};
	}
	return secondaryOut_->outOpts_.outName();

}

}  // namespace bibseq
