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
		switch (ioOptions_.outFormat_) {
			case SeqIOOptions::outFormats::FASTA:
				primaryOut_ = std::make_unique<std::ofstream>();
				openTextFile(*primaryOut_, ioOptions_.out_);
				break;
			case SeqIOOptions::outFormats::FASTQ:
				primaryOut_ = std::make_unique<std::ofstream>();
				if("" == ioOptions_.out_.outExtention_){
					bib::files::openTextFile(*primaryOut_, ioOptions_.out_.outFilename_, ".fastq",
							ioOptions_.out_.overWriteFile_, ioOptions_.out_.append_,
							ioOptions_.out_.exitOnFailureToWrite_);
				}else{
					openTextFile(*primaryOut_, ioOptions_.out_);
				}
				break;
			case SeqIOOptions::outFormats::FASTAQUAL:
				primaryOut_ = std::make_unique<std::ofstream>();
				bib::files::openTextFile(*primaryOut_, ioOptions_.out_.outFilename_, ".fasta",
						ioOptions_.out_.overWriteFile_, ioOptions_.out_.append_,
						ioOptions_.out_.exitOnFailureToWrite_);
				secondaryOut_ = std::make_unique<std::ofstream>();
				bib::files::openTextFile(*secondaryOut_, ioOptions_.out_.outFilename_,
						".fasta.qual", ioOptions_.out_.overWriteFile_, ioOptions_.out_.append_,
						ioOptions_.out_.exitOnFailureToWrite_);
				break;
			case SeqIOOptions::outFormats::FASTQPAIRED:
				primaryOut_ = std::make_unique<std::ofstream>();
				bib::files::openTextFile(*primaryOut_, ioOptions_.out_.outFilename_, "_R1.fastq",
						ioOptions_.out_.overWriteFile_, ioOptions_.out_.append_,
						ioOptions_.out_.exitOnFailureToWrite_);
				secondaryOut_ = std::make_unique<std::ofstream>();
				bib::files::openTextFile(*secondaryOut_, ioOptions_.out_.outFilename_,
						"_R2.fastq", ioOptions_.out_.overWriteFile_, ioOptions_.out_.append_,
						ioOptions_.out_.exitOnFailureToWrite_);
				break;
			case SeqIOOptions::outFormats::FLOW:
				primaryOut_ = std::make_unique<std::ofstream>();
				bib::files::openTextFile(*primaryOut_, ioOptions_.out_.outFilename_, ".dat",
						ioOptions_.out_.overWriteFile_, ioOptions_.out_.append_,
						ioOptions_.out_.exitOnFailureToWrite_);
				break;
			case SeqIOOptions::outFormats::FLOWMOTHUR:
				primaryOut_ = std::make_unique<std::ofstream>();
				bib::files::openTextFile(*primaryOut_, ioOptions_.out_.outFilename_, ".flow",
						ioOptions_.out_.overWriteFile_, ioOptions_.out_.append_,
						ioOptions_.out_.exitOnFailureToWrite_);
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
				"Error in readObjectIOOpt, attempted to write when out files aren't open" };
	}
	writeNoCheck(read);
}

void SeqOutput::writeNoCheck(const seqInfo & read) {
	switch (ioOptions_.outFormat_) {
	case SeqIOOptions::outFormats::FASTA:
		read.outPutSeq(*primaryOut_);
		break;
	case SeqIOOptions::outFormats::FASTQ:
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


	if (SeqIOOptions::outFormats::FASTQPAIRED == ioOptions_.outFormat_) {
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

void SeqOutput::seekpSec(size_t pos){
	secondaryOut_->seekp(pos);
}

}  // namespace bibseq
