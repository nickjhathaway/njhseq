#pragma once
/*
 * BioDataFileIO.hpp
 *
 *  Created on: Aug 1, 2016
 *      Author: nick
 */


#include "bibseq/IO/fileUtils.hpp"




namespace bibseq {


template<typename DATA>
class BioDataFileIO {
public:

	BioDataFileIO(const IoOptions & ioOpts):ioOpts_(ioOpts){}
	BioDataFileIO(const BioDataFileIO& other): ioOpts_(other.ioOpts_){

	}
	BioDataFileIO& operator=(const BioDataFileIO& other) = delete;

	const IoOptions ioOpts_;
	std::unique_ptr<std::ifstream> inFile_;
	std::unique_ptr<std::ofstream> outFile_;

	using value_type = DATA;

	void openIn(){
		if(!ioOpts_.in_.inExists()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, file " << ioOpts_.in_.inFilename_ << " doesn't exist" << "\n";
			throw std::runtime_error{ss.str()};
		}
		inFile_ = std::make_unique<std::ifstream>(ioOpts_.in_.inFilename_);
		if(!inFile_->is_open()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error in opening " << ioOpts_.in_.inFilename_ << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	void openOut(){
		outFile_ = std::make_unique<std::ofstream>();
		openTextFile(*outFile_, ioOpts_.out_);
	}

	bool isInOpen() const{
		return nullptr != inFile_ && inFile_->is_open();
	}
	bool isOutOpen() const{
		return nullptr != outFile_ && outFile_->is_open();
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
		func(record, *outFile_);
	}


	virtual ~BioDataFileIO(){
		inFile_ = nullptr;
		outFile_ = nullptr;
	}



};


}  // namespace bibseq
