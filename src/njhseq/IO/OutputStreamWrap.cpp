/*
 * OutputStreamWrap.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: nicholashathaway
 */


#include "OutputStreamWrap.hpp"

namespace njhseq {

OutputStreamWrap::OutputStreamWrap(const OutOptions & opts) :
		opts_(opts) {
}



void OutputStreamWrap::writeNoCheck(const std::string & line, bool flush){
	(*out_) << line << '\n';
	if(flush){
		out_->flush();
	}
}

void OutputStreamWrap::write(const std::string & line, bool flush){
	if(!outOpen_){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " trying to write when file not open"<< "\n";
		throw std::runtime_error{ss.str()};
	}
	writeNoCheck(line, flush);
}


bool OutputStreamWrap::outOpen() const{
	return outOpen_;
}

void OutputStreamWrap::openOut(){
	out_ = std::make_unique<OutputStream>(opts_);
	outOpen_ = true;
}

void OutputStreamWrap::closeOut(){
	outOpen_ = false;
	out_ = nullptr;
}
void OutputStreamWrap::closeOutForReopening(){
	opts_.append_ = true;
	closeOut();
}

}  // namespace njhseq

