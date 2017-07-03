#pragma once
/*
 * OutputFile.hpp
 *
 *  Created on: Jul 3, 2017
 *      Author: nick
 */

#include "bibseq/IO/IOOptions/OutOptions.hpp"
#include <bibcpp/files.h>

namespace bibseq {

class OutputStream {
public:
	OutputStream(const OutOptions & outOpts);

	const OutOptions outOpts_;


	std::unique_ptr<bib::GZSTREAM::ogzstream> outFileGz_;
	std::unique_ptr<std::ofstream> outFile_;
	std::unique_ptr<std::ostream> out_;


	void write(const char * data, size_t numBytes);

};

template<typename T>
OutputStream & operator <<(OutputStream & out, const T & val){
	(*out.out_) << val;
}

inline OutputStream & flush(OutputStream & out){
	out.out_->flush();
	return out;
}

inline OutputStream & endl(OutputStream & out){
	(*out.out_) << std::endl;
	return out;
}






} /* namespace bibseq */


