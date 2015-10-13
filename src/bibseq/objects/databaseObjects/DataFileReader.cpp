/*
 * DataFileReader.cpp
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */

#include "DataFileReader.hpp"

namespace bibseq {

DataFileReader::DataFileReader(const std::string & filename):filename_(filename){
	file_.open(filename, std::ios::in);
}

DataFileReader::operator bool() const {
	return file_.is_open();
}
DataFileReader::~DataFileReader(){
	file_.close();
}


} /* namespace bibseq */
