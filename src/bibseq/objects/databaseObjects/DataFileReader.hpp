#pragma once
/*
 *
 *  Created on: May 17, 2015
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"

namespace bibseq {




class DataFileReader {
public:

	DataFileReader(const std::string & filename);

	std::ifstream file_;
	std::string filename_;

	explicit operator bool() const ;

	template<typename DATA>
	bool readNextRecord(DATA & record){
		std::string line = "";
		if(std::getline(file_, line)){
			record = DATA(line);
			return true;
		}else{
			return false;
		}
	}

	virtual ~DataFileReader();


	// no-copy
	DataFileReader(const DataFileReader& other) = delete;
	DataFileReader() = delete;
	DataFileReader& operator=(const DataFileReader& other) = delete;
};

} /* namespace bibseq */



