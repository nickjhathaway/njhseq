#pragma once
//
//  IOUtils.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 9/14/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include <bibcpp/files.h>

class IoOptions {
public:

	IoOptions(const std::string & filename, bool out){
		if(out){
			setOutOptions(filename, "", "");
		}else{
			setInOptions(filename, "", "");
		}
	}

	IoOptions(const std::string & filename, const std::string & extention,
			const std::string & format, bool out) {
		if(out){
			setOutOptions(filename, extention, format);
		}else{
			setInOptions(filename, extention, format);
		}
	}

	IoOptions(const std::string & filename, const std::string & extention,
			const std::string & format, bool append, bool overWriteFile, bool exitOnFailureToWrite) {
		setOutOptions(filename, extention, format);
		setWritingOptions(append, overWriteFile, exitOnFailureToWrite);
	}

	std::string inFilename_;
	std::string inExtention_;
	std::string inFormat_;

  std::string outFilename_;
  std::string outExtention_;
  std::string outFormat_;

  bool append_ = false;
  bool overWriteFile_ = false;
  bool exitOnFailureToWrite_ = true;

  void setInOptions(const std::string & filename, const std::string & extention,
			const std::string & format) {
  	inFilename_ = filename;
  	inExtention_ = extention;
  	inFormat_ = format;
	}

  void setOutOptions(const std::string & filename, const std::string & extention,
			const std::string & format) {
  	outFilename_ = filename;
  	outExtention_ = extention;
  	outFormat_ = format;
	}

  void setWritingOptions(bool append, bool overWriteFile, bool exitOnFailureToWrite){
  	append_ = append;
  	overWriteFile_ = overWriteFile;
  	exitOnFailureToWrite_ = exitOnFailureToWrite;
  }

};


