#pragma once
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
//
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


