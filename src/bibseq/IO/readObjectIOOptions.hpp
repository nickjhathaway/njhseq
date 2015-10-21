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
//  readObjectIOOptions.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 2/07/15.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//


#include <bibcpp/jsonUtils.h>

namespace bibseq {


struct readObjectIOOptions {
  readObjectIOOptions() {}
	readObjectIOOptions(const std::string & firstName,
			const std::string & inFormat, bool processed);

  std::string firstName_;
  std::string secondName_;
  std::string inFormat_;

  std::string outFormat_;
  std::string outFilename_;
  std::string outExtention_;
  bool processed_ = false;
  std::string lowerCaseBases_;
  bool removeGaps_ = false;
  bool overWriteFile_ = false;
  bool exitOnFailureToWrite_ = false;
  bool append_ = false;
  bool includeWhiteSpaceInName_ = true;
  int32_t extra_ = 0;

  /**@b Create from a json string
   *
   * @param jsonStr The value stored in json
   */
  readObjectIOOptions(const std::string & jsonStr);
  /**@b output options as json
   *
   * @return options represented in json
   */
  Json::Value toJson()const;
};



}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "readObjectIOOptions.cpp"
#endif
