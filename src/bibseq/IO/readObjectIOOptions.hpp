#pragma once
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

  std::string firstName_;
  std::string secondName_;
  std::string thirdName_;
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
