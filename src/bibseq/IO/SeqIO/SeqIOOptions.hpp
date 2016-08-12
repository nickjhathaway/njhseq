#pragma once
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
//
//  readObjectIOOptions.hpp
//  sequenceTools
//
//  Created by Nick Hathaway on 2/07/15.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//


#include <bibcpp/jsonUtils.h>
#include "bibseq/IO/IOUtils.hpp"

namespace bibseq {


struct SeqIOOptions {
  enum class inFormats {
  	FASTQ,
		FASTQPAIRED,
		FASTQGZ,
		FASTA,
		FASTAGZ,
		FASTAQUAL,
		BAM,
		SFFTXT,
		SFFBIN,
		NOFORMAT
  };

  enum class outFormats {
  	FASTQ,
		FASTQPAIRED,
		FASTA,
		FASTAQUAL,
		FLOW,
		FLOWMOTHUR,
		NOFORMAT
  };

  SeqIOOptions();

  SeqIOOptions(const std::string & firstName,
  		inFormats inFormat, bool processed);

	SeqIOOptions(const OutOptions & out, outFormats outFormat);

	SeqIOOptions(const std::string & outFilename, outFormats outFormat, const OutOptions & out);

	/**@b Create from a json string
   *
   * @param jsonStr The value stored in json
   */
  explicit SeqIOOptions(const std::string & jsonStr);


  bool inExists() const;
	bool outExists() const;

  std::string firstName_;
  std::string secondName_;
  //std::string inFormat_;
  inFormats inFormat_ = inFormats::NOFORMAT;
  outFormats outFormat_ = outFormats::NOFORMAT;

  static inFormats getInFormat(const std::string & format);
  static outFormats getOutFormat(const std::string & format);
  static outFormats getOutFormat(inFormats format);

  static std::string getInFormat(inFormats format);
  static std::string getOutFormat(outFormats format);
  static std::string getOutExtension(outFormats format);

  std::string getOutExtension() const;

  OutOptions out_;

  bool complementMate_ = false;

  bool processed_ = false;
  std::string lowerCaseBases_;
  bool removeGaps_ = false;
  bool includeWhiteSpaceInName_ = true;
  int32_t extra_ = 0;



  static SeqIOOptions genFastqIn(const std::string & inFilename, bool processed = false);
  static SeqIOOptions genFastqOut(const std::string & outFilename);
  static SeqIOOptions genFastqInOut(const std::string & inFilename,
  		const std::string & outFilename, bool processed = false);

  static SeqIOOptions genFastaIn(const std::string & inFilename, bool processed = false);
  static SeqIOOptions genFastaOut(const std::string & outFilename);
  static SeqIOOptions genFastaInOut(const std::string & inFilename,
  		const std::string & outFilename, bool processed = false);

  /**@b output options as json
   *
   * @return options represented in json
   */
  Json::Value toJson()const;
};



}  // namespace bibseq


