#pragma once
//
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
//
//  readObjectIOOptions.hpp
//
//  Created by Nick Hathaway on 2/07/15.
//


#include <njhcpp/jsonUtils.h>
#include "njhseq/IO/IOUtils.hpp"

namespace njhseq {


struct SeqIOOptions {
  enum class inFormats {
  	FASTQ,
		FASTQPAIRED,
		FASTQPAIREDGZ,
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
		FASTQGZ,
		FASTA,
		FASTAGZ,
		FASTQPAIRED,
		FASTQPAIREDGZ,
		FASTAQUAL,
		FLOW,
		FLOWMOTHUR,
		NOFORMAT
  };

  SeqIOOptions();

  SeqIOOptions(const bfs::path & firstName,
  		inFormats inFormat, bool processed);

	SeqIOOptions(const OutOptions & out, outFormats outFormat);

	SeqIOOptions(const bfs::path & outFilename, outFormats outFormat);

	SeqIOOptions(const bfs::path & outFilename, outFormats outFormat, const OutOptions & out);

	/**@b Create from a json string
   *
   * @param jsonStr The value stored in json
   */
  explicit SeqIOOptions(const std::string & jsonStr);


  bool inExists() const;
	bool outExists() const;

  bfs::path firstName_;
  bfs::path secondName_;
  //std::string inFormat_;
  inFormats inFormat_ = inFormats::NOFORMAT;
  outFormats outFormat_ = outFormats::NOFORMAT;

  static inFormats getInFormat(const std::string & format);
  static outFormats getOutFormat(const std::string & format);
  static inFormats getInFormatFromFnp(const bfs::path & fnp);
  static outFormats getOutFormatFromFnp(const bfs::path & fnp);
  static outFormats getOutFormat(inFormats format);
  static inFormats getInFormat(outFormats format);

  static std::string getInFormat(inFormats format);
  static std::string getOutFormat(outFormats format);
  static std::string getOutExtension(outFormats format);
  static std::string getOutExtensionSecondary(outFormats format);

  std::string getOutExtension() const;

  OutOptions out_;

  bfs::path getPriamryOutName() const;
  bfs::path getSecondaryOutName() const;

  bool revComplMate_ = false;

  bool processed_ = false;
  std::string lowerCaseBases_;
  bool removeGaps_ = false;
  bool includeWhiteSpaceInName_ = true;
  int32_t extra_ = 0;

	bool isPairedIn() const;
	bool isPairedOut() const;

  //fastq
  static SeqIOOptions genFastqIn(const bfs::path & inFilename, bool processed =false);
	static SeqIOOptions genFastqOut(const bfs::path & outFilename);
	static SeqIOOptions genFastqInOut(const bfs::path & inFilename,const bfs::path & outFilename, bool processed = false);

	//gz fastq
  static SeqIOOptions genFastqInGz(const bfs::path & inFilename, bool processed =false);
	static SeqIOOptions genFastqOutGz(const bfs::path & outFilename);
	static SeqIOOptions genFastqInOutGz(const bfs::path & inFilename,const bfs::path & outFilename, bool processed = false);

	//fasta
	static SeqIOOptions genFastaIn(const bfs::path & inFilename, bool processed =false);
	static SeqIOOptions genFastaOut(const bfs::path & outFilename);
	static SeqIOOptions genFastaInOut(const bfs::path & inFilename,const bfs::path & outFilename, bool processed = false);

	//gz fasta
	static SeqIOOptions genFastaInGz(const bfs::path & inFilename, bool processed =false);
	static SeqIOOptions genFastaOutGz(const bfs::path & outFilename);
	static SeqIOOptions genFastaInOutGz(const bfs::path & inFilename, const bfs::path & outFilename, bool processed = false);

	//paired
	static SeqIOOptions genPairedOut(const bfs::path & outFilename);
	static SeqIOOptions genPairedIn(const bfs::path & r1reads, const bfs::path & r2reads, bool processed = false);
	static SeqIOOptions genPairedInOut(const bfs::path & r1reads, const bfs::path & r2reads, const bfs::path & outFilename, bool processed = false);

	//gz paired
	static SeqIOOptions genPairedOutGz(const bfs::path & outFilename);
	static SeqIOOptions genPairedInGz(const bfs::path & r1reads, const bfs::path & r2reads, bool processed = false);
	static SeqIOOptions genPairedInOutGz(const bfs::path & r1reads, const bfs::path & r2reads, const bfs::path & outFilename, bool processed = false);


	/**@b output options as json
	 *
	 * @return options represented in json
	 */
	Json::Value toJson() const;
};



}  // namespace njhseq


