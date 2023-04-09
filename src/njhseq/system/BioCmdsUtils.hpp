#pragma once

/*
 * BioCmdsUtils.hpp
 *
 *  Created on: Mar 13, 2017
 *      Author: nick
 */
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
#include <njhcpp/system.h>
#include "njhseq/utils.h"
#include "njhseq/IO/SeqIO/SeqIOOptions.hpp"

namespace njhseq {


class BioCmdsUtils {
public:
	BioCmdsUtils();
	BioCmdsUtils(bool verbose);
	bool verbose_ = false;
	std::string fastqDumpCmd_ = "fastq-dump";
	std::string fasterqDumpCmd_ = "fasterq-dump";

  njh::sys::RunOutput RunBowtie2Index(const bfs::path & genomeFnp) const;
  njh::sys::RunOutput RunMakeblastdb(const bfs::path & genomeFnp) const;
	njh::sys::RunOutput RunBwaIndex(const bfs::path & genomeFnp) const;
	njh::sys::RunOutput RunSamtoolsFastaIndex(const bfs::path & genomeFnp) const;
	njh::sys::RunOutput RunPicardFastaSeqDict(const bfs::path & genomeFnp) const;
	njh::sys::RunOutput RunFaToTwoBit(const bfs::path & genomeFnp) const;
	std::unordered_map<std::string, njh::sys::RunOutput> runAllPossibleIndexes(
			const bfs::path & genomeFnp) const;

	njh::sys::RunOutput bowtie2Align(const SeqIOOptions & opts,
			const bfs::path & genomeFnp, std::string additionalBowtie2Args = "") const;

	njh::sys::RunOutput bowtie2AlignNoSort(const SeqIOOptions & opts,
			const bfs::path & genomeFnp, std::string additionalBowtie2Args = "") const;


	struct LastZPars{
		double coverage = 90;
		double identity = 50;
		bfs::path genomeFnp = "";
		std::string outFormat = "SAM";
		std::string extraLastzArgs = "";
	};
	njh::sys::RunOutput lastzAlign(const SeqIOOptions & opts, const LastZPars & pars) const ;
	njh::sys::RunOutput lastzAlignNoSort(const SeqIOOptions & opts, const LastZPars & pars) const ;


	bool isSRAPairedEnd(const bfs::path & sraFnp) const;

	struct FastqDumpPars{
		bfs::path sraFnp_;

		std::string extraSraOptions_ = "";
		bfs::path outputDir_ = "./";
		bool exportBarCode_ = false;
		bool gzip_ = false;
		bool force_ = false;
	};

	struct FastqDumpResults{
		bfs::path firstMateFnp_;
		bfs::path barcodeFnp_;
		bfs::path secondMateFnp_;



		bool isGzipped_ = false;
		bool isPairedEnd_ = false;
		njh::sys::RunOutput output_;

		Json::Value toJson() const;
	};



	struct FasterqDumpResults{
		bfs::path firstMateFnp_;
		bfs::path secondMateFnp_;

		bool isPairedEnd_ = false;
		njh::sys::RunOutput output_;

		Json::Value toJson() const;
	};

	struct FasterqDumpPars{
		bfs::path sraFnp_;

		std::string extraSraOptions_ = "";
		bfs::path tempDir_ = "./";
		bfs::path outputDir_ = "./";
		bool force_ = false;

		uint32_t numThreads_ = 1;

		uint32_t readSubCacheSize_{100};
		uint32_t writeCacheSize{50};
		std::chrono::milliseconds writeWaitTime{std::chrono::milliseconds(1)};

	};

	FastqDumpResults runFastqDump(const FastqDumpPars & pars) const;
	FasterqDumpResults runFasterqDump(const FasterqDumpPars & pars) const;

	njh::sys::RunOutput runCmdCheck(const std::string & cmd,
			const bfs::path & input, const bfs::path & check) const;

	static void checkRunOutThrow(const njh::sys::RunOutput & runOut,
			const std::string & funcName);

};

} /* namespace njhseq */

