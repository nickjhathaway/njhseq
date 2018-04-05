#pragma once

/*
 * BioCmdsUtils.hpp
 *
 *  Created on: Mar 13, 2017
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include <bibcpp/system.h>
#include "bibseq/utils.h"
#include "bibseq/IO/SeqIO/SeqIOOptions.hpp"

namespace bibseq {


class BioCmdsUtils {
public:
	BioCmdsUtils();
	BioCmdsUtils(bool verbose);
	bool verbose_ = false;
	std::string fastqDumpCmd_ = "fastq-dump";

	bib::sys::RunOutput RunBowtie2Index(const bfs::path & genomeFnp) const;
	bib::sys::RunOutput RunBwaIndex(const bfs::path & genomeFnp) const;
	bib::sys::RunOutput RunSamtoolsFastaIndex(const bfs::path & genomeFnp) const;
	bib::sys::RunOutput RunPicardFastaSeqDict(const bfs::path & genomeFnp) const;
	bib::sys::RunOutput RunFaToTwoBit(const bfs::path & genomeFnp) const;
	std::unordered_map<std::string, bib::sys::RunOutput> runAllPossibleIndexes(
			const bfs::path & genomeFnp) const;

	bib::sys::RunOutput bowtie2Align(const SeqIOOptions & opts,
			const bfs::path & genomeFnp, std::string additionalBowtie2Args = "") const;

	bib::sys::RunOutput bowtie2AlignNoSort(const SeqIOOptions & opts,
			const bfs::path & genomeFnp, std::string additionalBowtie2Args = "") const;


	struct LastZPars{
		double coverage = 90;
		double identity = 50;
		bfs::path genomeFnp = "";
		std::string outFormat = "SAM";
		std::string extraLastzArgs = "";
	};
	bib::sys::RunOutput lastzAlign(const SeqIOOptions & opts, const LastZPars & pars) const ;


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
		bib::sys::RunOutput output_;

		Json::Value toJson() const;
	};

	FastqDumpResults runFastqDump(const FastqDumpPars & pars) const;

	bib::sys::RunOutput runCmdCheck(const std::string & cmd,
			const bfs::path & input, const bfs::path & check) const;

	static void checkRunOutThrow(const bib::sys::RunOutput & runOut,
			const std::string & funcName);

};

} /* namespace bibseq */

