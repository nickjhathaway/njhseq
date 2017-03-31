#pragma once

/*
 * BioCmdsUtils.hpp
 *
 *  Created on: Mar 13, 2017
 *      Author: nick
 */

#include <bibcpp/system.h>
#include "bibseq/utils.h"
#include "bibseq/IO/SeqIO/SeqIOOptions.hpp"

namespace bibseq {


class BioCmdsUtils {
public:
	BioCmdsUtils();
	BioCmdsUtils(bool verbose);
	bool verbose_ = false;

	bib::sys::RunOutput RunBowtie2Index(const bfs::path & genomeFnp);
	bib::sys::RunOutput RunBwaIndex(const bfs::path & genomeFnp);
	bib::sys::RunOutput RunSamtoolsFastaIndex(const bfs::path & genomeFnp);
	bib::sys::RunOutput RunPicardFastaSeqDict(const bfs::path & genomeFnp);
	bib::sys::RunOutput RunFaToTwoBit(const bfs::path & genomeFnp);


	bib::sys::RunOutput bowtie2Align(const SeqIOOptions & opts,
			const bfs::path & genomeFnp, std::string additionalBowtie2Args =
					"");

	bib::sys::RunOutput runCmdCheck(const std::string & cmd,
			const bfs::path & input, const bfs::path & check);

	static void checkRunOutThrow(const bib::sys::RunOutput & runOut,
			const std::string & funcName);

};

} /* namespace bibseq */
