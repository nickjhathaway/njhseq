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
//  readObjectIOOpt.hpp
//
//  Created by Nick Hathaway on 06/02/15.
//


#include "njhseq/IO/SeqIO/SeqIOOptions.hpp"
#include "njhseq/IO/SeqIO/SeqInput.hpp"
#include "njhseq/IO/SeqIO/SeqOutput.hpp"


namespace njhseq {

class SeqIO {
public:

	SeqIO(const SeqIOOptions & options);

	SeqIOOptions ioOptions_;

	SeqInput in_;
	SeqOutput out_;

	template<typename T>
	bool readNextRead(T & read){
		return in_.readNextRead(read);
	}

	void openIn();
	void closeIn();
	void openOut();
	void closeOut();

	template<typename T>
	void write(T & read){
		out_.write(read);
	}

	template<typename T>
	void openWrite(T & read){
		out_.openWrite(read);
	}

	static std::unordered_map<uint32_t, std::string>
	rewriteSeqsWithIndexAsName(const SeqIOOptions &inOpts, const SeqIOOptions &outOpts);

	static std::unordered_map<uint32_t, std::string>
	rewriteSeqsWithIndexAsName(const SeqIOOptions &inOpts, const SeqIOOptions &outOpts,
														 const OutOptions &seqKeyFileOutOpts);

	static std::unordered_map<uint32_t, std::string>
	rewriteSeqsWithIndexAsName(const std::vector<seqInfo> & seqs, const SeqIOOptions &outOpts);

	static std::unordered_map<uint32_t, std::string>
	rewriteSeqsWithIndexAsName(const std::vector<seqInfo> & seqs, const SeqIOOptions &outOpts,
														 const OutOptions &seqKeyFileOutOpts);

	std::mutex mut_;

	friend class MultiSeqIO;
};

}  // namespace njhseq

