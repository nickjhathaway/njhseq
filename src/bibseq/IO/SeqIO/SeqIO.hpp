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
//  readObjectIOOpt.hpp
//
//  Created by Nick Hathaway on 06/02/15.
//  Copyright (c) 2012 Nick Hathaway. All rights reserved.
//


#include "bibseq/IO/SeqIO/SeqIOOptions.hpp"
#include "bibseq/IO/SeqIO/SeqInput.hpp"
#include "bibseq/IO/SeqIO/SeqOutput.hpp"


namespace bibseq {

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


	std::mutex mut_;

	friend class MultiSeqIO;
};

}  // namespace bibseq

