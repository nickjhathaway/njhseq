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
/*
 * SeqOutput.hpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "njhseq/IO/cachedReader.hpp"
#include "njhseq/objects/seqObjects/readObject.hpp"
#include "njhseq/objects/seqObjects/Paired/PairedRead.hpp"
#include "njhseq/objects/seqObjects/sffObject.hpp"
#include "njhseq/IO/SeqIO/SeqIOOptions.hpp"
#include "njhseq/IO/OutputStream.hpp"

namespace njhseq {


class SeqOutput{
	std::function<bool(seqInfo &)> writerFunc_;
public:
	SeqOutput(const SeqIOOptions & options);


	SeqIOOptions ioOptions_;

	bool outOpen() const;
	void openOut();
	void closeOut();
	void closeOutForReopening();

	template<typename T>
	void write(const T & seq) {
		write(seq.seqBase_);
	}

	template<typename T>
	void openWrite(const T & seq) {
		openWrite(seq.seqBase_);
	}

	template<typename T>
	void write(const std::unique_ptr<T> & seq) {
		write(*seq);
	}

	template<typename T>
	void openWrite(const std::unique_ptr<T> & seq) {
		openWrite(*seq);
	}

	template<typename T>
	void write(const std::shared_ptr<T> & seq) {
		write(*seq);
	}

	template<typename T>
	void openWrite(const std::shared_ptr<T> & seq) {
		openWrite(*seq);
	}

	template<typename T>
	void writeNoCheck(const T & seq) {
		writeNoCheck(seq.seqBase_);
	}

	template<typename T>
	void writeNoCheck(const std::unique_ptr<T> & seq) {
		writeNoCheck(*seq);
	}

	template<typename T>
	void writeNoCheck(const std::shared_ptr<T> & seq) {
		writeNoCheck(*seq);
	}

	template<typename T>
	void write(const std::vector<T> & reads) {
		if (!outOpen_) {
			throw std::runtime_error {
					"Error in readObjectIOOpt, attempted to write when out files aren't open, out file: " + ioOptions_.out_.outName().string() };
		}
		for (const auto & seq : reads) {
			writeNoCheck(seq);
		}
	}

	template<typename T>
	void openWrite(const std::vector<T> & reads) {
		if (!outOpen_) {
			openOut();
		}
		for (const auto & seq : reads) {
			writeNoCheck(seq);
		}
	}

	template<typename T>
	void writeNoCheck(const std::vector<T> & reads) {
		for (const auto & seq : reads) {
			writeNoCheck(seq);
		}
	}

	template<typename T>
	static void write(const std::vector<T> & reads, const SeqIOOptions & opts){
		SeqOutput writer(opts);
		writer.openWrite(reads);
	}

	template<typename T>
	static void write(const std::vector<T> & reads, const std::string & filename, const SeqIOOptions & opts){
		SeqIOOptions outOpts(filename,opts.outFormat_, opts.out_);
		write(reads, outOpts);
	}

	void writeNoCheck(const seqInfo & seq);
	void write(const seqInfo & seq);
	void openWrite(const seqInfo & seq);

	void openWriteFlow(const sffObject & seq);
	void writeNoCheckFlow(const sffObject & seq);

	void writeNoCheck(const std::string & line);
	void openWrite(const std::string & line);

	void writeNoCheck(const PairedRead & seq);
	void openWrite(const PairedRead & seq);
	void write(const PairedRead & seq);

	size_t tellpPri();
	void seekpPri(size_t pos);

	size_t tellpSec();
	void seekpSec(size_t pos);

	std::mutex mut_;

	bfs::path getPrimaryOutFnp() const;
	bfs::path getSecondaryOutFnp() const;
private:
	bool outOpen_ = false;

	std::unique_ptr<OutputStream> primaryOut_;
	std::unique_ptr<OutputStream> secondaryOut_;

	std::unique_ptr<BamTools::BamWriter> bReader_;
	std::unique_ptr<BamTools::BamAlignment> aln_;

	friend class SeqIO;
};

}  // namespace njhseq



