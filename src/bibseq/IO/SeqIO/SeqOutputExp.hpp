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
/*
 * SeqOutputExp.hpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "bibseq/IO/cachedReader.hpp"
#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/seqObjects/Paired/PairedRead.hpp"
#include "bibseq/objects/seqObjects/sffObject.hpp"
#include "bibseq/IO/SeqIO/SeqIOOptions.hpp"


namespace bibseq {


class SeqOutputExp{
	std::function<bool(seqInfo &)> writerFunc_;
public:
	SeqOutputExp(const SeqIOOptions & options);


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
		SeqOutputExp writer(opts);
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
private:
	bool outOpen_ = false;

//	std::unique_ptr<std::ostream> primaryOut_;
//	std::unique_ptr<std::ostream> secondaryOut_;
//
//	std::unique_ptr<std::ofstream> primaryFileOut_;
//	std::unique_ptr<std::ofstream> secondaryFileOut_;

		std::unique_ptr<std::ofstream> primaryOut_;
		std::unique_ptr<std::ofstream> secondaryOut_;

	std::unique_ptr<bib::GZSTREAM::ogzstream> primaryGzOut_;
	std::unique_ptr<bib::GZSTREAM::ogzstream> secondaryGzOut_;

	std::unique_ptr<BamTools::BamWriter> bReader_;
	std::unique_ptr<BamTools::BamAlignment> aln_;

	friend class SeqIO;
};

}  // namespace bibseq



