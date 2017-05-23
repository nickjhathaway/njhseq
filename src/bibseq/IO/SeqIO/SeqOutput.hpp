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
 * SeqOutput.hpp
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


class SeqOutput{
public:
	SeqOutput(const SeqIOOptions & options);


	SeqIOOptions ioOptions_;

	bool outOpen() const;
	void openOut();
	void closeOut();
	void closeOutForReopening();

	template<typename T>
	void write(const T & read) {
		write(read.seqBase_);
	}

	template<typename T>
	void openWrite(const T & read) {
		openWrite(read.seqBase_);
	}

	template<typename T>
	void write(const std::unique_ptr<T> & read) {
		write(*read);
	}

	template<typename T>
	void openWrite(const std::unique_ptr<T> & read) {
		openWrite(*read);
	}

	template<typename T>
	void write(const std::shared_ptr<T> & read) {
		write(*read);
	}

	template<typename T>
	void openWrite(const std::shared_ptr<T> & read) {
		openWrite(*read);
	}

	template<typename T>
	void writeNoCheck(const T & read) {
		writeNoCheck(read.seqBase_);
	}

	template<typename T>
	void writeNoCheck(const std::unique_ptr<T> & read) {
		writeNoCheck(*read);
	}

	template<typename T>
	void writeNoCheck(const std::shared_ptr<T> & read) {
		writeNoCheck(*read);
	}

	template<typename T>
	void write(const std::vector<T> & reads) {
		if (!outOpen_) {
			throw std::runtime_error {
					"Error in readObjectIOOpt, attempted to write when out files aren't open, out file: " + ioOptions_.out_.outName().string() };
		}
		for (const auto & read : reads) {
			writeNoCheck(read);
		}
	}

	template<typename T>
	void openWrite(const std::vector<T> & reads) {
		if (!outOpen_) {
			openOut();
		}
		for (const auto & read : reads) {
			writeNoCheck(read);
		}
	}

	template<typename T>
	void writeNoCheck(const std::vector<T> & reads) {
		for (const auto & read : reads) {
			writeNoCheck(read);
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

	void writeNoCheck(const seqInfo & read);
	void write(const seqInfo & read);
	void openWrite(const seqInfo & read);

	void openWriteFlow(const sffObject & read);
	void writeNoCheckFlow(const sffObject & read);

	void writeNoCheck(const std::string & line);
	void openWrite(const std::string & line);

	void writeNoCheck(const PairedRead & read);
	void openWrite(const PairedRead & read);
	void write(const PairedRead & read);

	size_t tellpPri();
	void seekpPri(size_t pos);

	size_t tellpSec();
	void seekpSec(size_t pos);

	std::mutex mut_;
private:
	bool outOpen_ = false;
	std::unique_ptr<std::ofstream> primaryOut_;
	std::unique_ptr<std::ofstream> secondaryOut_;

	std::unique_ptr<BamTools::BamWriter> bReader_;
	std::unique_ptr<BamTools::BamAlignment> aln_;

	friend class SeqIO;
};

}  // namespace bibseq



