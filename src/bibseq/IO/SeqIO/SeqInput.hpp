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
 * SeqInput.hpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */

#include <api/BamReader.h>

#include <njhcpp/files/fileObjects.h>

#include "njhseq/IO/cachedReader.hpp"
#include "njhseq/IO/InputStream.hpp"
#include "njhseq/objects/seqObjects/readObject.hpp"
#include "njhseq/objects/seqObjects/Paired/PairedRead.hpp"
#include "njhseq/objects/seqObjects/sffObject.hpp"
#include "njhseq/IO/SeqIO/SeqIOOptions.hpp"
#include "njhseq/readVectorManipulation/readVectorOperations.h"


namespace njhseq {

/**@brief Class that handles the reading of several different sequence type files
 *
 * @todo Make a class for each type of reader to simplify things and to avoid having to check the input type for each readNextRead call
 */


class SeqInput{
	std::function<bool(seqInfo &)> firstTimeReaderFunc_;
	std::function<bool(seqInfo &)> readerFunc_;
public:

	SeqInput(const SeqIOOptions & options);
	SeqInput(const SeqInput& that);

	SeqIOOptions ioOptions_;

	void setReaderFunc();

	bool readNextRead(seqInfo & read);
	bool readNextRead(PairedRead & read);



	template<typename T>
	bool readNextRead(T & read) {
		return readNextRead(read.seqBase_);
	}

	template<typename T>
	bool readNextRead(std::shared_ptr<T> & read) {
		return readNextRead(*read);
	}
	template<typename T>
	bool readNextRead(std::unique_ptr<T> & read) {
		return readNextRead(*read);
	}

	template<typename T>
	std::shared_ptr<T> readNextRead() {
		auto ret = std::make_shared<T>();
		auto status = readNextRead(ret);
		if (status) {
			return ret;
		} else {
			return nullptr;
		}
	}

	template<typename T>
	std::vector<T> readAllReads(){
		if(!inOpen_){
			openIn();
		}
		T read;
		std::vector<T> ret;
		while(readNextRead(read)){
			ret.emplace_back(read);
		}
	  readVec::handelLowerCaseBases(ret, ioOptions_.lowerCaseBases_);
	  if (ioOptions_.removeGaps_) {
	    readVec::removeGapsFromReads(ret);
	  }
		return ret;
	}

	template<typename T>
	std::vector<std::shared_ptr<T>> readAllReadsPtrs(){
		if(!inOpen_){
			openIn();
		}
		std::vector<std::shared_ptr<T>> ret;
		auto read = readNextRead<T>();
		while(read){
			ret.emplace_back(read);
			read = readNextRead<T>();
		}
	  readVec::handelLowerCaseBases(ret, ioOptions_.lowerCaseBases_);
	  if (ioOptions_.removeGaps_) {
	    readVec::removeGapsFromReads(ret);
	  }
		return ret;
	}

	size_t tellgPri();
	void seekgPri(size_t pos);

	size_t tellgSec();
	void seekgSec(size_t pos);

	bool inOpen() const;
	void openIn();
	void closeIn();
	void reOpenIn();

	std::mutex mut_;
	std::unique_ptr<sffObject> lastSffRead_;

	template<typename T>
	static std::vector<T> getSeqVec(const SeqIOOptions & seqOptions, uint64_t& maxLength){
		std::vector<T> ret;
		SeqInput reader(seqOptions);
		reader.openIn();
		T seq;
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxLength);
			ret.emplace_back(seq);
		}
		return ret;
	}

	template<typename T>
	static std::vector<T> getSeqVec(const SeqIOOptions & seqOptions){
		uint64_t maxlen = 0;
		return getSeqVec<T>(seqOptions, maxlen);
	}

	static std::vector<readObject> getReferenceSeq(
			const SeqIOOptions & refOptions, uint64_t& maxLength);

	bool loadIndex();


	template<typename T>
	void getReadNoCheck(T & read, size_t pos){
		seekgPri(index_[pos]);
		readNextRead(read);
	}

	template<typename T>
	void getRead(T & read, size_t pos) {
		loadIndex();
		seekToSeqIndex(pos);
		readNextRead(read);
	}

	template<typename T>
	std::vector<T> getReads(size_t pos, size_t number) {
		openIn();
		loadIndex();
		seekToSeqIndex(pos);
		std::vector<T> ret;
		uint32_t count = 0;
		T read;
		while (readNextRead(read) && count < number) {
			ret.emplace_back(read);
			++count;
		}
		return ret;
	}

	template<typename T>
	std::vector<std::shared_ptr<T>> getReadsPtrs(size_t pos, size_t number) {
		openIn();
		loadIndex();
		seekToSeqIndex(pos);
		std::vector<std::shared_ptr<T>> ret;
		uint32_t count = 0;
		std::shared_ptr<T> read = std::make_shared<T>();
		while (readNextRead(read) && count < number) {
			ret.emplace_back(read);
			++count;
		}
		return ret;
	}

	void seekToSeqIndex(size_t seqIndex);
	void seekToSeqIndexNoCheck(size_t seqIndex);
	void seekToSeqIndexNoCheckNoPairedCheck(size_t seqIndex);

	static void buildIndex(const SeqIOOptions & ioOpts);


	std::vector<unsigned long long> randomlySampleIndex(
			njh::randomGenerator & gen, const std::string& sample) const;

	bool isFirstEmpty() const;

	//will throw if second is blank
	bool isSecondEmpty() const;

private:

	std::vector<unsigned long long> index_;
	//only created for files with secondary files fasta  with qual and R1 and R2 paired reads
	std::vector<unsigned long long> secIndex_;
	bool indexLoad_ = false;

	std::unique_ptr<BamTools::BamReader> bReader_;
	std::unique_ptr<BamTools::BamAlignment> aln_;


	std::unique_ptr<InputStream> priReader_;
	std::unique_ptr<InputStream> secReader_;


public:
	std::unique_ptr<VecStr> sffTxtHeader_;
	std::unique_ptr<sffBinaryHeader> sffBinHeader_;
	bool readNextFastaStream(std::istream& fastaFile, seqInfo& read,
			bool processed);
	bool readNextFastqStream(std::istream& fastqFile, uint32_t offSet, seqInfo& read,
			bool processed);
private:
	bool readNextQualStream(std::istream& qualFile, std::vector<uint32_t>& quals,
			std::string & name);
	bool readNextFastaQualStream(std::istream& fastaFile,
			std::istream& qualFile, seqInfo & read, bool processed);

	bool readNextFastqStream(const VecStr & data, const uint32_t lCount, uint32_t offSet, seqInfo& read,
			bool processed);

	bool readNextBam(BamTools::BamReader & bReader, seqInfo& read,
			BamTools::BamAlignment & aln, bool processed);

	VecStr readSffTxtHeader(std::istream & inFile);
	bool readNextSff(std::istream & inFile, sffObject & read);
	bool readNextSff(std::istream & inFile, seqInfo & read);

	bool readNextSffBin(std::istream& in, sffObject& read, int32_t numFlowReads);
	void readHeader(std::istream& in, sffBinaryHeader& header);

	bool inOpen_ = false;

	friend class SeqIO;
};

template<typename T>
std::vector<T> getSeqs(const bfs::path & fnp) {
	SeqIOOptions popOpts;
	popOpts.firstName_ = fnp;
	popOpts.inFormat_ = SeqIOOptions::getInFormat(njh::files::getExtension(fnp));
	SeqInput reader(popOpts);
	return reader.readAllReads<T>();
}


}  // namespace njhseq


