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

#include "SeqInput.hpp"
#include "njhseq/seqToolsUtils/seqToolsUtils.hpp"
#include <TwoBit.h>

#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"


namespace njhseq {


void SeqInput::seekToSeqIndex(size_t pos){
	if (pos >= index_.size()) {
		std::stringstream ss;
		ss << "In " << __PRETTY_FUNCTION__ << " pos: " << pos
				<< " is greater than max: " << index_.size() << std::endl;
		throw std::out_of_range(ss.str());
	}
	seekgPri(index_[pos]);
	if(ioOptions_.isPairedIn()){
		seekgSec(secIndex_[pos]);
	}
}

void SeqInput::seekToSeqIndexNoCheck(size_t pos){
	seekgPri(index_[pos]);
	if(ioOptions_.isPairedIn()){
		seekgSec(secIndex_[pos]);
	}
}

void SeqInput::seekToSeqIndexNoCheckNoPairedCheck(size_t pos){
	seekgPri(index_[pos]);
	seekgSec(secIndex_[pos]);
}



void SeqInput::buildIndex(const SeqIOOptions & ioOpts) {
	if(ioOpts.isPairedIn()){
		SeqInput reader(ioOpts);
		bfs::path outIndexName = ioOpts.firstName_.string() + ".idx.gz";
		bfs::path outSecIndexName = ioOpts.secondName_.string() + ".idx.gz";
		PairedRead info;
		reader.openIn();
		unsigned long long pos = 0;
		unsigned long long secpos = 0;
		std::vector<unsigned long long> index;
		std::vector<unsigned long long> secIndex;
		while (reader.readNextRead(info)) {
			index.emplace_back(pos);
			pos = reader.tellgPri();

			secIndex.emplace_back(secpos);
			secpos = reader.tellgSec();
		}
		njh::files::writePODvectorGz(outIndexName, index);
		njh::files::writePODvectorGz(outSecIndexName, secIndex);
	}else{
		SeqInput reader(ioOpts);
		bfs::path outIndexName = ioOpts.firstName_.string() + ".idx.gz";
		seqInfo info;
		reader.openIn();
		unsigned long long pos = 0;
		std::vector<unsigned long long> index;
		while (reader.readNextRead(info)) {
			index.emplace_back(pos);
			pos = reader.tellgPri();
		}
		njh::files::writePODvectorGz(outIndexName, index);
	}

}

bool SeqInput::loadIndex() {
	if(ioOptions_.isPairedIn()){
		bfs::path outIndexName = ioOptions_.firstName_.string() + ".idx.gz";
		bfs::path secOutIndexName = ioOptions_.secondName_.string() + ".idx.gz";
		bool indexNeedsUpdate = true;
		if (njh::files::bfs::exists (outIndexName)) {
			auto indexTime = njh::files::last_write_time(outIndexName);
			auto fileTime = njh::files::last_write_time(
					ioOptions_.firstName_);
			if (fileTime < indexTime) {
				indexNeedsUpdate = false;
			}
		}
		if (njh::files::bfs::exists (secOutIndexName)) {
			auto indexTime = njh::files::last_write_time(secOutIndexName);
			auto fileTime = njh::files::last_write_time(
					ioOptions_.secondName_);
			if (fileTime < indexTime) {
				indexNeedsUpdate = false;
			}
		}else{
			indexNeedsUpdate = true;
		}
		if (indexNeedsUpdate) {
			buildIndex(ioOptions_);
			index_ = njh::files::readPODvectorGz<unsigned long long>(outIndexName);
			secIndex_ = njh::files::readPODvectorGz<unsigned long long>(secOutIndexName);
			indexLoad_ = true;
		} else {
			if(!indexLoad_){
				index_ = njh::files::readPODvectorGz<unsigned long long>(outIndexName);
				secIndex_ = njh::files::readPODvectorGz<unsigned long long>(secOutIndexName);
				indexLoad_ = true;
			}
		}
		return indexNeedsUpdate;
	}else{
		bfs::path outIndexName = ioOptions_.firstName_.string() + ".idx.gz";
		bool indexNeedsUpdate = true;
		if (njh::files::bfs::exists (outIndexName)) {
			auto indexTime = njh::files::last_write_time(outIndexName);
			auto fileTime = njh::files::last_write_time(
					ioOptions_.firstName_);
			if (fileTime < indexTime) {
				indexNeedsUpdate = false;
			}
		}
		if (indexNeedsUpdate) {
			buildIndex(ioOptions_);
			index_ = njh::files::readPODvectorGz<unsigned long long>(outIndexName);
			indexLoad_ = true;
		} else {
			if(!indexLoad_){
				index_ = njh::files::readPODvectorGz<unsigned long long>(outIndexName);
				indexLoad_ = true;
			}
		}
		return indexNeedsUpdate;
	}

}

std::vector<readObject> SeqInput::getReferenceSeq(
    const SeqIOOptions & refOptions, uint64_t& maxLength) {
	SeqInput refReader(refOptions);
	std::vector<readObject> reads = refReader.readAllReads<readObject>();
  readVec::getMaxLength(reads, maxLength);
  return reads;
}


SeqInput::SeqInput(const SeqIOOptions & options) :
				ioOptions_(options) {

}

SeqInput::SeqInput(const SeqInput& that){
	ioOptions_ = that.ioOptions_;
}

bool SeqInput::inOpen() {
	std::lock_guard<std::mutex> lock(mut_);
	return inOpenLockFree();
}

bool SeqInput::inOpenLockFree() const {
	return inOpen_;
}


size_t SeqInput::tellgPri(){
	std::lock_guard<std::mutex> lock(mut_);
	return priReader_->tellg();
}
void SeqInput::seekgPri(size_t pos){
	std::lock_guard<std::mutex> lock(mut_);
	priReader_->seekg(pos);
}

size_t SeqInput::tellgSec(){
	std::lock_guard<std::mutex> lock(mut_);
	return secReader_->tellg();
}
void SeqInput::seekgSec(size_t pos){
	std::lock_guard<std::mutex> lock(mut_);
	secReader_->seekg(pos);
}


std::vector<unsigned long long> SeqInput::randomlySampleIndex(
		njh::randomGenerator & gen, const std::string& sample) const {
	uint32_t sampleNum = processRunCutoff(sample, index_.size());
	return gen.unifRandSelectionVec(index_, sampleNum, false);
}

bool SeqInput::isFirstEmpty() const {
	if (0 == bfs::file_size(ioOptions_.firstName_)) {
		return true;
	}
	return false;
}

//will throw if second is blank
bool SeqInput::isSecondEmpty() const{
	if (0 == bfs::file_size(ioOptions_.secondName_)) {
		return true;
	}
	return false;
}


void SeqInput::setReaderFunc(){
	std::stringstream ssFormatCheck;
	switch (ioOptions_.inFormat_) {
		case SeqIOOptions::inFormats::FASTQ:
		case SeqIOOptions::inFormats::FASTQGZ:
			readerFunc_ = [this](seqInfo & seq){
				return readNextFastqStream(*priReader_, SangerQualOffset, seq,
						ioOptions_.processed_);
			};
			break;
		case SeqIOOptions::inFormats::FASTA:
		case SeqIOOptions::inFormats::FASTAGZ:
			readerFunc_ = [this](seqInfo & seq) {
				return readNextFastaStream(*priReader_, seq, ioOptions_.processed_);
			};
			break;
		case SeqIOOptions::inFormats::FASTQPAIRED:
		case SeqIOOptions::inFormats::FASTQPAIREDGZ:
			break;
		case SeqIOOptions::inFormats::FASTAQUAL:
			readerFunc_ = [this](seqInfo & seq) {
				return readNextFastaQualStream(*priReader_, *secReader_, seq,
						ioOptions_.processed_);;
			};
			break;
		case SeqIOOptions::inFormats::BAM:
			readerFunc_ = [this](seqInfo & seq) {
				return readNextBam(*bReader_, seq, *aln_, ioOptions_.processed_);
			};
			break;
		case SeqIOOptions::inFormats::SFFTXT:
			readerFunc_ = [this](seqInfo & seq) {
				return readNextSff(*priReader_, seq);
			};
			break;
		case SeqIOOptions::inFormats::SFFBIN:
			readerFunc_ = [this](seqInfo & seq) {
				bool wasAbleToRead = false;
				lastSffRead_ = std::make_unique<sffObject>();
				readNextSffBin(*priReader_, *lastSffRead_,
						sffBinHeader_->numFlowsPerRead);
				if(lastSffRead_->sanityCheck()){
					seq = lastSffRead_->seqBase_;
					wasAbleToRead = true;
				}
				return wasAbleToRead;
			};
			break;
		default:
			ssFormatCheck << "Unrecognized file type : " << " in " << __PRETTY_FUNCTION__
					<< ", not reading " << ioOptions_.firstName_ << std::endl;
			ssFormatCheck << "Acceptable types are fasta,fastaQual,fastq,fastqPaired, bam, fastagz, sff, and sffbin" << std::endl;
			throw std::runtime_error { ssFormatCheck.str() };
			break;
	}
	firstTimeReaderFunc_ = [this](seqInfo & seq) {
		if (!inOpen_) {
			throw std::runtime_error { "Error in " + std::string(__PRETTY_FUNCTION__)
					+ ", attempted to read when in file " + ioOptions_.firstName_.string()
					+ " hasn't been opened" };
		}
		bool wasAbleToRead = readerFunc_(seq);
		firstTimeReaderFunc_ = readerFunc_;
		return wasAbleToRead;
	};
}

void SeqInput::openInLockFree() {
	if (inOpen_) {
		return;
	}
	if (njh::strToLowerRet(ioOptions_.firstName_.string()) != "stdin"
			&& !njh::files::bfs::exists(ioOptions_.firstName_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error file: "
				<< njh::bashCT::boldRed(ioOptions_.firstName_.string())
				<< " doesn't exist\n";
		if ("" != ioOptions_.secondName_.string()
				&& !njh::files::bfs::exists(ioOptions_.secondName_.string())) {
			ss << "and file: "
					<< njh::bashCT::boldRed(ioOptions_.secondName_.string())
					<< " doesn't exist\n";
		}
		throw std::runtime_error { ss.str() };
	}
	bool failedToOpen = false;

	auto openPrim =
			[this,&failedToOpen]() {
				priReader_ = std::make_unique<InputStream>(InOptions(ioOptions_.firstName_));
				if (!(*priReader_)) {
					failedToOpen = true;
				}
			};
	auto openPrimSec =
			[this,&failedToOpen]() {
				priReader_ = std::make_unique<InputStream>(InOptions(ioOptions_.firstName_));
				if (!(*priReader_)) {
					failedToOpen = true;
				}
				secReader_ = std::make_unique<InputStream>(InOptions(ioOptions_.secondName_));
				if (!(*secReader_)) {
					failedToOpen = true;
				}
			};
	readerFunc_ = [](seqInfo & seq) {
		std::stringstream ss;
		ss << "Error in seqInput, readerFunc_ not set" << "\n";
		throw std::runtime_error {ss.str()};
		return false;
	};
	std::stringstream ssFormatCheck;
	switch (ioOptions_.inFormat_) {
	case SeqIOOptions::inFormats::FASTQ:
	case SeqIOOptions::inFormats::FASTQGZ:
	case SeqIOOptions::inFormats::FASTA:
	case SeqIOOptions::inFormats::FASTAGZ:
		openPrim();
		break;
	case SeqIOOptions::inFormats::FASTQPAIREDGZ:
	case SeqIOOptions::inFormats::FASTQPAIRED:
	case SeqIOOptions::inFormats::FASTAQUAL:
		openPrimSec();
		break;
	case SeqIOOptions::inFormats::BAM:
		bReader_ = std::make_unique<BamTools::BamReader>();
		aln_ = std::make_unique<BamTools::BamAlignment>();
		bReader_->Open(ioOptions_.firstName_.string());
		if (!bReader_->IsOpen()) {
			failedToOpen = true;
		}
		break;
	case SeqIOOptions::inFormats::SFFTXT:
		openPrim();
		sffTxtHeader_ = std::make_unique<VecStr>(readSffTxtHeader(*priReader_));
		break;
	case SeqIOOptions::inFormats::SFFBIN:
		//priReader_ = std::make_unique<std::ifstream>(ioOptions_.firstName_.string(),std::ios::binary);
		//			if (!(*priReader_)) {
		//				failedToOpen = true;
		//			}
		openPrim();
		sffBinHeader_ = std::make_unique<sffBinaryHeader>();
		readHeader(*priReader_, *sffBinHeader_);
		break;
	default:
		ssFormatCheck << "Unrecognized file type : " << " in "
				<< __PRETTY_FUNCTION__ << ", not reading " << ioOptions_.firstName_
				<< std::endl;
		ssFormatCheck
				<< "Acceptable types are fasta,fastaQual,fastq,fastqPaired, bam, fastagz, sff, and sffbin"
				<< std::endl;
		throw std::runtime_error { ssFormatCheck.str() };
		break;
	}

	if (failedToOpen) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error in opening : "
				<< ioOptions_.firstName_;
		if ("" != ioOptions_.secondName_) {
			ss << " or " << ioOptions_.secondName_;
		}
		ss << "\n";
		throw std::runtime_error { ss.str() };
	}
	setReaderFunc();
	inOpen_ = true;
}

void SeqInput::openIn() {
	std::lock_guard<std::mutex> lock(mut_);
	openInLockFree();
}
void SeqInput::closeInLockFree() {
	if(inOpen_){
		priReader_ = nullptr;
		secReader_ = nullptr;
		if (bReader_ && bReader_->IsOpen()) {
			bReader_->Close();
		}
		inOpen_ = false;
		firstTimeReaderFunc_ = [this](seqInfo & seq) {
			if (!inOpen_) {
				throw std::runtime_error { "Error in " + std::string(__PRETTY_FUNCTION__)
						+ ", attempted to read when in file " + ioOptions_.firstName_.string()
						+ " hasn't been opened" };
			}
			bool wasAbleToRead = readerFunc_(seq);
			firstTimeReaderFunc_ = readerFunc_;
			return wasAbleToRead;
		};
	}
}
void SeqInput::closeIn() {
	std::lock_guard<std::mutex> lock(mut_);
	closeInLockFree();
}

void SeqInput::reOpenIn(){
	std::lock_guard<std::mutex> lock(mut_);
	reOpenInLockFree();
}

void SeqInput::reOpenInLockFree(){
	closeInLockFree();
	openInLockFree();
}

bool SeqInput::readNextReadLock(seqInfo & seq) {
	std::lock_guard<std::mutex> lock(mut_);
	return readNextRead(seq);
}

bool SeqInput::readNextRead(seqInfo & seq) {
	return firstTimeReaderFunc_(seq);
}

bool SeqInput::readNextReadLock(PairedRead & seq) {
	std::lock_guard<std::mutex> lock(mut_);
	return readNextRead(seq);
}

bool SeqInput::readNextRead(PairedRead & seq) {
	if (!inOpen_) {
		throw std::runtime_error { "Error in " + std::string(__PRETTY_FUNCTION__)
				+ ", attempted to read when in file " + ioOptions_.firstName_.string()
				+ " hasn't been opened" };
	}
	if (SeqIOOptions::inFormats::FASTQPAIRED == ioOptions_.inFormat_
			|| SeqIOOptions::inFormats::FASTQPAIREDGZ == ioOptions_.inFormat_) {
		bool firstMate = readNextFastqStream(*priReader_, SangerQualOffset, seq.seqBase_,
				ioOptions_.processed_);
		bool secondMate = readNextFastqStream(*secReader_, SangerQualOffset, seq.mateSeqBase_,
				ioOptions_.processed_);
		readVec::handelLowerCaseBases(seq.seqBase_, ioOptions_.lowerCaseBases_);
		readVec::handelLowerCaseBases(seq.mateSeqBase_, ioOptions_.lowerCaseBases_);
		if(secondMate && ioOptions_.revComplMate_){
			seq.mateSeqBase_.reverseComplementRead(false, true);
		}

		seq.mateRComplemented_ = ioOptions_.revComplMate_;

		if((firstMate && secondMate) || (!firstMate && !secondMate)){
			return firstMate && secondMate;
		}else{
			std::stringstream ss;
			ss << "Error in " << __PRETTY_FUNCTION__;
			if (firstMate && !secondMate) {
				ss << " read from " << ioOptions_.firstName_ << " but not from "
						<< ioOptions_.secondName_ << "\n";
			} else {
				ss << " did not read from " << ioOptions_.firstName_
						<< " but read from " << ioOptions_.secondName_ << "\n";
			}
			throw std::runtime_error{ss.str()};
		}
	} else {
		std::stringstream ss;
		ss << "Unrecognized in format in " << __PRETTY_FUNCTION__ << ", format: " << "\n";
		throw std::runtime_error{ss.str()};
	}
	return false;
}

bool SeqInput::readNextBam(BamTools::BamReader & bReader, seqInfo& read,
		BamTools::BamAlignment & aln, bool processed) {
	bool succes = bReader.GetNextAlignment(aln);
	if (succes) {
		//read = seqInfo(aln.Name, aln.QueryBases, aln.Qualities, SangerQualOffset);
		read = bamAlnToSeqInfo(aln, false);
	}
	return succes;
}

bool SeqInput::readNextFastaStream(std::istream & fastaFile, seqInfo& read,
		bool processed) {
	if (!fastaFile.good()) {
		return false;
	}
	std::string name = "";
	std::string buildingSeq = "";
	std::string line = "";
	if ('>' == fastaFile.peek()) {
		njh::files::crossPlatGetline(fastaFile, name);
		while (fastaFile.peek() != std::ifstream::eofbit && fastaFile.good() && fastaFile.peek() != '>') {
			njh::files::crossPlatGetline(fastaFile, line);
			buildingSeq.append(line);
		}
		if (!ioOptions_.includeWhiteSpaceInName_ && name.find(" ") != std::string::npos) {
			//not really safe if name starts with space but hopefully no would do that
			read = seqInfo(name.substr(1, name.find(" ") - 1), buildingSeq);
		} else {
			read = seqInfo(name.substr(1), buildingSeq);
		}
		if (processed) {
			read.processRead(processed);
		}
		return true;
	} else {
		std::stringstream ss;
		ss << "error in reading fasta file, " << ioOptions_.firstName_ << " in " << __PRETTY_FUNCTION__ << ", line doesn't begin with >, starts with: "
			 << std::endl;
		ss << fastaFile.peek() << std::endl;
		throw std::runtime_error { ss.str() };
		return false;
	}
}

bool SeqInput::readNextQualStream(std::istream & qualFile,
		std::vector<uint8_t>& quals, std::string & name) {
	name.clear();
	quals.clear();

	if (!qualFile.good()) {
		return false;
	}
	std::string buildingQual = "";
	std::string line = "";
	if ('>' == qualFile.peek()) {
		njh::files::crossPlatGetline(qualFile, name);
		while (qualFile.peek() != std::ifstream::eofbit && qualFile.good()
				&& qualFile.peek() != '>') {
			njh::files::crossPlatGetline(qualFile, line);
			if ("" != buildingQual && ' ' != buildingQual.back()
					&& ' ' != line.front()) {
				buildingQual.push_back(' ');
			}
			buildingQual.append(line);
		}
		if (!ioOptions_.includeWhiteSpaceInName_
				&& name.find(" ") != std::string::npos) {
			//not really safe if name starts with space but hopefully no would do that
			name = name.substr(1, name.find(" ") - 1);
		} else {
			name = name.substr(1);
		}
		quals = stringToVector<uint8_t>(buildingQual);
		return true;
	} else {
		std::stringstream ss;
		ss << "error in reading fasta file in " << __PRETTY_FUNCTION__
				<< ", line doesn't begin with >, starts with: " << std::endl;
		ss << qualFile.peek() << std::endl;
		throw std::runtime_error { ss.str() };
		return false;
	}
}


bool SeqInput::readNextFastaQualStream(std::istream& fastaReader,
		std::istream& qualReader, seqInfo & read, bool processed) {
	if (!readNextFastaStream(fastaReader, read, processed)) {
		return false;
	}
	std::vector<uint8_t> quals;
	std::string name;
	bool qualStatus = readNextQualStream(qualReader, quals, name);
	if (qualStatus) {
		if (name == read.name_) {
			read.addQual(quals);
			return true;
		} else {
			std::stringstream ss;
			ss
					<< "Error in reading fasta and qual files, name in qual did not match fasta file, check ordering"
					<< std::endl;
			ss << "Fasta seq name: " << read.name_ << std::endl;
			ss << "Qual seq name: " << name << std::endl;
			throw std::runtime_error { ss.str() };
			return false;
		}
	} else {
		std::stringstream ss;
		ss << "Error in reading qual files in readNextFastaQualStream" << std::endl;
		throw std::runtime_error { ss.str() };
		return false;
	}
}


bool SeqInput::readNextFastqStream(const VecStr & data, const uint32_t lCount, uint8_t offSet, seqInfo& seq,
		bool processed){
	if (lCount == 4) {
		if (data[0].size() < 1 || data[0][0] != '@' || data[2].size() < 1
				|| data[2][0] != '+') {
			std::stringstream ss;
			ss << "Error in processing record:\n" << vectorToString(data, "\n")
					<< std::endl;
			if (data[0].size() == 0) {
				ss << "First line is empty" << std::endl;
			} else if (data[0][0] != '@') {
				ss << "First line did not begin with @" << std::endl;
			}
			if (data[2].size() == 0) {
				ss << "Third line is empty" << std::endl;
			} else if (data[2][0] != '+') {
				ss << "Third line did not begin with +" << std::endl;
			}
			throw std::runtime_error { ss.str() };
		}
		seq = seqInfo(data[0].substr(1), data[1], data[3], offSet);
		if (!ioOptions_.includeWhiteSpaceInName_ && seq.name_.find(" ") != std::string::npos) {
			seq.name_ = seq.name_.substr(0, seq.name_.find_first_of(" "));
		}
		seq.processRead(processed);
		return true;
	} else if (lCount > 0 && lCount < 4) {
		bool allBlanks = true;
		for (uint32_t i = 0; i < lCount; ++i) {
			if (data[i] != "") {
				allBlanks = false;
				break;
			}
		}
		if (allBlanks) {
			return false;
		}
		std::stringstream ss;
		ss << njh::bashCT::bold << "Incomplete sequence, read only"
				<< njh::bashCT::reset << std::endl;
		for (uint32_t i = 0; i < lCount; ++i) {
			ss << "!" << data[i] << "!" << std::endl;
		}
		ss << njh::bashCT::bold << "exiting" << njh::bashCT::reset << std::endl;
		throw std::runtime_error { ss.str() };
		return false;
	} else {
		return false;
	}
}

bool SeqInput::readNextFastqStream(std::istream& is, uint8_t offSet,
		seqInfo& seq, bool processed) {
	// assumes that there is no wrapping of lines, lines go name, seq, comments,
	// qual
	VecStr data(4, "");
	uint32_t count = 0;
	while (!is.eof()) {
		if (count > 3) {
			break;
		} else {
			njh::files::crossPlatGetline(is, data[count]);
		}
		++count;
	}
	return readNextFastqStream(data, count, offSet, seq, processed);
}



bool SeqInput::readNextSff(std::istream & inFile, sffObject & read) {
	if (!inFile.good()) {
		return false;
	}
	std::unordered_map<std::string, std::string> info;
	std::string line = "";
	if (inFile.peek() == '>') {
		std::getline(inFile, line);
		info["Name"] = line.substr(1);
	} else {
		std::cerr << "Error in reading in: readNextSff" << std::endl;
	}
	while (inFile.good() && '>' != inFile.peek()) {
		std::getline(inFile, line);
		if (line.empty() || allWhiteSpaceStr(line)) {
			continue;
		}
		auto toks = tokenizeString(line, ":");
		if (toks.size() != 2) {
			std::cerr << "!" << line << "!" << std::endl;
			throw std::runtime_error { "readNextSff: Error in parsing " + line
					+ " , should be two items sep by :" };
		}
		info[trimEndWhiteSpaceReturn(toks[0])] = trimEndWhiteSpaceReturn(toks[1]);
	}
	read.addSffInfo(info);
	return true;
}

bool SeqInput::readNextSff(std::istream & inFile, seqInfo & read) {
	lastSffRead_ = std::make_unique<sffObject>();
	bool succes = readNextSff(inFile, *lastSffRead_);
	if (succes) {
		read = lastSffRead_->seqBase_;
		return true;
	}
	return false;
}

VecStr SeqInput::readSffTxtHeader(std::istream & inFile) {
	VecStr header;
	std::string line = "";
	while ('>' != inFile.peek()) {
		std::getline(inFile, line);
		header.emplace_back(line);
	}
	return header;
}


bool SeqInput::readNextSffBin(std::istream& in, sffObject& read,
		int32_t numFlowReads) {
	try {
		//uint64_t startSpotInFile = in.tellg();
		in.tellg();
		if (!in.eof()) {
			/*****************************************/
			//read header
			//read header length
			in.read(reinterpret_cast<char*>(&read.headerLength_),
					sizeof(read.headerLength_));
			read.headerLength_ = swap_uint16(read.headerLength_);

			//read name length
			in.read(reinterpret_cast<char*>(&read.nameLength_),
					sizeof(read.nameLength_));
			read.nameLength_ = swap_uint16(read.nameLength_);

			//read num bases
			in.read(reinterpret_cast<char*>(&read.numBases_), sizeof(read.numBases_));
			read.numBases_ = swap_uint32(read.numBases_);

			//read clip qual left
			in.read(reinterpret_cast<char*>(&read.clipQualLeft_),
					sizeof(read.clipQualLeft_));
			read.clipQualLeft_ = swap_uint16(read.clipQualLeft_);
			read.clipQualLeft_ = 5;

			//read clip qual right
			in.read(reinterpret_cast<char*>(&read.clipQualRight_),
					sizeof(read.clipQualRight_));
			read.clipQualRight_ = swap_uint16(read.clipQualRight_);

			//read clipAdapterLeft
			in.read(reinterpret_cast<char*>(&read.clipAdapterLeft_),
					sizeof(read.clipAdapterLeft_));
			read.clipAdapterLeft_ = swap_uint16(read.clipAdapterLeft_);

			//read clipAdapterRight
			in.read(reinterpret_cast<char*>(&read.clipAdapterRight_),
					sizeof(read.clipAdapterRight_));
			read.clipAdapterRight_ = swap_uint16(read.clipAdapterRight_);

			//read name
			char* nameBuffer = new char[read.nameLength_];
			in.read(&(*nameBuffer), read.nameLength_);
			read.seqBase_.name_ = nameBuffer;
			if (read.seqBase_.name_.length() > read.nameLength_) {
				read.seqBase_.name_ = read.seqBase_.name_.substr(0, read.nameLength_);
			}
			delete[] nameBuffer;

			//extract info from name
			read.decodeName();

			// Pad to 8 chars
			uint64_t spotInFile = in.tellg();
			uint64_t spot = (spotInFile + 7) & ~7;
			in.seekg(spot);

			//read seq
			read.flowValues.resize(numFlowReads);
			for (int32_t i = 0; i < numFlowReads; i++) {
				uint16_t tempFlow = 0;
				in.read(reinterpret_cast<char*>(&tempFlow), sizeof(tempFlow));
				read.flowValues[i] = swap_uint16(tempFlow) / 100.0;
			}
			//read baseFlowIndex
			read.baseFlowIndex_.resize(read.numBases_);
			for (uint32_t i = 0; i < read.numBases_; i++) {
				in.read(reinterpret_cast<char*>(&read.baseFlowIndex_[i]), 1);
			}
			//set read flowIndexes
			read.flowIndex_.resize(read.baseFlowIndex_.size());
			uint32_t sum = 0;
			for(uint32_t i = 0; i < read.baseFlowIndex_.size(); ++i){
				sum +=  read.baseFlowIndex_[i];
				read.flowIndex_[i] = sum;
			}
			//read bases
			char* seqBuffer = new char[read.numBases_];
			in.read(&(*seqBuffer), read.numBases_);
			read.seqBase_.seq_ = seqBuffer;
			if (read.seqBase_.seq_.length() > read.numBases_) {
				read.seqBase_.seq_ = read.seqBase_.seq_.substr(0, read.numBases_);
			}
			delete[] seqBuffer;

			//read qual scores
			read.seqBase_.qual_.resize(read.numBases_);
			for (uint32_t i = 0; i < read.numBases_; i++) {
				in.read(reinterpret_cast<char*>(&read.seqBase_.qual_[i]), 1);
			}

			// Pad to 8 chars
			spotInFile = in.tellg();
			spot = (spotInFile + 7) & ~7;
			in.seekg(spot);
		} else {
			std::stringstream ss;
			ss << "Error in reading sff binary object" << std::endl;
			throw std::runtime_error { ss.str() };
		}

		if (in.eof()) {
			return true;
		}
		return false;
	} catch (std::exception& e) {
		std::stringstream ss;
		ss << "Error in reading during readSeqData" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}

void SeqInput::readHeader(std::istream& in, sffBinaryHeader& header) {
	try {
		if (!in.eof()) {
			//read magic number
			in.read(reinterpret_cast<char*>(&header.magicNumber),
					sizeof(header.magicNumber));
			header.magicNumber = swap_uint32(header.magicNumber);

			//read version
			for (int i = 0; i < 4; i++) {
				int32_t tempNum = 0;
				in.read(reinterpret_cast<char*>(&tempNum), 1);
				header.version += estd::to_string(tempNum);
			}

			//read offset
			in.read(reinterpret_cast<char*>(&header.indexOffset),
					sizeof(header.indexOffset));
			header.indexOffset = swap_uint64(header.indexOffset);

			//read index length
			in.read(reinterpret_cast<char*>(&header.indexLength),
					sizeof(header.indexLength));
			header.indexLength = swap_uint32(header.indexLength);

			//read num reads
			in.read(reinterpret_cast<char*>(&header.numReads),
					sizeof(header.numReads));
			header.numReads = swap_uint32(header.numReads);

			//read header length
			in.read(reinterpret_cast<char*>(&header.headerLength),
					sizeof(header.headerLength));
			header.headerLength = swap_uint16(header.headerLength);

			//read key length
			in.read(reinterpret_cast<char*>(&header.keyLength),
					sizeof(header.keyLength));
			header.keyLength = swap_uint16(header.keyLength);

			//read number of flow reads
			in.read(reinterpret_cast<char*>(&header.numFlowsPerRead),
					sizeof(header.numFlowsPerRead));
			header.numFlowsPerRead = swap_uint16(header.numFlowsPerRead);

			//read format code
			in.read(reinterpret_cast<char*>(&header.flogramFormatCode), 1);

			//read flow chars
			char* tempBuffer = new char[header.numFlowsPerRead];
			in.read(&(*tempBuffer), header.numFlowsPerRead);
			header.flowChars = tempBuffer;
			if (header.flowChars.length() > header.numFlowsPerRead) {
				header.flowChars = header.flowChars.substr(0, header.numFlowsPerRead);
			}
			delete[] tempBuffer;

			//read key
			char* tempBuffer2 = new char[header.keyLength];
			in.read(&(*tempBuffer2), header.keyLength);
			header.keySequence = tempBuffer2;
			if (header.keySequence.length() > header.keyLength) {
				header.keySequence = header.keySequence.substr(0, header.keyLength);
			}
			delete[] tempBuffer2;

			/* Pad to 8 chars */
			uint64_t spotInFile = in.tellg();
			uint64_t spot = (spotInFile + 7) & ~7;  // ~ inverts
			in.seekg(spot);
		} else {
			std::stringstream ss;
			ss << "Error reading sff common header" << std::endl;
			throw std::runtime_error { ss.str() };
		}
	} catch (std::exception& e) {
		std::stringstream ss;
		ss << "Error in reading during readHeader" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}

}  // namespace njhseq

