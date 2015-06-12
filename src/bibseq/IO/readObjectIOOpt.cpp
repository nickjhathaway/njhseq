#include "readObjectIOOpt.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {


readObjectIOOpt::readObjectIOOpt(const readObjectIOOptions & options) :
		ioOptions_(options) {

}

bool readObjectIOOpt::inOpen()const{
	return inOpen_;
}
void readObjectIOOpt::openIn(){
	if(!bib::files::bfs::exists(ioOptions_.firstName_)){
  	std::stringstream ss;
    ss << "Error file: " << ioOptions_.firstName_ << "doesn't exist\n";
    throw std::runtime_error{ss.str()};
	}
	bool failedToOpen = false;
  if (ioOptions_.inFormat_ == "fastq" || ioOptions_.inFormat_ == "fq" || ioOptions_.inFormat_ == "fnq") {
  	fqReader_ = std::make_unique<std::ifstream>(ioOptions_.firstName_);
  	if(!(*fqReader_)){
  		failedToOpen = true;
  	}
  } else if (ioOptions_.inFormat_ == "fastqgz" ) {
  	fqReaderBackup_ = std::make_unique<std::ifstream>(ioOptions_.firstName_, std::ios_base::in | std::ios_base::binary);
    fqInputBuffer_ = std::make_unique<boost::iostreams::filtering_streambuf<boost::iostreams::input>>();
    (*fqInputBuffer_).push(boost::iostreams::gzip_decompressor());
    (*fqInputBuffer_).push(*fqReaderBackup_);
    //Convert streambuf to istream
    fqReader_ = std::make_unique<std::istream>(&(*fqInputBuffer_));
  	if(!(*fqReader_)){
  		failedToOpen = true;
  	}
  } else if (ioOptions_.inFormat_ == "fasta") {
  	fStreamReader_ = std::make_unique<std::ifstream>(ioOptions_.firstName_);
  	if(!(*fStreamReader_)){
  		failedToOpen = true;
  	}
  	fReader_ = std::make_unique<cachedReader>(*fStreamReader_);
  } else if (ioOptions_.inFormat_ == "bam") {
  	bReader_.Open(ioOptions_.firstName_);
  	if(!bReader_.IsOpen()){
  		failedToOpen = true;
  	}
  } else if (ioOptions_.inFormat_ == "fastaQual") {
  	fStreamReader_ = std::make_unique<std::ifstream>(ioOptions_.firstName_);
  	if(!(*fStreamReader_)){
  		failedToOpen = true;
  	}
  	fReader_ = std::make_unique<cachedReader>(*fStreamReader_);

  	qStreamReader_ = std::make_unique<std::ifstream>(ioOptions_.secondName_);
  	if(!(*qStreamReader_)){
  		failedToOpen = true;
  	}
  	qReader_ = std::make_unique<cachedReader>(*qStreamReader_);
  } else if (ioOptions_.inFormat_ == "sff") {
  	fStreamReader_ = std::make_unique<std::ifstream>(ioOptions_.firstName_);
  	if(!(*fStreamReader_)){
  		failedToOpen = true;
  	}
  	sffTxtHeader_ = readSffTxtHeader(*fStreamReader_);
  }  else {
  	std::stringstream ss;
    ss << "Unrecognized file type : " << ioOptions_.inFormat_ << ", not reading "
              << ioOptions_.firstName_ << std::endl;
    ss << "Acceptable types are fasta, qual,fastq, stub, sff,  "
              << std::endl;
    throw std::runtime_error{ss.str()};
  }

  if(failedToOpen){
  	std::stringstream ss;
    ss << "Error in opening : " << ioOptions_.firstName_ << "\n";

    throw std::runtime_error{ss.str()};
  }
  inOpen_ = true;
}

void readObjectIOOpt::closeIn(){
	fReader_ = nullptr;
	qReader_ = nullptr;
	fStreamReader_ = nullptr;
	qStreamReader_ = nullptr;

  fqInputBuffer_ = nullptr;
  fqReaderBackup_ = nullptr;
  fqReader_ = nullptr;
  if(bReader_.IsOpen()){
  	bReader_.Close();
  }
  inOpen_ = false;
}

bool readObjectIOOpt::readNextRead(readObject & read){
  if (ioOptions_.inFormat_ == "fastq" || ioOptions_.inFormat_ == "fq" || ioOptions_.inFormat_ == "fnq") {
  	return readNextFastqStream(*fqReader_, SangerQualOffset, read, ioOptions_.processed_);
  } else if (ioOptions_.inFormat_ == "fastqgz" ) {
  	return readNextFastqStream(*fqReader_, SangerQualOffset, read, ioOptions_.processed_);
  } else if (ioOptions_.inFormat_ == "fasta") {
  	return readNextFastaStream(*fReader_, read, ioOptions_.processed_);
  } else if (ioOptions_.inFormat_ == "bam") {
  	return readNextBam(bReader_, read, aln_, ioOptions_.processed_);
  } else if (ioOptions_.inFormat_ == "fastaQual") {
  	return readNextFastaQualStream(*fReader_,*qReader_,read, ioOptions_.processed_);
  } else if (ioOptions_.inFormat_ == "sff") {
  	return readNextSff(*fStreamReader_,read);
  } else {
  	std::cerr << "shouldn't be happening" << ", readnextRead in readIoOpt" << "\n";
  }
  return false;
}

bool readObjectIOOpt::readNextBam(BamTools::BamReader & bReader, readObject& read,
		BamTools::BamAlignment & aln, bool processed){
	bool succes = bReader.GetNextAlignment(aln);
	if(succes){
		read = readObject(seqInfo(aln.Name, aln.QueryBases,
				aln.Qualities,
        readObjectIOOpt::SangerQualOffset));
	}
	return succes;
}

bool readObjectIOOpt::readNextFastaStream(cachedReader & cReader,
		readObject& read, bool processed) {
	if(cReader.done()){
		return false;
	}
	//std::cout << "reading next" << std::endl;
  std::string name = "";
  std::string buildingSeq = "";
  if(cReader.currentLine().front() =='>'){
  	name = cReader.currentLine();
  	//std::cout << "name: " << name << std::endl;
  	while( cReader.setNextLine() && cReader.currentLine().front() !='>'){
  		//std::cout <<"seq:" << cReader.currentLine() << std::endl;
			buildingSeq.append(cReader.currentLine());
		}
  	read = readObject(seqInfo(name.substr(1), buildingSeq));
  	if(processed){
  		read.processRead(processed);
  	}
  	return true;
  }else{
  	std::stringstream ss;
  	ss << "error in reading fasta file, line doesn't begin with >, starts with: " << std::endl;
  	ss << cReader.currentLine().front() << std::endl;
  	throw std::runtime_error{ss.str()};
  	return false;
  }
}


bool readObjectIOOpt::readNextFastaQualStream(cachedReader& fastaReader, cachedReader& qualReader, readObject & read,bool processed){
	if(!readNextFastaStream(fastaReader, read, processed)){
		return false;
	}
	std::vector<uint32_t> quals;
	std::string name;
	bool qualStatus = readNextQualStream(qualReader, quals, name);
	if(qualStatus){
		if(name == read.seqBase_.name_){
			read.seqBase_.addQual(quals);
			return true;
		}else{
			std::stringstream ss;
			ss << "Error in reading fasta and qual files, name in qual did not match fasta file, check ordering" << std::endl;
			ss << "Fasta seq name: " << read.seqBase_.name_ << std::endl;
			ss << "Qual seq name: " << name << std::endl;
			throw std::runtime_error{ss.str()};
			return false;
		}
	}else{
		std::stringstream ss;
		ss << "Error in reading qual files in readNextFastaQualStream" << std::endl;
		throw std::runtime_error{ss.str()};
		return false;
	}
}


bool readObjectIOOpt::readNextQualStream(cachedReader& cReader, std::vector<uint32_t>& quals, std::string & name){
	if(cReader.done()){
		return false;
	}
	quals.clear();
  if(cReader.currentLine().front() =='>'){
  	name = cReader.currentLine().substr(1);
  	while( cReader.setNextLine() && cReader.currentLine().front() !='>'){
  		addOtherVec(quals, stringToVector<uint32_t>(cReader.currentLine()));
		}
  	return true;
  }else{
  	std::stringstream ss;
  	ss << "error in reading qual file, line doesn't begin with >, starts with: " << std::endl;
  	ss << cReader.currentLine().front() << std::endl;

  	throw std::runtime_error{ss.str()};
  	return false;
  }
}



bool readObjectIOOpt::readNextFastqStream(std::istream& is, uint32_t offSet,
                                       readObject& read, bool processed) {
  // assumes that there is no wrapping of lines, lines go name, seq, comments,
  // qual
  VecStr data(4);
  uint32_t count = 0;
  while (!is.eof()) {
    if (count > 3) {
      break;
    } else {
      std::getline(is, data[count]);
    }
    ++count;
  }
  if (count == 4) {
    data[0].erase(0, 1);
    read.seqBase_ = seqInfo(data[0], data[1], data[3], offSet);
    read.processRead(processed);
    return true;
  } else if (count > 0 && count < 4) {
    bool allBlanks = true;
    for (uint32_t i = 0; i < count; ++i) {
      if (data[i] != "") {
        allBlanks = false;
        break;
      }
    }
    if (allBlanks) {
      return false;
    }
    std::stringstream ss;
    ss << bib::bashCT::bold
    		<< "Incomplete sequence, read only"
    		<< bib::bashCT::reset << std::endl;
    for (uint32_t i = 0; i < count; ++i) {
      ss << "!" << data[i] << "!" << std::endl;
    }
    ss << bib::bashCT::bold
    		<< "exiting" << bib::bashCT::reset << std::endl;
    throw std::runtime_error{ss.str()};
    return false;
  } else {
    return false;
  }
}

bool readObjectIOOpt::outOpen()const {
	return outOpen_;
};

void readObjectIOOpt::openOut() {
	if (ioOptions_.outFormat_ == "fastq") {
		fOut_ = std::make_unique<std::ofstream>();
		bib::files::openTextFile(*fOut_, ioOptions_.outFilename_,
				ioOptions_.outExtention_, ioOptions_.overWriteFile_, ioOptions_.append_,
				ioOptions_.exitOnFailureToWrite_);
	} else if (ioOptions_.outFormat_ == "fastq") {
		fOut_ = std::make_unique<std::ofstream>();
		bib::files::openTextFile(*fOut_, ioOptions_.outFilename_,
				ioOptions_.outExtention_, ioOptions_.overWriteFile_, ioOptions_.append_,
				ioOptions_.exitOnFailureToWrite_);
	} else if (ioOptions_.outFormat_ == "fastaQual") {
		fOut_ = std::make_unique<std::ofstream>();
		bib::files::openTextFile(*fOut_, ioOptions_.outFilename_,
				".fasta", ioOptions_.overWriteFile_, ioOptions_.append_,
				ioOptions_.exitOnFailureToWrite_);
		qOut_ = std::make_unique<std::ofstream>();
		bib::files::openTextFile(*qOut_, ioOptions_.outFilename_,
				".fasta.qual", ioOptions_.overWriteFile_, ioOptions_.append_,
				ioOptions_.exitOnFailureToWrite_);
	} else if (ioOptions_.outFormat_ == "mothurData") {
		fOut_ = std::make_unique<std::ofstream>();
		bib::files::openTextFile(*fOut_, ioOptions_.outFilename_,
				".dat", ioOptions_.overWriteFile_, ioOptions_.append_,
				ioOptions_.exitOnFailureToWrite_);
		(*fOut_) << ioOptions_.extra_ << std::endl;
	}  else {
		std::stringstream ss;
		ss << "Unrecognized out file type : " << ioOptions_.outFormat_ << std::endl;
		ss << "Acceptable types are fasta, fastaQual" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	outOpen_ = true;
}

void readObjectIOOpt::closeOut() {
	fOut_ = nullptr;
	qOut_ = nullptr;
	outOpen_ = false;
}



bool readObjectIOOpt::readNextSff(std::ifstream & inFile, sffObject & read){
	if(!inFile.good()){
		return false;
	}
	std::unordered_map<std::string, std::string> info;
	std::string line = "";
	if(inFile.peek() == '>'){
		std::getline(inFile, line);
		info["Name"] = line.substr(1);
	}else{
		std::cerr << "Error in reading in: readNextSff" << std::endl;
	}
	while(inFile.good() && '>' != inFile.peek()){
		std::getline(inFile, line);
		if(line.empty() || allWhiteSpaceStr(line)){
			continue;
		}
		auto toks = tokenizeString(line, ":");
		if(toks.size() != 2){
			std::cerr << "!" << line << "!" << std::endl;
			throw std::runtime_error{"readNextSff: Error in parsing " + line + " , should be two items sep by :"};
		}
		info[trimEndWhiteSpaceReturn(toks[0])] = trimEndWhiteSpaceReturn(toks[1]);
	}
	read.addSffInfo(info);
	return true;
}

bool readObjectIOOpt::readNextSff(std::ifstream & inFile, readObject & read){
	sffObject sffRead;
	bool succes = readNextSff(inFile, sffRead);
	if(succes){
		read = convertSffObject(sffRead);
		return true;
	}
	return false;
}

VecStr readObjectIOOpt::readSffTxtHeader(std::ifstream & inFile){
	VecStr header;
	std::string line = "";
	while('>' != inFile.peek()){
		std::getline(inFile, line);
		header.emplace_back(line);
	}
	return header;
}


const uint32_t readObjectIOOpt::IlluminaQualOffset = 64;
const uint32_t readObjectIOOpt::SangerQualOffset = 33;



}  // namespace bib
