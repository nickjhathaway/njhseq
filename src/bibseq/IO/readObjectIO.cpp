#include "readObjectIO.hpp"
#include <bibcpp/bashUtils.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
namespace bibseq {


void readObjectIO::read(const readObjectIOOptions& options) {
  includeSpaceInNames = options.includeWhiteSpaceInName_;
  read(options.inFormat_, options.firstName_, options.secondName_,
       options.thirdName_, options.processed_);
  readVec::handelLowerCaseBases(reads, options.lowerCaseBases_);
  if (options.removeGaps_) {
    std::cout << "mark removing gaps" << std::endl;
    readVec::removeGapsFromReads(reads);
  }
}

// main reader
void readObjectIO::read(const std::string& fileType, std::string filename,
                        std::string secondName, std::string thirdName,
                        bool processed) {
	if(!fexists(filename)){
		std::stringstream ss;
		ss << "error, file " << filename << " doesn't exist " << std::endl;
		throw std::runtime_error{ss.str()};
	}
	std::ifstream inFile(filename);
  if (fileType == "clustal") {
    readClustal(filename, processed);
  } else if (fileType == "clustal-ng") {
    readClustalNg(filename, processed);
  } else if (fileType == "sff") {
    readSff(filename);
  } else if (fileType == "sffbin") {
    readSffbin(filename);
  }else if (fileType == "shorah") {
    readShorah(filename);
  } else if (fileType == "shorahOld") {
    readShorahOld(filename);
  } else if (fileType == "fastq" || fileType == "fq" || fileType == "fnq") {
  	readFastqStream(inFile, SangerQualOffset ,processed, false);
  } else if (fileType == "fastqgz" ) {
  	readFastqStreamGz(filename, SangerQualOffset ,processed, false);
  } else if (fileType == "fasta") {
  	readFastaStream(inFile, processed, false);
  } else if (fileType == "bam") {
    readBam(filename, processed);
  } else if (fileType == "fastaQual") {
    readFastaQual(filename, secondName, processed);
  } else {
  	std::stringstream ss;
    ss << "Unrecognized file type : " << fileType << ", not reading "
              << filename << std::endl;
    ss << "Acceptable types are fasta, qual,fastq, stub, sff,  "
              << std::endl;
    throw std::runtime_error{ss.str()};
  }
}


// single readers
void readObjectIO::read(const std::string& fileType, std::string filename,
                        bool processed) {
  read(fileType, filename, "noSecondName", "noThirdName", processed);
}

// double readers
void readObjectIO::read(const std::string& fileType, std::string filename,
                        std::string secondName, bool processed) {
  read(fileType, filename, secondName, "noThirdName", processed);
}

void readObjectIO::readShorahOld(const std::string& shorahFilename) {
	if(!fexists(shorahFilename)){
		std::stringstream ss;
		ss << "file: " << shorahFilename << ", doesn't exist or can't be open" << std::endl;
		throw std::runtime_error{ss.str()};
	}else{
		std::ifstream inFile(shorahFilename);
	  readFastaStream(inFile, false, false);
	  // allRemoveGaps(reads);
	  for (auto & read :reads) {
	    VecStr toks = tokenizeString(read.seqBase_.name_, "_");
	    read.seqBase_.frac_ = std::stod(toks[toks.size() - 1]);
	    read.seqBase_.cnt_ = read.seqBase_.frac_ * 100000;
	  }
	  readVec::allUpdateName(reads);
	}

}

void readObjectIO::readShorah(const std::string& shorahFilename) {
	if(!fexists(shorahFilename)){
		std::stringstream ss;
		ss << "file: " << shorahFilename << ", doesn't exist or can't be open" << std::endl;
		throw std::runtime_error{ss.str()};
	}else{
		std::ifstream inFile(shorahFilename);
	  readFastaStream(inFile, false, false);
	  for (auto & read : reads) {
	    VecStr toks = tokenizeString(read.seqBase_.name_, "|");
	    read.seqBase_.name_ = toks[0];
	    VecStr secondToks = tokenizeString(toks[1], "=");
	    read.seqBase_.cnt_ = std::stod(secondToks.back());
	  }
	  readVec::allUpdateName(reads);
	}
}

void readObjectIO::readBam(std::string filename, bool processed) {
  reads.clear();
  BamTools::BamReader bReader;
  bReader.Open(filename);
  BamTools::BamAlignment aln;
  readObject read;

  while (readNextBam(bReader, read, aln, processed)) {
    reads.emplace_back(read);
  }
}

bool readObjectIO::readNextBam(BamTools::BamReader & bReader, readObject& read,
		BamTools::BamAlignment & aln, bool processed){
	bool succes = bReader.GetNextAlignment(aln);
	if(succes){
		read = readObject(seqInfo(aln.Name, aln.QueryBases,
				aln.Qualities,
        readObjectIO::SangerQualOffset));
	}
	return succes;
}




void readObjectIO::readFastaQual(std::string fastaName,
                                 std::string qualFilename,
																 bool processed) {
	//std::cout << "here1" << std::endl;


  reads.clear();
  std::stringstream ss;
  bool fail = false;
	if(!fexists(fastaName)){
		ss << "file: " << fastaName << ", doesn't exist or can't be open" << std::endl;
		fail = true;
	}
	if (!fexists(qualFilename)){
		ss << "file: " << qualFilename << ", doesn't exist or can't be open" << std::endl;
		fail = true;
	}
	if(fail){
		throw std::runtime_error{ss.str()};
	} else {
		std::ifstream inFastaFile(fastaName);
		std::ifstream inQualFile(qualFilename);
		if(!inFastaFile){
			ss << "Error in opening in fastaFile: " << fastaName << std::endl;
			fail = true;
		}
		if(!inQualFile){
			ss << "Error in opening in qualFile: " << qualFilename << std::endl;
			fail = true;
		}
		if(fail){
			throw std::runtime_error{ss.str()};
		}else{
			//std::cout << "here2" << std::endl;
			readFastaQualStream(inFastaFile, inQualFile, processed, false);
		}
	}
}


bool readObjectIO::readNextFastaStream(cachedReader & cReader,
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
/*
bool readObjectIO::readNextFastaStream(std::istream& is,
		readObject& read, bool processed) {
  cachedReader cReader(is);
  return readNextFastaStream(cReader, read, processed);
}

bool readObjectIO::readNextQualStream(std::istream& is, std::vector<uint32_t>& quals, std::string & name){
  cachedReader cReader(is);
  return readNextQualStream(cReader, quals, name );
}*/

void readObjectIO::readFastaQualStream(std::istream& fastaReader,std::istream& qualReader, bool processed, bool add){
	if(!add){
		reads.clear();
	}
  readObject tempObj;
  cachedReader fastaCacheReader(fastaReader);
  cachedReader qualCacheReader(qualReader);
  //std::cout << "here3" << std::endl;
  while(readNextFastaQualStream(fastaCacheReader, qualCacheReader, tempObj, processed)){
  	reads.emplace_back(tempObj);
  }
}
bool readObjectIO::readNextFastaQualStream(cachedReader& fastaReader, cachedReader& qualReader, readObject & read,bool processed){
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
/*
bool readObjectIO::readNextFastaQualStream(std::istream& fastaReader, std::istream& qualReader, readObject & read,bool processed){
  cachedReader fastaCacheReader(fastaReader);
  cachedReader qualCacheReader(qualReader);
	return readNextFastaQualStream(fastaCacheReader, qualCacheReader, read, processed);
}*/

bool readObjectIO::readNextQualStream(cachedReader& cReader, std::vector<uint32_t>& quals, std::string & name){
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

void readObjectIO::readFastaStream(std::istream& is, bool processed, bool add) {
	if(!add){
		reads.clear();
	}
  readObject tempObj;
  cachedReader cReader(is);
  while(readNextFastaStream(cReader, tempObj, processed)){
  	reads.emplace_back(tempObj);
  }
}


bool readObjectIO::readNextFastqStream(std::istream& is, uint32_t offSet,
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

void readObjectIO::readFastqStream(std::istream& is, uint32_t offSet,
                                   bool processed, bool add) {
	if(!add){
	  reads.clear();
	}
  readObject read;
  while (readNextFastqStream(is, offSet, read, processed)) {
    reads.emplace_back(read);
  }
}

void readObjectIO::readFastqStreamGz(std::string filename, uint32_t offSet,
                                   bool processed, bool add){
  if(!fexists(filename)){
  	throw std::runtime_error{"file " + filename + " doesn't exist"};
  }
	std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);

  boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
  inbuf.push(boost::iostreams::gzip_decompressor());
  inbuf.push(file);
  //Convert streambuf to istream
  std::istream instream(&inbuf);
  readFastqStream(instream, offSet, processed, add);
  file.close();
}





void readObjectIO::readClustal(std::string filename, bool processed) {
	reads.clear();
	table inTab(filename);
	std::vector<std::pair<std::string, readObject>> readMap;
	std::vector<std::pair<std::string, readObject>>::iterator readMapIter;
	for (const auto& row : inTab.content_) {
		if (row.size() != 2 || row[0][0] == '*') {
		} else {
			bool foundMatch = false;
			for (readMapIter = readMap.begin(); readMapIter != readMap.end();
					++readMapIter) {
				if (readMapIter->first == row[0]) {
					foundMatch = true;
					readMapIter->second.seqBase_.seq_.append(row[1]);
					break;
				}
			}
			if (!foundMatch) {
				readMap.emplace_back(row[0], readObject(seqInfo(row[0], row[1])));
			}
		}
	}
	for (readMapIter = readMap.begin(); readMapIter != readMap.end();
			readMapIter++) {
		reads.emplace_back(readMapIter->second);
	}
}
void readObjectIO::readClustalNg(std::string filename, bool processed) {
  reads.clear();
  readClustal(filename, processed);
  for (auto& read : reads) {
    seqUtil::removeGaps(read.seqBase_.seq_);
  }
}

void readObjectIO::readSff(std::string filename) {
  reads.clear();
  sffReads.clear();
  std::fstream sffFile;
  sffFile.open(filename.c_str());
  if (!sffFile) {
  	std::stringstream ss;
    ss << "Error in opening " << filename;
    throw std::runtime_error{ss.str()};
  }
  std::string currentLine;
  getline(sffFile, currentLine);
  int counter = 1;
  bool endWhile = false;
  while (!endWhile) {
    if (sffFile.eof()) {
      endWhile = true;
      break;
    }
    if (counter % 5000 == 0) {
      std::cout << counter << std::endl;
    }
    sffObject nextRead = readNextSff(currentLine, sffFile);
    if (nextRead.seqBase_.seq_.length() != 0) {
      sffReads.emplace_back(nextRead);
    }
    ++counter;
  }
  int convertCount = 0;
  // std::cout<<"Converting sff reads into regular read objects"<<std::endl;
  for (const auto& read : sffReads) {
    reads.emplace_back(convertSffObject(read));
    ++convertCount;
  }
}


///////////////////////reading the files


sffObject readObjectIO::readNextSff(std::string& currentLine,
                                    std::fstream& sffFile) {
  sffObject tempReadObject = sffObject();
  std::unordered_map<std::string, std::string> info;
  if (currentLine[0] == '>') {
    if (currentLine[0] == '>') {
      currentLine.erase(currentLine.begin());
      info["Name"] = currentLine;
      //info.insert(std::make_pair("Name", currentLine));
    }
    std::string tempTitleString;
    std::string tempInfoString;
    getline(sffFile, currentLine, ':');
    // gather the info about the read
    bool endInfoGather = false;
    while (!endInfoGather) {
      if (currentLine[0] == 'Q') {
        endInfoGather = true;
      }
      tempTitleString = currentLine;
      getline(sffFile, currentLine);
      tempInfoString = currentLine;
      trimEndWhiteSpace(tempTitleString);
      trimEndWhiteSpace(tempInfoString);

      //info.insert(make_pair(tempTitleString, tempInfoString));
      info[tempTitleString] = tempInfoString;
      tempTitleString.clear();
      tempInfoString.clear();
      if (endInfoGather) {
        getline(sffFile, currentLine);
        getline(sffFile, currentLine);
      } else {
        getline(sffFile, currentLine, ':');
      }
    }
    tempReadObject.addSffInfo(info);
  } else {
    getline(sffFile, currentLine);
  }
  return tempReadObject;
}

const uint32_t readObjectIO::IlluminaQualOffset = 64;
const uint32_t readObjectIO::SangerQualOffset = 33;

std::vector<readObject> readObjectIO::getReferenceSeq(
    const std::string& refFilename, const std::string& refFormat,
    bool refProcessed, uint64_t& maxLength) {
  readObjectIO refReader;
  refReader.read(refFormat, refFilename, refProcessed);
  readVec::getMaxLength(refReader.reads, maxLength);
  return refReader.reads;
}


bool readObjectIO::readNextSffBin(std::ifstream& in, sffObject& read, int32_t numFlowReads){
	try {
		//uint64_t startSpotInFile = in.tellg();
		in.tellg();
		if (!in.eof()) {

			/*****************************************/
			//read header

			//read header length
			in.read(reinterpret_cast<char*>(&read.headerLength_), sizeof(read.headerLength_));
			read.headerLength_ = swap_uint16(read.headerLength_);

			//read name length
			in.read(reinterpret_cast<char*>(&read.nameLength_), sizeof(read.nameLength_));
			read.nameLength_ = swap_uint16(read.nameLength_);

			//read num bases
			in.read(reinterpret_cast<char*>(&read.numBases_), sizeof(read.numBases_));
			read.numBases_ = swap_uint32(read.numBases_);

			//read clip qual left
			in.read(reinterpret_cast<char*>(&read.clipQualLeft_), sizeof(read.clipQualLeft_));
			read.clipQualLeft_ = swap_uint16(read.clipQualLeft_);
			read.clipQualLeft_ = 5;


			//read clip qual right
			in.read(reinterpret_cast<char*>(&read.clipQualRight_), sizeof(read.clipQualRight_));
			read.clipQualRight_ = swap_uint16(read.clipQualRight_);

			//read clipAdapterLeft
			in.read(reinterpret_cast<char*>(&read.clipAdapterLeft_), sizeof(read.clipAdapterLeft_));
			read.clipAdapterLeft_ = swap_uint16(read.clipAdapterLeft_);

			//read clipAdapterRight
			in.read(reinterpret_cast<char*>(&read.clipAdapterRight_), sizeof(read.clipAdapterRight_));
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
			uint64_t spot = (spotInFile + 7)& ~7;
			in.seekg(spot);

			//read seq
			read.flowValues.resize(numFlowReads);
			for (int32_t i = 0; i < numFlowReads; i++) {
				uint16_t tempFlow = 0;
				in.read(reinterpret_cast<char*>(&tempFlow), sizeof(tempFlow));
				read.flowValues[i] = swap_uint16(tempFlow) / 100.0;
			}
			//read flowIndex
			read.flowIndex_.resize(read.numBases_);
			for (uint32_t i = 0; i < read.numBases_; i++) {
				in.read(reinterpret_cast<char*>(&read.flowIndex_[i]), 1);
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
			spot = (spotInFile + 7)& ~7;
			in.seekg(spot);
		}else{
			std::stringstream ss;
			ss << "Error in reading sff binary object" << std::endl;
			throw std::runtime_error{ss.str()};
		}

    if (in.eof()) {
    	return true;
    }

		return false;
	}
	catch(std::exception& e) {
		std::stringstream ss;
		ss << "Error in reading during readSeqData" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error{ss.str()};
	}
}

void readObjectIO::readHeader(std::ifstream& in, sffBinaryHeader& header){
	try {

		if (!in.eof()) {
			//read magic number
			in.read(reinterpret_cast<char*>(&header.magicNumber), sizeof(header.magicNumber));
			header.magicNumber = swap_uint32(header.magicNumber);

			//read version
			for (int i = 0; i < 4; i++) {
				int32_t tempNum = 0;
				in.read(reinterpret_cast<char*>(&tempNum), 1);
				header.version += to_string(tempNum);
			}

			//read offset
			in.read(reinterpret_cast<char*>(&header.indexOffset), sizeof(header.indexOffset));
			header.indexOffset = swap_uint64(header.indexOffset);

			//read index length
			in.read(reinterpret_cast<char*>(&header.indexLength), sizeof(header.indexLength));
			header.indexLength = swap_uint32(header.indexLength);

			//read num reads
			in.read(reinterpret_cast<char*>(&header.numReads), sizeof(header.numReads));
			header.numReads = swap_uint32(header.numReads);

			//read header length
			in.read(reinterpret_cast<char*>(&header.headerLength), sizeof(header.headerLength));
			header.headerLength = swap_uint16(header.headerLength);

			//read key length
			in.read(reinterpret_cast<char*>(&header.keyLength), sizeof(header.keyLength));
			header.keyLength = swap_uint16(header.keyLength);

			//read number of flow reads
			in.read(reinterpret_cast<char*>(&header.numFlowsPerRead), sizeof(header.numFlowsPerRead));
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
			if (header.keySequence.length() > header.keyLength) { header.keySequence = header.keySequence.substr(0, header.keyLength);  }
			delete[] tempBuffer2;

			/* Pad to 8 chars */
			uint64_t spotInFile = in.tellg();
			uint64_t spot = (spotInFile + 7)& ~7;  // ~ inverts
			in.seekg(spot);
			}else{
				std::stringstream ss;
				ss << "Error reading sff common header" << std::endl;
				throw std::runtime_error{ss.str()};
		}
	}
	catch(std::exception& e) {
		std::stringstream ss;
		ss << "Error in reading during readHeader" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error{ss.str()};
	}
}

void readObjectIO::readSffbin(const std::string & filename){
	std::ifstream inFile;
	if(!fexists(filename)){
		std::stringstream ss;
		ss << "File: " << filename << " doesn't exist" << std::endl;
		throw std::runtime_error{ss.str()};
	}
	inFile.open(filename, std::ios::binary);
	sffBinaryHeader header;
	readHeader(inFile, header);
	uint32_t count = 0;
	while (!inFile.eof()) {
		//read data
		sffObject read;
		readNextSffBin(inFile, read, header.numFlowsPerRead);
		bool okay = read.sanityCheck();
		if (!okay) {
			break;
		}
		/**@todo find out why numberOfFlows isn't being set properly*/
		reads.emplace_back(convertSffObject(read));

		++count;
		if((count + 1) % 10000 == 0){
			std::cout << count + 1 << std::endl;
		}
		if (count >= header.numReads) {
			break;
		}
	}
}


}  // namespace bib
