//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "readObjectIO.hpp"
#include <bibcpp/bashUtils.h>

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
		std::cout << "error, file " << filename << " doesn't exist " << std::endl;
		exit(1);
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
  } else if (fileType == "fasta") {
  	readFastaStream(inFile, processed, false);
  } else if (fileType == "flow") {
    readFlow(filename, processed, false);
  } else if (fileType == "qual") {
    readQual(filename, processed, false);
  } else if (fileType == "raw") {
    readRaw(filename, processed, false);
  } else if (fileType == "bam") {
    readBam(filename, processed);
  } else if (fileType == "pyroData") {
    readPyroData(filename, processed);
  } else if (fileType == "fastaQual") {
    readFastaQual(filename, secondName, processed);
  } else if (fileType == "fastqFlow") {
    readFastqFlow(filename, secondName, processed);
  } else if (fileType == "fastaFlow") {
    readFastaFlow(filename, secondName, processed);
  } else if (fileType == "flowQual") {
    readFlowQual(filename, secondName, processed);
  } else if (fileType == "fastaQualFlow") {
    readFastaQualFlow(filename, secondName, thirdName, processed);
  } else {
    std::cout << "Unrecognized file type : " << fileType << ", not reading "
              << filename << std::endl;
    std::cout << "Acceptable types are fasta, qual,fastq, stub, sff,  "
              << std::endl;
    exit(1);
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
		std::cout << "file: " << shorahFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
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
		std::cout << "file: " << shorahFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
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

void readObjectIO::readFlow(std::string filename, bool processed, bool add) {
  if (!add) {
    reads.clear();
  }
  FILE* dataFile;
  dataFile = fopen(filename.c_str(), "r");
  if (!dataFile) {
    std::cout << "Error in opening " << filename << std::endl;
    exit(1);
  }
  char c;
  c = fgetc(dataFile);
  std::string name = "";
  std::vector<double> flows;
  while (c != EOF) {
    if (c == '>' && c != EOF) {
      readNextData(dataFile, c, name, flows);
      if (add) {
        readVec::getReadByName(reads, name).flowValues = flows;
      } else {
        readObject nextRead = readObject(
            seqInfo(name, seqUtil::getSeqFromFlow(flows)), processed);
        nextRead.flowValues = flows;
        reads.emplace_back(nextRead);
      }
    }
  }
  fclose(dataFile);
}
void readObjectIO::readRaw(std::string filename, bool processed, bool add) {
  if (!add) {
    reads.clear();
  }
  FILE* dataFile;
  dataFile = fopen(filename.c_str(), "r");
  if (!dataFile) {
    std::cout << "Error in opening " << filename << std::endl;
    exit(1);
  }

  char c;
  c = fgetc(dataFile);
  std::string firstLine = "";
  while (c != '>') {
    firstLine.push_back(c);
    c = fgetc(dataFile);
  }
  // std::cout<<"firstLine: "<<firstLine<<std::endl;
  std::string name = "";
  std::vector<double> flows;
  while (c != EOF) {
    if (c == '>' && c != EOF) {
      readNextData(dataFile, c, name, flows);
      if (add) {
        readVec::getReadByName(reads, name).flowValues = flows;
      } else {
        readObject nextRead = readObject(
            seqInfo(name, seqUtil::getSeqFromFlow(flows)), processed);
        nextRead.flowValues = flows;
        reads.emplace_back(nextRead);
      }
    }
  }
  fclose(dataFile);
}

void readObjectIO::readQual(std::string filename, bool processed, bool add) {
  if (!add) {
    reads.clear();
  }
  FILE* qualFile;
  qualFile = fopen(filename.c_str(), "r");
  if (!qualFile) {
    std::cout << "Error in opening " << filename << std::endl;
    exit(1);
  }
  char c;
  c = fgetc(qualFile);
  std::string name = "";
  std::string qual = "";
  // std::vector<double> flows;
  while (c != EOF) {
    if (c == '>' && c != EOF) {
      readNextQual(qualFile, c, name, qual);
      if (add) {
        readVec::getReadByName(reads, name).addQual(qual);
      } else {
        reads.emplace_back(readObject(seqInfo(name, "", qual), processed));
      }
    }
  }
  fclose(qualFile);
}

void readObjectIO::readFastaQual(std::string fastaName,
                                 std::string qualFilename, bool processed) {
  reads.clear();
	if(!fexists(fastaName)){
		std::cout << "file: " << fastaName << ", doesn't exist or can't be open" << std::endl;
		exit(1);
	} else if (!fexists(qualFilename)){
		std::cout << "file: " << qualFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
	} else {
		std::ifstream inFile(fastaName);
	  readFastaStream(inFile, processed, false);
	  readQual(qualFilename, processed, true);
	}

}

void readObjectIO::readFastqFlow(std::string fastqFilename,
                                 std::string flowFilename, bool processed) {
  reads.clear();
	if(!fexists(fastqFilename)){
		std::cout << "file: " << fastqFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
	}else{
		std::ifstream inFile(fastqFilename);
	  readFastqStream(inFile,SangerQualOffset, processed, false);
	  readFlow(flowFilename, processed, true);
	}
}

void readObjectIO::readFastaFlow(std::string fastaFilename,
                                 std::string flowFilename, bool processed) {
  reads.clear();
	if(!fexists(fastaFilename)){
		std::cout << "file: " << fastaFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
	}else{
		std::ifstream inFile(fastaFilename);
	  readFastaStream(inFile, processed, false);
	  readFlow(flowFilename, processed, true);
	}
}

void readObjectIO::readFlowQual(std::string flowFilename,
                                std::string qualFilename, bool processed) {
  reads.clear();
  readFlow(flowFilename, processed, false);
  readQual(qualFilename, processed, true);
}

void readObjectIO::readPyroData(std::string filename, bool processed) {
  reads.clear();
  textFileReader txtReader = textFileReader();
  txtReader.readFile(filename);
  int count = 0;
  for (auto& fileIter : txtReader.fileContent.content_) {
    if (count == 0) {
      std::string name = fileIter[0];
      int numberOfFlows = atoi(fileIter[0].c_str());
      fileIter.erase(fileIter.begin());
      fileIter.erase(fileIter.begin());
      reads.emplace_back(readObject(seqInfo(
          name, seqUtil::getSeqFromFlow(
                    stringToVector<double>(vectorToString(fileIter))))));
      reads[reads.size() - 1].numberOfFlows = numberOfFlows;
      reads[reads.size() - 1].flowValues =
          stringToVector<double>(vectorToString(fileIter));
    }
    ++count;
  }
}
bool readObjectIO::refillBuffer(std::istream & is){
	lineNum_ = 0;
	while(lineNum_ < bufferMax_ && !is.eof() ){
		std::getline(is, lineBuffer_[lineNum_]);
		++lineNum_;
	}

	bufferPos_ = 0;
	if(lineNum_ == 0){
		//didn't read in anything, probably at end of the file
		return false;
	}else{
		// did read in some in the buffer
		return true;
	}
}
bool readObjectIO::readNextFastaStream(std::istream& is,
		readObject& read, bool processed) {
	bool fileNotEmpty = true;
	if(bufferPos_ == UINT32_MAX || bufferPos_ == bufferMax_ || bufferPos_ == lineNum_){
		// if the buffer is initialize yet or it reached the
		// end of the buffer read in the next ten line
		fileNotEmpty = refillBuffer(is);
	}
	if(!fileNotEmpty){
		return false;
	}
  std::string name = "";
  std::string buildingSeq = "";
  if(lineBuffer_[bufferPos_][0] =='>'){
  	name = lineBuffer_[bufferPos_];
  	//std::cout << boldText("name: ", "31") << name << std::endl;
  	++bufferPos_;
		if(bufferPos_ == bufferMax_){
			refillBuffer(is);
		}
  	while( bufferPos_ < lineNum_ && lineBuffer_[bufferPos_][0] !='>'){
			buildingSeq.append(lineBuffer_[bufferPos_]);
			++bufferPos_;
			if(bufferPos_ == bufferMax_){
				refillBuffer(is);
			}
			//std::cout << boldText("buildingSeq: ", "31") << buildingSeq << std::endl;
		}
  	read = readObject(seqInfo(name.substr(1), buildingSeq));
  	if(processed){
  		read.processRead(processed);
  	}
  	return true;
  }else{
  	std::cout << "error in reading fasta file, line doesn't begin with >" << std::endl;
  	std::cout << lineBuffer_[bufferPos_][0] << std::endl;
  	std::cout << "bufferPos " << bufferPos_ << std::endl;
  	printVector(lineBuffer_, "\n");
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
    std::cerr << bib::bashCT::bold
    		<< "Incomplete sequence, read only"
    		<< bib::bashCT::reset << std::endl;
    for (uint32_t i = 0; i < count; ++i) {
      std::cout << "!" << data[i] << "!" << std::endl;
    }
    std::cerr << bib::bashCT::bold
    		<< "exiting" << bib::bashCT::reset << std::endl;
    exit(1);
    return false;
  } else {
    return false;
  }
}

void readObjectIO::readFastaStream(std::istream& is, bool processed, bool add) {
	lineBuffer_.clear();
	if(!add){
		reads.clear();
	}
  readObject tempObj;
  while(readNextFastaStream(is, tempObj, processed)){
  	reads.emplace_back(tempObj);
  }
}

void readObjectIO::readClustal(std::string filename, bool processed) {
  reads.clear();
  textFileReader txtReader;
  txtReader.readFile(filename);
  std::vector<std::pair<std::string, readObject>> readMap;
  std::vector<std::pair<std::string, readObject>>::iterator readMapIter;
  for (const auto& fileIter : txtReader.fileContent.content_) {
    if (fileIter.size() != 2 || fileIter[0][0] == '*') {
    } else {
      bool foundMatch = false;
      for (readMapIter = readMap.begin(); readMapIter != readMap.end();
           ++readMapIter) {
        if (readMapIter->first == fileIter[0]) {
          foundMatch = true;
          readMapIter->second.seqBase_.seq_.append(fileIter[1]);
          break;
        }
      }
      if (!foundMatch) {
        readMap.emplace_back(std::make_pair(
            fileIter[0], readObject(seqInfo(fileIter[0], fileIter[1]))));
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
    std::cout << "Problem opening " << filename << std::endl;
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
  for (const auto& sIter : sffReads) {
    // readObject tempRead(sIter->seqBase_.name_,sIter->seq,sIter->qual);
    // tempRead.flowValues=sIter->flowValues;
    reads.emplace_back(convertSffObject(sIter));
    ++convertCount;
  }
}

void readObjectIO::readFastaQualFlow(std::string fastaFilename,
                                     std::string qualFilename,
                                     std::string flowFilename, bool processed) {
  reads.clear();
	if(!fexists(fastaFilename)){
		std::cout << "file: " << fastaFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
	} else if (!fexists(qualFilename)){
		std::cout << "file: " << qualFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
	} else if (!fexists(flowFilename)){
		std::cout << "file: " << flowFilename << ", doesn't exist or can't be open" << std::endl;
		exit(1);
	} else {
		std::ifstream inFile(fastaFilename);
	  readFastaStream(inFile, processed, false);
	  readQual(qualFilename, processed, true);
	  readFlow(flowFilename, processed, true);
	}

}

///////////////////////reading the files
void readObjectIO::readNextName(FILE* file, std::string& name, char& c) {
  while ((includeSpaceInNames || c != 9) && (includeSpaceInNames || c != 32) &&
         c != 10 && c != 13 && c != EOF) {
    name.push_back(c);
    c = fgetc(file);
  }
  while (c != 10 && c != 13 && c != EOF) {
    c = fgetc(file);
  }
}

void readObjectIO::readNextSeq(FILE* fastaFile, char& c, std::string& name,
                               std::string& seq) {
  name = "";
  seq = "";
  c = fgetc(fastaFile);
  readNextName(fastaFile, name, c);
  c = fgetc(fastaFile);
  while (c != '>' && c != EOF) {
    if (c == 10 || c == 13) {
    } else {
      seq.push_back(c);
    }
    c = fgetc(fastaFile);
  }
  // return readObject(name,seq,processed);
}

void readObjectIO::readNextQual(FILE* qualFile, char& c, std::string& name,
                                std::string& qual) {
  name = "";
  qual = "";
  c = fgetc(qualFile);
  readNextName(qualFile, name, c);
  /*
    while (c!=9 && c!=32 && c!=10 && c!=13 && c!=EOF) {
    name.push_back(c);
    c=fgetc(qualFile);
    }
    while (c!=10 && c!=13 && c!=EOF) {
    c=fgetc(qualFile);
    }*/

  c = fgetc(qualFile);
  while (c != '>' && c != EOF) {
    if (c == 10 || c == 13) {
      qual.push_back(' ');
    } else {
      qual.push_back(c);
    }
    c = fgetc(qualFile);
  }
  // return qual;
}

void readObjectIO::readNextData(FILE* dataFile, char& c, std::string& name,
                                std::vector<double>& flows) {
  name = "";
  std::string data = "";
  c = fgetc(dataFile);
  readNextName(dataFile, name, c);
  /*
    while (c!=9 && c!=32 && c!=10 && c!=13 && c!=EOF) {
    name.push_back(c);
    c=fgetc(dataFile);
    }
    while (c!=10 && c!=13 && c!=EOF) {
    c=fgetc(dataFile);
    }*/
  c = fgetc(dataFile);
  while (c != '>' && c != EOF) {
    if (c == 10 || c == 13) {
      data.push_back(' ');
    } else {
      data.push_back(c);
    }
    c = fgetc(dataFile);
  }
  flows = stringToVector<double>(data);
  // need to find a better way to do this checking, nick 7.29.2013
  if (flows[0] > 50) {
    flows.erase(flows.begin());
  }
  // return stringToVector<double>(data);
}

void readObjectIO::readNextFastq(FILE* file, char& c, std::string& name,
                                 std::string& seq, std::string& qual) {
  name = "";
  seq = "";
  qual = "";
  c = fgetc(file);
  // id before first white space
  readNextName(file, name, c);
  /*
    while (c!=9 && c!=32 && c!=10 && c!=13 && c!=EOF) {
    name.push_back(c);
    c=fgetc(file);
    }
    while (c!=10 && c!=13 && c!=EOF) {
    c=fgetc(file);
    }*/
  // seq
  c = fgetc(file);
  while (c != '+' && c != EOF) {
    if (c == 10 || c == 13) {
    } else {
      seq.push_back(c);
    }
    c = fgetc(file);
  }
  // qual
  c = fgetc(file);
  c = fgetc(file);
  while (c != 10 && c != 13 && c != EOF) {
    qual.push_back(c);
    c = fgetc(file);
  }
  if (c != EOF) {
    c = fgetc(file);
  }
  // return readObject(name,seq, qual, 33, processed);
}

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
		uint64_t startSpotInFile = in.tellg();

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
			std::cout << "Error in reading sff binary object" << std::endl;
		}

    if (in.eof()) {
    	return true;
    }

		return false;
	}
	catch(std::exception& e) {
		std::cout << "Error in reading during readSeqData" << std::endl;
		exit(1);
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
				std::cout << "Error reading sff common header." << std::endl;
		}
	}
	catch(std::exception& e) {
		//m->errorOut(e, "SffInfoCommand", "readCommonHeader");
		exit(1);
	}
}

void readObjectIO::readSffbin(const std::string & filename){
	std::ifstream inFile;
	if(!fexists(filename)){
		std::cout << "File: " << filename << " doesn't exist" << std::endl;
		exit(1);
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
