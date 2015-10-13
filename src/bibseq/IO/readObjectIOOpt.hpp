#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include <api/BamReader.h>
#include <bibcpp/jsonUtils.h>
#include "bibseq/IO/readObjectIOOptions.hpp"
#include "bibseq/IO/cachedReader.hpp"
#include "bibseq/objects/seqObjects/readObject.hpp"
#include "bibseq/objects/seqObjects/sffObject.hpp"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace bibseq {

class readObjectIOOpt {
public:
	readObjectIOOpt(const readObjectIOOptions & options);
  static const uint32_t IlluminaQualOffset;
  static const uint32_t SangerQualOffset;
  const readObjectIOOptions ioOptions_;
  bool includeSpaceInNames = true;
  bool readNextRead(readObject & read);
  bool readNextRead(std::shared_ptr<readObject> & read);
  bool readNextRead(std::unique_ptr<readObject> & read);
  std::shared_ptr<readObject> readNextRead();


  bool inOpen()const;
  void openIn();
  void closeIn();

  bool outOpen()const;
  void openOut();
  void closeOut();

  template<typename T>
  void write(const T & read){
  	if(!outOpen_){
  		throw std::runtime_error{"Error in readObjectIOOpt, attempted to write when out files aren't open"};
  	}
  	if(ioOptions_.outFormat_ == "fastaQual"){
  		read.seqBase_.outPut(*fOut_, *qOut_, ioOptions_);
  	}else if(ioOptions_.outFormat_ == "mothurData"){
  		read.outPutPyroData(*fOut_);
  	}else{
  		read.seqBase_.outPut(*fOut_, ioOptions_);
  	}
  }

private:

  BamTools::BamReader bReader_;
  BamTools::BamAlignment aln_;

  std::unique_ptr<cachedReader> fReader_;
  std::unique_ptr<cachedReader> qReader_;
  std::unique_ptr<std::ifstream> fStreamReader_;
  std::unique_ptr<std::ifstream> qStreamReader_;

  VecStr sffTxtHeader_;

  std::unique_ptr<boost::iostreams::filtering_streambuf<boost::iostreams::input>> fqInputBuffer_;
  std::unique_ptr<std::ifstream> fqReaderBackup_;
  std::unique_ptr<std::istream> fqReader_;
  bool inOpen_ = false;
  bool outOpen_ = false;
  std::unique_ptr<std::ofstream> fOut_;
  std::unique_ptr<std::ofstream> qOut_;
  bool readNextFastaStream(cachedReader& cReader, readObject& read, bool processed);
  bool readNextFastaQualStream(cachedReader& fastaReader, cachedReader& qualReader, readObject & read,bool processed);
  bool readNextQualStream(cachedReader& cReader, std::vector<uint32_t>& quals, std::string & name);
  bool readNextFastqStream(std::istream& is, uint32_t offSet, readObject& read, bool processed);
  bool readNextBam(BamTools::BamReader & bReader, readObject& read,
  		BamTools::BamAlignment & aln, bool processed);

  bool readNextSff(std::ifstream & inFile, sffObject & read);
  bool readNextSff(std::ifstream & inFile, readObject & read);

  VecStr readSffTxtHeader(std::ifstream & inFile);


};


}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "readObjectIOOpt.cpp"
#endif
