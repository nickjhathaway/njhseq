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
#include "seqToolsUtils.hpp"
#include "bibseq/objects/seqObjects/cluster.hpp"
#include "bibseq/seqToolsUtils/distCalc.hpp"

//////tools for dealing with multipleSampleCollapse
namespace bibseq {





void getReadCnt(const std::string & filename, std::string & format,
		bool processed, uint64_t & currentCount) {
	readObjectIO reader;
	std::ifstream inFile(filename);
	if (!inFile) {
		std::stringstream ss;
		ss << bib::bashCT::red << bib::bashCT::bold
				<< "Error in opening file: " << filename << bib::bashCT::reset;
		throw std::runtime_error{ss.str()};
	}
	if (format == "fastq") {
		readObject read;
		while (reader.readNextFastqStream(inFile, readObjectIO::SangerQualOffset,
				read, processed)) {
			currentCount += read.seqBase_.cnt_;
		}
	} else if (format == "bam") {
		BamTools::BamAlignment aln;
		BamTools::BamReader bReader;
		readObject read;
		bReader.Open(filename);
		if (!bReader.IsOpen()) {
			std::stringstream ss;
			ss << bib::bashCT::red << bib::bashCT::bold
					<< "Error in opening file: " << filename << bib::bashCT::reset;
			throw std::runtime_error{ss.str()};
		}
		while (reader.readNextBam(bReader, read, aln, processed)) {
			currentCount += read.seqBase_.cnt_;
		}
	} else if (format == "fasta") {
		readObject read;
		cachedReader cinFile(inFile);
		while (reader.readNextFastaStream(cinFile, read, processed)) {
			currentCount += read.seqBase_.cnt_;
		}
	} else {
		reader.read(format, filename, processed);
		for (const auto & read : reader.reads) {
			currentCount += read.seqBase_.cnt_;
		}
	}
}




void setUpSampleDirs(
    const std::string& barcodeFilename,
		const std::string& mainDirectoryName) {
	auto topDir = bib::replaceString(mainDirectoryName, "./", "");
	table inTab(barcodeFilename, "whitespace", false);
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> sampleDirWithSubDirs;
	for(const auto & rowPos : iter::range(inTab.content_.size())){
		const auto & row = inTab.content_[rowPos];
		if(row.empty() || row[0].front() == '#'){
			continue;
		}
		if (row.size() < 3) {
			throw std::runtime_error { bib::err::F()
					<< "setUpSampleDirs: rows should have at least 3 columns, row: "
					<< rowPos << "has " << row.size() };
		}
		for(const auto & colPos : iter::range<uint32_t>(2,row.size())){
			if(row[colPos] == "" || allWhiteSpaceStr(row[colPos])){
				continue;
			}
			sampleDirWithSubDirs[row[1]].emplace_back(row[colPos], row[0]);
		}
	}
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> indexsWithFullSampPathNames;
	auto cwd = get_cwd();
	for(const auto & sDirs : sampleDirWithSubDirs){
		std::cout << bib::bashCT::bold << "Making Sample Dir: " << bib::bashCT::green << topDir + sDirs.first
				<< bib::bashCT::reset << std::endl;
		std::string sampDir = "";
		try {
			sampDir = bib::files::makeDir(topDir, sDirs.first);
		} catch (std::exception & e) {
			std::stringstream ss;
			ss << bib::bashCT::boldRed(e.what()) << std::endl;
			throw std::runtime_error{ss.str()};
		}
		for(const auto & midNames : sDirs.second){
			std::cout << bib::bashCT::bold << "Making Mid Dir: " << bib::bashCT::blue << sampDir + midNames.first
					<< bib::bashCT::reset << std::endl;
			std::string midDir = "";
			try {
				midDir = bib::files::makeDir(sampDir,midNames.first);
			} catch (std::exception & e) {
				std::stringstream ss;
				std::cerr << bib::bashCT::boldRed(e.what()) << std::endl;
				throw std::runtime_error{ss.str()};
			}
			bib::files::appendAsNeeded(cwd, "/");
			indexsWithFullSampPathNames[midNames.second].emplace_back(midNames.first,cwd + midDir);
		}
	}

  std::string indexDir = bib::files::makeDir(mainDirectoryName, "locationByIndex");
  for (const auto& indexIter : indexsWithFullSampPathNames) {
    std::ofstream indexFile;
    openTextFile(indexFile, indexDir + replaceString(indexIter.first),
                 ".tab.txt", false, false);
    for (const auto& spIter : indexIter.second) {
      indexFile << spIter.first << "\t"
                << replaceString(spIter.second, "./", "") << std::endl;
    }
  }
}
std::string genHtmlStrForPsuedoMintree(std::string jsonFileName){
	std::string ret = "<!DOCTYPE html>"
	"<meta charset=\"utf-8\">"
	"<body>"
	"<script src=\"http://d3js.org/d3.v3.min.js\"></script>"
	"<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js\"></script>"
	"<script src=\"http://bib4.umassmed.edu/~hathawan/js/psuedoMinTree.js\"></script>"
	"<button id=\"save_svg\">Save as Svg</button>"
	"<svg id = \"main\"></svg>"
	"<script>"
	"var jsonDat = \""+ jsonFileName + "\";"
	"var add = \"#main\";"
	"drawPsuedoMinTree(jsonDat, add);"
	"var bName = \"#save_svg\";"
	"addSvgSaveButton(bName, add);"
	"</script>"
	"</body>";
	return ret;
}
simulation::errorProfile getErrors(const simulation::profile & allProfile,
		const std::vector<char> & alphabet){
	simulation::errorProfile eProfile(alphabet);
  for(auto & firstChar : alphabet){
  	for(auto & secondChar : alphabet){
  		if(firstChar != secondChar){
  			eProfile.increaseCountAmount(firstChar, secondChar, allProfile.counts_[firstChar][secondChar]);
  		}
  	}
  }
  eProfile.setFractions();
  return eProfile;
}


substituteMatrix getMatrixFromProfile(simulation::profile allProfile, const std::vector<char> & alphabet){
	std::stringstream oddsOut;
	allProfile.printLogOddsMatrix(alphabet, true, oddsOut);
	substituteMatrix ans (seqUtil::readScoringMatrix(oddsOut));
	return ans;
}




void processRunCutoff(uint32_t& runCutOff, const std::string& runCutOffString,
                      int counter) {
	auto toks = tokenizeString(runCutOffString, ",");
	if(toks.size() == 1){
	  if (runCutOffString.back() == '%') {
	    runCutOff = std::round(
	        std::stof(runCutOffString.substr(0, runCutOffString.length() - 1)) *
	        counter / 100);
	  } else {
	    runCutOff = std::stoi(runCutOffString);
	  }
	}else{
	  if (toks[0].back() == '%') {
	    runCutOff = std::round(
	        std::stod(toks[0].substr(0, toks[0].length() - 1)) *
	        counter / 100.0);
	  } else {
	    runCutOff = std::stoi(toks[0]);
	  }
	  int32_t hardCutOff = std::stoi(toks[1]);
	  if(hardCutOff > runCutOff){
	  	runCutOff = hardCutOff;
	  }
	}
}

uint32_t processRunCutoff(const std::string& runCutOffString,
                      uint64_t counter){
	uint32_t ret = 0;
	processRunCutoff(ret, runCutOffString, counter);
  return ret;
}







/////////tools for finding additional location output
std::string makeIDNameComparable(const std::string& idName) {
  std::string ans = replaceString(idName, "ID0", "");
  ans = replaceString(ans, "M0", "M");
  return replaceString(ans, "ID", "");
}
std::string findIdNameFromFileName(const std::string& filename) {
  auto periodPos = filename.rfind(".");
  if (periodPos != std::string::npos && !isdigit(filename[periodPos - 1])) {
    periodPos = filename.rfind("_");
  }
  auto lastMPos = filename.rfind("M");
  if(lastMPos == std::string::npos){
  	return filename;
  }
  return filename.substr(lastMPos, periodPos - lastMPos);
}
std::string processFileNameForID(const std::string& fileName) {
  return makeIDNameComparable(findIdNameFromFileName(fileName));
}

std::string findAdditonalOutLocation(const std::string& locationFile,
                                     const std::string& fileName) {
	table inTab(locationFile, "\t");
  MapStrStr additionalOutNames;
  for (const auto& fIter : inTab.content_) {
    additionalOutNames[makeIDNameComparable(fIter[0])] = fIter[1];
  }

  return additionalOutNames[processFileNameForID(fileName)];
}




void processAlnInfoInput(aligner& alignerObj,
                         const std::string& alnInfoDirName) {
  if (alnInfoDirName != "") {
  	auto fullPath = bib::files::normalize(alnInfoDirName).string();
  	auto search = alignment::alnCacheDirLocks.find(fullPath);
  	if(search == alignment::alnCacheDirLocks.end()){
  		alignment::alnCacheDirLocks.emplace(std::piecewise_construct, std::make_tuple(fullPath), std::tuple<>());
  	}
  	auto realSearch = alignment::alnCacheDirLocks.find(fullPath);
  	std::lock_guard<std::mutex> lock(realSearch->second);

    int directoryStatus =
        mkdir(alnInfoDirName.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    if (directoryStatus != 0) {
    	auto readInHoler = alnInfoMasterHolder(alnInfoDirName, alignerObj.parts_.gapScores_,
    			alignerObj.parts_.scoring_);

    	//add local alignments
    	//std::cout << "Read in previous alns infos " << std::endl;
    	for(const auto & lHolder : readInHoler.localHolder_){
    		//std::cout << lHolder.first << std::endl;
    		alignerObj.alnHolder_.localHolder_[lHolder.first] = lHolder.second;
    	}
    	//add global alignments
    	for(const auto & gHolder : readInHoler.globalHolder_){
    		//std::cout << gHolder.first << std::endl;
    		alignerObj.alnHolder_.globalHolder_[gHolder.first] = gHolder.second;
    	}
      //alignerObj.alnHolders_ = alnInfoMasterHolder(alnInfoDirName);
    } else {
      chmod(alnInfoDirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
  }
}




std::vector<uint32_t> getWindowQuals(const std::vector<uint32_t>& qual,
                                     uint32_t qWindowSize, uint32_t pos) {
  uint32_t lowerBound = 0;
  uint32_t higherBound = qual.size();
  if (pos > qWindowSize) {
    lowerBound = pos - qWindowSize;
  }
  if (pos + qWindowSize + 1 < higherBound) {
    higherBound = pos + qWindowSize + 1;
  }
  return getSubVector(qual, lowerBound, higherBound - lowerBound);
}
std::vector<double> likelihoodForBaseQ(
    const std::vector<uint32_t>& qual,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    ans.emplace_back(likelihoods.at(qual[i]));
  }
  return ans;
}
std::vector<double> likelihoodForMeanQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    double currentMean = vectorMean(getWindowQuals(qual, qWindowSize, i));
    ans.emplace_back(likelihoods.at(roundDecPlaces(currentMean, 2)));
  }
  return ans;
}
std::vector<double> likelihoodForMedianQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    double currentMedian = vectorMedian(getWindowQuals(qual, qWindowSize, i));
    ans.emplace_back(likelihoods.at(roundDecPlaces(currentMedian, 2)));
  }
  return ans;
}
std::vector<double> likelihoodForMinQ(
    const std::vector<uint32_t>& qual, uint32_t qWindowSize,
    std::unordered_map<double, double>& likelihoods) {
  std::vector<double> ans;
  for (const auto& i : iter::range(qual.size())) {
    double currentMin = vectorMinimum(getWindowQuals(qual, qWindowSize, i));
    ans.emplace_back(likelihoods.at(roundDecPlaces(currentMin, 2)));
  }
  return ans;
}

cluster createClusterFromSimiliarReads(std::vector<readObject> similiarReads,
                                       aligner& alignerObj) {
  if (similiarReads.size() > 0) {
    sort(similiarReads);
    cluster output(similiarReads.front());
    for (const auto& readPos : iter::range<uint32_t>(1, similiarReads.size())) {
      output.addRead(similiarReads[readPos]);
    }
    output.calculateConsensus(alignerObj, true);
    return output;
  } else {
    return cluster();
  }
}





table getSeqPosTab(const std::string & str){
	table out(VecStr{"char", "pos"});
	for(const auto & pos : iter::range(str.size())){
		out.content_.emplace_back(toVecStr(str[pos], pos));
	}
	return out;
}

uint32_t processCutOffStr(const std::string& runCutOffString,
  uint64_t readCount){
	uint32_t runCutOff;
  if (runCutOffString.back() == '%') {
    runCutOff = std::round(
        std::stod(runCutOffString.substr(0, runCutOffString.length() - 1)) *
        readCount / 100.0);
  } else {
    runCutOff = std::stoi(runCutOffString);
  }
  return runCutOff;
}




uint32_t getAlnPos(const std::string & seq, uint32_t realSeqPos) {
  uint32_t offSet = 0;
  for (uint32_t i = 0; i < seq.size(); i++) {
    if ((i - offSet) == realSeqPos) {
      return i;
    }
    if (seq[i] == '-') {
      ++offSet;
    }
  }
  return seq.size();
}

uint32_t getRealPos(const std::string & seq, uint32_t seqAlnPos) {
  uint32_t offSet = 0;
  for (uint32_t i = 0; i < seqAlnPos; i++) {
    if (seq[i] == '-') {
      ++offSet;
    }
  }
  return seqAlnPos - offSet;
}
}  // namespace bib
