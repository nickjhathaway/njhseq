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

#include "bibseq/IO/SeqIO.h"
#include "seqToolsUtils.hpp"
#include "bibseq/objects/seqObjects/Clusters/cluster.hpp"
#include "bibseq/seqToolsUtils/distCalc.hpp"
#include "bibseq/objects/kmer/kmerCalculator.hpp"
#include "bibseq/objects/dataContainers/graphs/UndirWeightedGraph.hpp"


namespace bibseq {

void processRunCutoff(uint32_t& runCutOff, const std::string& runCutOffString,
		int counter) {
	auto toks = tokenizeString(runCutOffString, ",");
	if (toks.size() == 1) {
		if (runCutOffString.back() == '%') {
			runCutOff = std::round(
					std::stof(runCutOffString.substr(0, runCutOffString.length() - 1))
							* counter / 100);
		} else {
			runCutOff = estd::stou(runCutOffString);
		}
	} else {
		if (toks[0].back() == '%') {
			runCutOff = std::round(
					std::stod(toks[0].substr(0, toks[0].length() - 1)) * counter / 100.0);
		} else {
			runCutOff = estd::stou(toks[0]);
		}
		uint32_t hardCutOff = estd::stou(toks[1]);
		if (hardCutOff > runCutOff) {
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



std::string genHtmlStrForPsuedoMintree(std::string jsonFileName){
	std::string ret = "<!DOCTYPE html>"
	"<meta charset=\"utf-8\">"
	"<body>"
	"<script src=\"http://d3js.org/d3.v3.min.js\"></script>"
	"<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js\"></script>"
	"<script src=\"http://bib8.umassmed.edu/~hathawan/js/psuedoMinTree.js\"></script>"
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

std::string genHtmlStrForPsuedoMintree(std::string jsonFileName, std::string javaScriptFile){
	std::string ret = "<!DOCTYPE html>"
	"<meta charset=\"utf-8\">"
	"<body>"
	"<script src=\"http://d3js.org/d3.v3.min.js\"></script>"
	"<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js\"></script>"
	"<script src=\"" + javaScriptFile + "\"></script>"
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

Json::Value genMinTreeData(const std::vector<readObject> & reads,
		aligner & alignerObj,
		std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
		std::mutex & alignerLock, uint32_t numThreads) {

	std::function<
			uint32_t(const readObject &, const readObject &,
					std::unordered_map<std::string, std::unique_ptr<aligner>>&, aligner&)> getMismatchesFunc =
			[&alignerLock](const readObject & read1, const readObject & read2,std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners, aligner &alignerObj) {
				alignerLock.lock();
				auto threadId = estd::to_string(std::this_thread::get_id());
				//std::cout << threadId<< std::endl;
				if(aligners.find(threadId) == aligners.end()) {
					aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
				}
				alignerLock.unlock();
				aligners.at(threadId)->alignCache(read1.seqBase_,read2.seqBase_, false);
				aligners.at(threadId)->profilePrimerAlignment(read1.seqBase_, read2.seqBase_);
				return aligners.at(threadId)->comp_.hqMismatches_;
			};
	auto distances = getDistanceNonConst(reads, numThreads, getMismatchesFunc,
			aligners, alignerObj);
	readDistGraph<uint32_t> graphMis(distances, reads);
	std::vector<std::string> popNames;
	for (const auto & n : graphMis.nodes_) {
		popNames.emplace_back(n->name_);
	}
	auto nameColors = getColorsForNames(popNames);
	auto ret = graphMis.toJsonMismatchGraphAll(bib::color("#000000"), nameColors);
	return ret;
}


Json::Value genMinTreeData(const std::vector<readObject> & reads, aligner & alignerObj){
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	uint32_t numThreads = 2;
	std::mutex alignerLock;
	std::function<uint32_t(const readObject &, const readObject &,std::unordered_map<std::string, std::unique_ptr<aligner>>& , aligner& )> getMismatchesFunc =
			[&alignerLock](const readObject & read1, const readObject & read2,std::unordered_map<std::string, std::unique_ptr<aligner>>&  aligners, aligner &alignerObj) {
		alignerLock.lock();
		auto threadId = estd::to_string(std::this_thread::get_id());
		//std::cout << threadId<< std::endl;
		if(aligners.find(threadId) == aligners.end()){
			aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
		}
		alignerLock.unlock();
		aligners.at(threadId)->alignCache(read1.seqBase_,read2.seqBase_, false);
		aligners.at(threadId)->profilePrimerAlignment(read1.seqBase_, read2.seqBase_);
				return aligners.at(threadId)->comp_.hqMismatches_;
			};

	auto distances = getDistanceNonConst(reads, numThreads, getMismatchesFunc, aligners, alignerObj);
	readDistGraph<uint32_t> graphMis(distances, reads);
	std::vector<std::string> popNames;
	for (const auto & n : graphMis.nodes_) {
		popNames.emplace_back(n->name_);
	}
	bool debug = false;
	auto nameColors = getColorsForNames(popNames);
	for(const auto & alnObj : aligners){
		if(debug){
			std::cout << alnObj.first << ": " << alnObj.second->numberOfAlingmentsDone_ << std::endl;
		}
		if(alnObj.second->numberOfAlingmentsDone_ > 0){
			alignerObj.alnHolder_.mergeOtherHolder(alnObj.second->alnHolder_);
		}
	}
	return graphMis.toJsonMismatchGraphAll(bib::color("#000000"), nameColors);
}

Json::Value genMinTreeData(const std::vector<readObject> & reads) {
	uint64_t maxLength = 0;
	readVec::getMaxLength(reads, maxLength);
	aligner alignerObj(maxLength, gapScoringParameters(5, 1),
			substituteMatrix::createDegenScoreMatrix(2, -2));
	std::function<uint32_t(const readObject &, const readObject &, aligner)> misFun =
			getMismatches<readObject>;
	return genMinTreeData(reads, alignerObj);
}

uint32_t countSeqs(const SeqIOOptions & opts, bool verbose) {
	uint32_t ret = 0;
	if (bib::files::bfs::exists(opts.firstName_)) {
		SeqIO reader(opts);
		reader.openIn();
		if(!reader.in_.inOpen()){
			std::stringstream ss;
			ss << "Error in " << __PRETTY_FUNCTION__ << " in opening file with options" << std::endl;
			ss << opts.toJson() << std::endl;
			throw std::runtime_error{ss.str()};
		}
		if(SeqIOOptions::inFormats::FASTQPAIRED == opts.inFormat_){
			PairedRead read;
			while (reader.readNextRead(read)) {
				ret += ::round(read.seqBase_.cnt_);
			}
		}else{
			seqInfo read;
			while (reader.readNextRead(read)) {
				ret += ::round(read.cnt_);
			}
		}

	} else if (verbose) {
		std::cout << "File: " << opts.firstName_ << "doesn't exist, returning 0"
				<< std::endl;
	}
	return ret;
}

ExtractionInfo collectExtractionInfo(const std::string & dirName, const std::string & indexToDir, const std::string & sampNames){
	std::string nameDelim = "_extractor_";
	auto dirs = getNewestDirs(dirName, nameDelim);
	table indexNamesTab(indexToDir,"\t", true);
  table mainTableExtractionProfile;
  table mainTableExtractionStats;
	std::unordered_map<std::string, std::string> nameToIndex;
	for (const auto & row : indexNamesTab.content_) {
		if (row.size() != 2) {
			std::cerr << "Error in parsing " << indexToDir
					<< ", should have two columns, not " << row.size() << std::endl;
			std::cerr << vectorToString(row, "\t") << std::endl;
			exit(1);
		}
		auto lastPeriod = row[1].rfind(".");
		nameToIndex[row[1].substr(0, lastPeriod)] = row[0];
	}
	uint32_t count = 0;
	VecStr oldProfileColNames;

	VecStr oldStatsColNames;
	table inTab(sampNames, "whitespace", false);
	//goes sample name -> pairs of midname - index name
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
			//Discrepancy  between people naming MID01 and MID1 so replacing MID0 with MID
			sampleDirWithSubDirs[row[1]].emplace_back(bib::replaceString(row[colPos], "MID0", "MID"), row[0]);
		}
	}

	//extraction info by index by mid;
	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> extractionInfo;

	for (const auto &dir : dirs) {
		table inProfileTab(dir + "/extractionProfile.tab.txt", "\t", true);
		//i have accidentally added one more tab than is needed to
		//extractionProfile which makes it so all rows have an empty
		//at the end that correspond to any column names
		for (auto & row : inProfileTab.content_) {
			if (row.back() == "") {
				row.erase(row.begin() + row.size() - 1);
			}
		}
		table inStatsTab(dir + "/extractionStats.tab.txt", "\t", true);
		std::string dirName = bib::files::getFileName(dir);
		auto pos = dirName.find(nameDelim);
		std::string indexName = nameToIndex[dirName.substr(0, pos)];
		if (count == 0) {
			oldStatsColNames = inStatsTab.columnNames_;
			oldProfileColNames = inProfileTab.columnNames_;
		}
		inStatsTab.addColumn(VecStr { indexName }, "IndexName");
		inProfileTab.addColumn(VecStr { indexName }, "IndexName");
		for(const auto & row : inProfileTab.content_){
			auto midName = row[inProfileTab.getColPos("name")];
			auto midPos = midName.rfind("MID");
			midName = bib::replaceString(midName.substr(midPos), "MID0", "MID");
			extractionInfo[row[inProfileTab.getColPos("IndexName")]][midName] = row;
		}
		if (count == 0) {
			mainTableExtractionStats = inStatsTab;
			mainTableExtractionProfile = inProfileTab;
		} else {
			mainTableExtractionStats.rbind(inStatsTab, false);
			mainTableExtractionProfile.rbind(inProfileTab, false);
		}
		++count;
	}

	table outSampleInfo(catenateVectors(VecStr{"Sample"}, mainTableExtractionProfile.columnNames_));
	for(const auto & samp : sampleDirWithSubDirs){
		for(const auto & indexMidNames : samp.second){
			auto addingRow = extractionInfo[indexMidNames.second][indexMidNames.first];
			if(addingRow.empty()){
				addingRow = std::vector<std::string>(mainTableExtractionProfile.content_.front().size(), "0");
			}
			addingRow[mainTableExtractionProfile.getColPos("IndexName")] = indexMidNames.second;
			addingRow[mainTableExtractionProfile.getColPos("name")] = indexMidNames.first;
			outSampleInfo.content_.emplace_back(catenateVectors(VecStr{samp.first}, addingRow));
		}
	}

	auto outSampleColName = catenateVectors(VecStr{"Sample"}, catenateVectors(VecStr{"IndexName"}, oldProfileColNames));
	auto profileColNames = catenateVectors(VecStr{"IndexName"}, oldProfileColNames);
	auto statsColName = catenateVectors(VecStr{"IndexName"}, oldStatsColNames);


	outSampleInfo = outSampleInfo.getColumns(outSampleColName);
	outSampleInfo.sortTable("Sample", false);
	mainTableExtractionProfile = mainTableExtractionProfile.getColumns(profileColNames);
	mainTableExtractionStats = mainTableExtractionStats.getColumns(statsColName);
	outSampleInfo.columnNames_[outSampleInfo.getColPos("name")] = "MidName";
	outSampleInfo.setColNamePositions();
	return ExtractionInfo(mainTableExtractionStats, mainTableExtractionProfile, outSampleInfo);
}

ExtractionInfo collectExtractionInfoDirectName(const std::string & dirName, const std::string & indexToDir, const std::string & sampNames){
	std::string nameDelim = "_extractor_";
	VecStr dirs;
	table indexNamesTab(indexToDir,"\t", true);
  table mainTableExtractionProfile;
  table mainTableExtractionStats;
	std::unordered_map<std::string, std::string> nameToIndex;
	for (const auto & row : indexNamesTab.content_) {
		if (row.size() != 2) {
			std::cerr << "Error in parsing " << indexToDir
					<< ", should have two columns, not " << row.size() << std::endl;
			std::cerr << vectorToString(row, "\t") << std::endl;
			exit(1);
		}
		dirs.emplace_back(row[1]);
		nameToIndex[row[1]] = row[0];
	}
	uint32_t count = 0;
	VecStr oldProfileColNames;

	VecStr oldStatsColNames;
	table inTab(sampNames, "whitespace", false);
	//goes sample name -> pairs of midname - index name
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
			//Discrepancy  between people naming MID01 and MID1 so replacing MID0 with MID
			sampleDirWithSubDirs[row[1]].emplace_back(bib::replaceString(row[colPos], "MID0", "MID"), row[0]);
		}
	}

	//extraction info by index by mid;
	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> extractionInfo;

	for (const auto &dir : dirs) {
		auto fullDirName = bib::appendAsNeededRet(dirName, "/") + dir;
		table inProfileTab(fullDirName + "/extractionProfile.tab.txt", "\t", true);
		//i have accidentally added one more tab than is needed to
		//extractionProfile which makes it so all rows have an empty
		//at the end that correspond to any column names
		for (auto & row : inProfileTab.content_) {
			if (row.back() == "") {
				row.erase(row.begin() + row.size() - 1);
			}
		}
		table inStatsTab(fullDirName + "/extractionStats.tab.txt", "\t", true);
		std::string indexName = nameToIndex[dir];
		if (count == 0) {
			oldStatsColNames = inStatsTab.columnNames_;
			oldProfileColNames = inProfileTab.columnNames_;
		}
		inStatsTab.addColumn(VecStr { indexName }, "IndexName");
		inProfileTab.addColumn(VecStr { indexName }, "IndexName");
		for(const auto & row : inProfileTab.content_){
			auto midName = row[inProfileTab.getColPos("name")];
			auto midPos = midName.rfind("MID");
			midName = bib::replaceString(midName.substr(midPos), "MID0", "MID");
			extractionInfo[row[inProfileTab.getColPos("IndexName")]][midName] = row;
		}
		if (count == 0) {
			mainTableExtractionStats = inStatsTab;
			mainTableExtractionProfile = inProfileTab;
		} else {
			mainTableExtractionStats.rbind(inStatsTab, false);
			mainTableExtractionProfile.rbind(inProfileTab, false);
		}
		++count;
	}

	table outSampleInfo(catenateVectors(VecStr{"Sample"}, mainTableExtractionProfile.columnNames_));
	for(const auto & samp : sampleDirWithSubDirs){
		for(const auto & indexMidNames : samp.second){
			auto addingRow = extractionInfo[indexMidNames.second][indexMidNames.first];
			if(addingRow.empty()){
				addingRow = std::vector<std::string>(mainTableExtractionProfile.content_.front().size(), "0");
			}
			addingRow[mainTableExtractionProfile.getColPos("IndexName")] = indexMidNames.second;
			addingRow[mainTableExtractionProfile.getColPos("name")] = indexMidNames.first;
			outSampleInfo.content_.emplace_back(catenateVectors(VecStr{samp.first}, addingRow));
		}
	}

	auto outSampleColName = catenateVectors(VecStr{"Sample"}, catenateVectors(VecStr{"IndexName"}, oldProfileColNames));
	auto profileColNames = catenateVectors(VecStr{"IndexName"}, oldProfileColNames);
	auto statsColName = catenateVectors(VecStr{"IndexName"}, oldStatsColNames);


	outSampleInfo = outSampleInfo.getColumns(outSampleColName);
	outSampleInfo.sortTable("Sample", false);
	mainTableExtractionProfile = mainTableExtractionProfile.getColumns(profileColNames);
	mainTableExtractionStats = mainTableExtractionStats.getColumns(statsColName);
	outSampleInfo.columnNames_[outSampleInfo.getColPos("name")] = "MidName";
	outSampleInfo.setColNamePositions();
	return ExtractionInfo(mainTableExtractionStats, mainTableExtractionProfile, outSampleInfo);
}





void setUpSampleDirs(
		const std::string& sampleNamesFilename,
		const std::string& mainDirectoryName,
		bool separatedDirs,
		bool verbose) {
	auto topDir = bib::replaceString(mainDirectoryName, "./", "");

	bib::appendAsNeeded(topDir, "/");
	bib::files::makeDirP(bib::files::MkdirPar(topDir));
	table inTab(sampleNamesFilename, "whitespace", false);

	//first key is target/index, second key is samp, value is vector of rep names and then the full path name for that rep
	std::unordered_map<std::string,
			std::unordered_map<std::string,
					std::vector<std::pair<std::string, std::string>>> >sampleDirWithSubDirs;
	for (const auto & rowPos : iter::range(inTab.content_.size())) {
		const auto & row = inTab.content_[rowPos];
		if (row.empty() || row[0].front() == '#') {
			continue;
		}
		if (row.size() < 3) {
			throw std::runtime_error { bib::err::F() << __PRETTY_FUNCTION__
					<< ": rows should have at least 3 columns, row: " << rowPos << "has "
					<< row.size() };
		}
		VecStr repNames;
		for (const auto & colPos : iter::range<uint32_t>(2, row.size())) {
			if (row[colPos] == "" || allWhiteSpaceStr(row[colPos])) {
				continue;
			}
			repNames.emplace_back(row[colPos]);
		}
		for (const auto & rep : repNames) {
			sampleDirWithSubDirs[row[0]][row[1]].emplace_back(
					std::pair<std::string, std::string> { rep, "" });
		}
	}
	auto cwd = bib::files::get_cwd();
	try {
		if (separatedDirs) {
			for (auto & targetDirs : sampleDirWithSubDirs) {
				if (verbose) {
					std::cout << bib::bashCT::bold << "Making Target Dir: "
							<< bib::bashCT::purple << topDir + targetDirs.first
							<< bib::bashCT::reset << std::endl;
				}
				std::string targetDir = bib::files::makeDir(topDir,
						bib::files::MkdirPar(targetDirs.first, false));
				for (auto & sampDirs : targetDirs.second) {
					if (verbose) {
						std::cout << bib::bashCT::bold << "Making Samp Dir: "
								<< bib::bashCT::green << targetDir + sampDirs.first
								<< bib::bashCT::reset << std::endl;
					}

					std::string sampDir = bib::files::makeDir(targetDir,
							bib::files::MkdirPar(sampDirs.first, false));
					for (auto & rep : sampDirs.second) {
						if (verbose) {
							std::cout << bib::bashCT::bold << "Making Rep Dir: "
									<< bib::bashCT::blue << sampDir + rep.first
									<< bib::bashCT::reset << std::endl;
						}

						std::string repDir = bib::files::makeDir(sampDir,
								bib::files::MkdirPar(rep.first, false));
						rep.second = bib::files::join(cwd, repDir);
					}
				}
			}
		} else {
			for (auto & targetDirs : sampleDirWithSubDirs) {
				for (auto & sampDirs : targetDirs.second) {
					std::string sampDir = bib::files::join(topDir, sampDirs.first);
					if (!bib::files::bfs::exists(sampDir)) {
						if (verbose) {
							std::cout << bib::bashCT::bold << "Making Samp Dir: "
									<< bib::bashCT::green << topDir + sampDirs.first
									<< bib::bashCT::reset << std::endl;
						}

						bib::files::makeDir(topDir,
								bib::files::MkdirPar(sampDirs.first, false));
					}
					for (auto & rep : sampDirs.second) {
						if (verbose) {
							std::cout << bib::bashCT::bold << "Making Rep Dir: "
									<< bib::bashCT::blue
									<< bib::appendAsNeeded(sampDir, "/") + rep.first
									<< bib::bashCT::reset << std::endl;
						}

						std::string repDir = bib::files::makeDir(sampDir,
								bib::files::MkdirPar(rep.first, false));
						rep.second = bib::files::join(cwd, repDir);
					}
				}
			}
		}
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << bib::bashCT::boldRed(e.what()) << std::endl;
		throw std::runtime_error { ss.str() };
	}
	//log the locations
	std::string indexDir = bib::files::makeDir(mainDirectoryName,
			bib::files::MkdirPar("locationByIndex", false));
	for (const auto & targetDirs : sampleDirWithSubDirs) {
		std::ofstream indexFile;
		openTextFile(indexFile, indexDir + targetDirs.first, ".tab.txt", false,
				false);
		for (const auto & sampDirs : targetDirs.second) {
			for (auto & rep : sampDirs.second) {
				indexFile << rep.first << "\t" << rep.second << std::endl;
			}
		}
	}
}




std::pair<uint32_t, uint32_t> getMaximumRounds(double cnt){
	uint32_t count = 0;
	double currentCnt = cnt;
	while((currentCnt/3.0) > 1){
		currentCnt = currentCnt/3.0;
		++count;
	}

	return {count,  static_cast<uint32_t>(currentCnt) };
}




std::vector<char> getAAs(bool noStopCodon){
	std::vector<char> ans = getVectorOfMapKeys(aminoAcidInfo::infos::allInfo);
	if(noStopCodon){
		removeElement(ans, '*');
	}
	return ans;
}





void makeMultipleSampleDirectory(const std::string& barcodeFilename,
                                 const std::string& mainDirectoryName) {
  seqUtil seqHelper = seqUtil();
  VecStr barcodes = seqHelper.getBarcodesInOrderTheyAppear(barcodeFilename);
  if (barcodes.size() % 2 != 0) {
  	std::stringstream ss;
    ss << "need to have an even amount of barcodes, currently have "
              << barcodes.size() << std::endl;
    throw std::runtime_error{ss.str()};
  }
  auto longestSharedName = seqUtil::findLongestSharedSubString(barcodes);
  VecStr combinedNames;
  std::string directoryName = bib::files::makeDir("./", bib::files::MkdirPar(mainDirectoryName, false));
  for (size_t i = 0; i < barcodes.size(); i += 2) {
    auto firstName = bib::replaceString(barcodes[i], longestSharedName, "");
    auto secondName = bib::replaceString(barcodes[i + 1], longestSharedName, "");
    auto combinedName = longestSharedName + firstName + secondName;
    combinedNames.push_back(combinedName);
    auto combinedDirectoryName = bib::files::makeDir(directoryName, bib::files::MkdirPar(combinedName, false));
    bib::files::makeDir(combinedDirectoryName, bib::files::MkdirPar(barcodes[i], false));
    bib::files::makeDir(combinedDirectoryName, bib::files::MkdirPar(barcodes[i + 1], false));
  }
  // std::cout<<vectorToString(combinedNames)<<std::endl;
}
void makeSampleDirectoriesWithSubDirectories(
    const std::string& barcodeFilename, const std::string& mainDirectoryName) {

  std::map<std::pair<std::string, std::string>, VecStr> directoryNames;
  table inTab(barcodeFilename, "\t");
  for (const auto& fIter : inTab.content_) {

    auto indexAndSampleName = std::make_pair(fIter[0], fIter[1]);
    int count = 0;
    for (const auto& lIter : fIter) {
      if (count == 0) {
        // TODO: this shouldn't be needed, as [] always return valid object
        directoryNames[indexAndSampleName] = {};
      } else if (count == 1) {

      } else {
        directoryNames[indexAndSampleName].emplace_back(lIter);
      }
      ++count;
    }
  }

  VecStr actualDirectoryNames;
  std::map<std::string, std::vector<std::pair<std::string, std::string>>>
      indexMap;
  for (const auto& dIter : directoryNames) {
    actualDirectoryNames.push_back(dIter.first.second);
    if (indexMap.find(dIter.first.first) == indexMap.end()) {
      // TODO: this shouldn't be needed, as [] always return valid object
      indexMap[dIter.first.first] = {};
    }
    for (const auto& sIter : dIter.second) {
      auto currentWorkingDirectory = bib::files::get_cwd();
      indexMap[dIter.first.first]
          .push_back({sIter, currentWorkingDirectory + "/" + mainDirectoryName +
                                 dIter.first.second + "/" + sIter + "/"});
      actualDirectoryNames.push_back(dIter.first.second + "/" + sIter);
    }
  }

  std::cout << "Making following directories" << std::endl;
  std::cout << vectorToString(actualDirectoryNames, "\n");
  std::cout << std::endl;
  std::string indexDir = bib::files::makeDir(mainDirectoryName, bib::files::MkdirPar("locationByIndex", false));
  for (const auto& indexIter : indexMap) {
    std::ofstream indexFile;
    openTextFile(indexFile, indexDir + bib::replaceString(indexIter.first, " ", "_"),
                 ".tab.txt", false, false);
    for (const auto& spIter : indexIter.second) {
      indexFile << spIter.first << "\t"
                << bib::replaceString(spIter.second, "./", "") << std::endl;
    }
  }
  for (const auto& e : actualDirectoryNames) {
  	bib::files::makeDir(mainDirectoryName, bib::files::MkdirPar(e, false));
  }
}

void processKrecName(readObject& read, bool post) {
  if (post) {
    auto toks = tokenizeString(read.seqBase_.name_, "|");
    read.seqBase_.cnt_ = 20000 * (atof(toks[1].c_str()) / 100.00);
  } else {
    auto toks = tokenizeString(read.seqBase_.name_, "_");
    read.seqBase_.cnt_ = atoi(toks[1].c_str());
  }
  read.updateName();
}



/////////tools for finding additional location output
std::string makeIDNameComparable(const std::string& idName) {
  std::string ans = bib::replaceString(idName, "MID0", "M");
  ans = bib::replaceString(ans, "M0", "M");
  return bib::replaceString(ans, "MID", "M");
}
std::string findIdNameFromFileName(const std::string& filename) {
  auto periodPos = filename.rfind(".");
  auto lastMPos = filename.rfind("MID");

  if (periodPos != std::string::npos && !isdigit(filename[periodPos - 1])) {
    auto underScorePos = filename.rfind("_");
    if(std::string::npos == underScorePos || lastMPos > underScorePos){
    	periodPos = filename.rfind(".");
    }else{
    	periodPos = filename.rfind("_");
    }
  }
  //std::cout << "lastMPos: " << lastMPos << std::endl;
  if(lastMPos == std::string::npos){
  	return filename;
  }
  //std::cout << filename.substr(lastMPos, periodPos - lastMPos) << std::endl;
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

VecStr getPossibleDNASubstringsForProtein(const std::string& seq,
                                          const std::string& protein,
                                          const std::string& seqType) {

  VecStr ans;
  std::vector<std::unordered_map<size_t, size_t>> positions;
  std::map<char, VecStr> codonMap;
  if (seqType == "DNA") {
  	for(const auto & aa : aminoAcidInfo::infos::allInfo){
  		codonMap[aa.first] = aa.second.dnaCodons_;
  	}
  } else if (seqType == "RNA") {
  	for(const auto & aa : aminoAcidInfo::infos::allInfo){
  		codonMap[aa.first] = aa.second.rnaCodons_;
  	}
  } else {
  	std::stringstream ss;
    ss << "Unrecognized seqType : " << seqType << std::endl;
    throw std::runtime_error{ss.str()};
  }
  std::vector<size_t> firstCodonPositions;
  for (const auto& codon : codonMap.at(protein.front())) {
    addOtherVec(firstCodonPositions, findOccurences(seq, codon));
  }
  for (const auto& c : protein) {
    if (c == protein.front()) {
      continue;
    }
    std::vector<size_t> currentAminoAcid;
    for (const auto& codon : codonMap.at(c)) {
      addOtherVec(currentAminoAcid, findOccurences(seq, codon));
    }
    std::unordered_map<size_t, size_t> currentPositons;
    for (const auto& pos : currentAminoAcid) {
      currentPositons[pos] = pos;
    }
    positions.push_back(currentPositons);
  }
  std::vector<std::vector<size_t>> possiblePositions;
  for (const auto& firstPos : firstCodonPositions) {
    if (positions.back().find(firstPos + 3 * positions.size()) ==
        positions.back().end()) {
      continue;
    }
    std::vector<size_t> currentPositions;
    currentPositions.push_back(firstPos);
    bool madeItToTheEnd = true;
    size_t count = 1;
    for (const auto& codonPositions : positions) {
      if (codonPositions.find(firstPos + count * 3) != codonPositions.end()) {
        currentPositions.push_back(firstPos + count * 3);
      } else {
        madeItToTheEnd = false;
        break;
      }
      ++count;
    }

    if (madeItToTheEnd) {
      possiblePositions.push_back(currentPositions);
    }
  }
  for (const auto& possible : possiblePositions) {
    ans.emplace_back(getStringFromSubstrings(seq, possible, 3));
  }
  return ans;
}
VecStr findPossibleDNA(const std::string& seq, const std::string& protein,
                       const std::string& seqType, bool checkComplement) {
  VecStr ans = getPossibleDNASubstringsForProtein(seq, protein, seqType);
  if (checkComplement) {
    std::string complementSeq = seqUtil::reverseComplement(seq, "DNA");
    VecStr complementAns =
        getPossibleDNASubstringsForProtein(complementSeq, protein, seqType);
    for (auto& str : complementAns) {
      str = seqUtil::reverseComplement(str, "DNA");
    }
    addOtherVec(ans, complementAns);
  }
  return ans;
}
VecStr getAllCycloProteinFragments(const std::string& protein) {
  VecStr ans;
  for (auto i : iter::range(1, (int)protein.length())) {
    for (auto j : iter::range(protein.length())) {
      std::string currentFragment = "";
      if (j + i > protein.size()) {
        currentFragment = protein.substr(j, protein.size() - j) +
                          protein.substr(0, i - (protein.size() - j));
      } else {
        currentFragment = protein.substr(j, i);
      }
      ans.push_back(currentFragment);
    }
  }
  ans.push_back(protein);
  return ans;
}

VecStr getAllLinearProteinFragments(const std::string& protein) {
  VecStr ans;
  ans.push_back("");
  for (auto i : iter::range(1, (int)protein.length())) {
    for (auto j : iter::range(protein.length())) {
      std::string currentFragment = "";
      if (j + i > protein.size()) {

      } else {
        currentFragment = protein.substr(j, i);
      }
      ans.push_back(currentFragment);
    }
  }
  ans.push_back(protein);
  return ans;
}
std::multimap<int, std::string> getProteinFragmentSpectrum(
    const std::string& protein) {
  std::multimap<int, std::string, std::less<int>> ans;
  ans.insert({0, ""});
  VecStr fragments = getAllCycloProteinFragments(protein);
  for (const auto frag : fragments) {
    ans.insert({seqUtil::calculateWeightOfProteinInt(frag), frag});
  }
  return ans;
}

std::vector<int> getRealPossibleWeights(const std::vector<int>& spectrum) {
  std::vector<int> possibleWeights;
  for (auto i : spectrum) {
    if (i < 187) {
      possibleWeights.push_back(i);
    }
  }
  removeDuplicates(possibleWeights);
  // possibleWeights = getUnique(possibleWeights);
  std::vector<int> realPossibleWeights;
  for (const auto& pw : possibleWeights) {
    if (aminoAcidInfo::infos::weightIntToAminoAcid.find(pw) !=
        aminoAcidInfo::infos::weightIntToAminoAcid.end()) {
      realPossibleWeights.push_back(pw);
    }
  }
  return realPossibleWeights;
}


VecStr organizeLexicallyKmers(const VecStr& input, int colNum) {
  VecStr secondInputAsNumbers;

  MapStrStr stone;
  char count = 'A';
  for (const auto& sIter : input) {
    stone[std::string(1, count)] = sIter;
    secondInputAsNumbers.push_back(std::string(1, count));
    ++count;
  }
  uint64_t reserveSize = 0;
  for (auto i = 1; i <= colNum; ++i) {
    reserveSize += Factorial(i);
  }
  std::vector<VecStr> secondAns;
  secondAns.reserve(reserveSize);

  for (auto i = 1; i <= colNum; ++i) {
    auto temp = permuteVector(secondInputAsNumbers, i);
    secondAns.insert(secondAns.end(), temp.begin(), temp.end());
  }
  VecStr thirdAns;
  for (const auto& iter : secondAns) {
    thirdAns.push_back(vectorToString(iter, ""));
  }
  std::sort(thirdAns.begin(), thirdAns.end());
  for (auto& iter : thirdAns) {
    translateStringWithKey(iter, stone);
  }
  return thirdAns;
}
// doesn't actually work don't use
VecStr organizeLexicallyKmers(const std::string& input, size_t colNum) {
  std::string transformedString = "";

  MapStrStr stone;
  char count = 'A';
  for (const auto& c : input) {
    stone[std::string(1, count)] = c;
    transformedString.push_back(count);
    ++count;
  }
  uint64_t reserveSize = 0;
  for (size_t i = 1; i <= colNum; ++i) {
    reserveSize += Factorial(i);
  }
  VecStr secondAns;
  secondAns.reserve(reserveSize);
  VecStr fragments = getAllCycloProteinFragments(transformedString);
  VecStr trimmedFragments;
  for (const auto& frag : fragments) {
    if (frag.size() <= colNum) {
      trimmedFragments.push_back(frag);
    }
  }
  std::sort(trimmedFragments.begin(), trimmedFragments.end());
  printVector(trimmedFragments);
  for (const auto& frag : trimmedFragments) {
    auto temp = fastPermuteVectorOneLength(frag);
    printVector(temp);
    addOtherVec(secondAns, temp);
  }
  std::sort(secondAns.begin(), secondAns.end());
  for (auto& iter : secondAns) {
    translateStringWithKey(iter, stone);
  }
  return secondAns;
  /*
   exit(1);
   for (auto i = 1; i <= colNum; ++i) {
   auto temp = permuteVector(secondInputAsNumbers, i);
   secondAns.insert(secondAns.end(), temp.begin(), temp.end());
   }
   VecStr thirdAns;
   for (const auto& iter : secondAns) {
   thirdAns.push_back(vectorToString(iter, ""));
   }
   std::sort(thirdAns.begin(), thirdAns.end());

   return thirdAns;*/
}
uint64_t smallestSizePossible(uint64_t weight) {
  return std::floor(weight / 186.00);
}
int64_t getPossibleNumberOfProteins(
    int64_t weight, std::unordered_map<int64_t, int64_t>& cache) {
  if (0 == weight) {
    return 1;
  }
  if (0 > weight) {
    return 0;
  }
  int64_t count = 0;
  for (auto w : aminoAcidInfo::infos::weightIntToAminoAcid) {
    auto q = weight - w.first;
    // std::cout << q << std::endl;
    if (cache.find(q) == cache.end()) {
      auto t = getPossibleNumberOfProteins(q, cache);
      cache[q] = t;
    }
    count += cache[q];
  }
  return count;
}

probabilityProfile randomlyFindBestProfile(const VecStr& dnaStrings,
                                           const std::vector<VecStr>& allKmers,
                                           int numberOfKmers,
                                           bib::randomGenerator& gen) {
  VecStr randomMers;
  // bool needToSeed=true;

  std::vector<int> randomKmersNums =
      gen.unifRandVector(0, numberOfKmers, dnaStrings.size());
  uint32_t pos = 0;
  for (const auto& i : randomKmersNums) {
    randomMers.push_back(allKmers[pos][i]);
    ++pos;
  }
  probabilityProfile bestProfile(randomMers);
  while (true) {
    VecStr currentMers;
    for (const auto& dString : dnaStrings) {
      currentMers.push_back(bestProfile.mostProbableKmers(dString)[0].k_);
    }
    probabilityProfile currentProfile(currentMers);
    if (currentProfile.score_ < bestProfile.score_) {
      bestProfile = currentProfile;
    } else {
      return bestProfile;
    }
  }
  return bestProfile;
}

probabilityProfile randomMotifSearch(const VecStr& dnaStrings, int kLength,
                                     int numberOfRuns, bool gibs, int gibsRuns,
																		 bib::randomGenerator& gen) {
  std::vector<VecStr> allKmers;
  int numOfKmers = 0;
  for (const auto& dString : dnaStrings) {
    allKmers.push_back(kmerCalculator::getAllKmers(dString, kLength));
    numOfKmers = allKmers.back().size();
  }
  probabilityProfile bestProfile =
      randomlyFindBestProfile(dnaStrings, allKmers, numOfKmers, gen);
  for (int i = 0; i < numberOfRuns; ++i) {
    if (gibs) {
    	probabilityProfile currentProfile = randomlyFindBestProfileGibs(dnaStrings, allKmers,
                                                   numOfKmers, gibsRuns, gen);
      if (currentProfile.score_ < bestProfile.score_) {
        bestProfile = currentProfile;
      }
    } else {
    	probabilityProfile currentProfile =
          randomlyFindBestProfile(dnaStrings, allKmers, numOfKmers, gen);
      if (currentProfile.score_ < bestProfile.score_) {
        bestProfile = currentProfile;
      }
    }

  }
  return bestProfile;
}
probabilityProfile randomlyFindBestProfileGibs(
    const VecStr& dnaStrings, const std::vector<VecStr>& allKmers,
    int numberOfKmers, int runTimes, bib::randomGenerator& gen) {
  VecStr randomMers;
  std::vector<int> randomKmersNums =
      gen.unifRandVector(0, numberOfKmers, dnaStrings.size());
  uint32_t pos = 0;
  for (const auto& i : randomKmersNums) {
    randomMers.push_back(allKmers[pos][i]);
    ++pos;
  }
  probabilityProfile bestProfile(randomMers);
  // run random selection j times
  for (int j = 0; j < runTimes; ++j) {
    // pick random kmer to replace and get new profile without that kmer
    size_t randomStringNum =
        gen.unifRand(0, (int)bestProfile.dnaStrings_.size());
    VecStr currentMotifs;
    for (auto i : iter::range(bestProfile.dnaStrings_.size())) {
      if (i != randomStringNum) {
        currentMotifs.push_back(bestProfile.dnaStrings_[i]);
      }
    }
    probabilityProfile currentProfile(currentMotifs);
    // now select a new kmer by preferentially picking a more probable kmer but
    // still randomly
    std::multimap<double, std::string> kmersByProb;
    double cumProb = 0.0;
    for (const auto& ks : allKmers[randomStringNum]) {
      double currentProb = currentProfile.getProbabilityOfKmer(ks);
      kmersByProb.insert({currentProb, ks});
      cumProb += currentProb;
    }
    double randomProb = gen.unifRand(0.0, cumProb);
    double sumOfProbs = 0;
    std::string randomKmer = "";
    for (const auto& kByProb : kmersByProb) {
      sumOfProbs += kByProb.first;
      if (sumOfProbs > randomProb) {
        currentProfile.add(kByProb.second);
        currentProfile.updateScore();
        randomKmer = kByProb.second;
        break;
      }
    }
    // if new randomly replaced kmer creates a better profile repalce
    // best profile but keep order of strings so that they are correctly
    // replaced latter
    if (currentProfile.score_ < bestProfile.score_) {
      bestProfile.dnaStrings_[randomStringNum] = randomKmer;
      bestProfile = probabilityProfile(bestProfile.dnaStrings_);
    }
  }
  return bestProfile;
}
std::vector<std::vector<int>> growNextCycle(
    const std::vector<std::vector<int>>& previousCycle,
    const std::vector<int>& possibleWeights) {
  std::vector<std::vector<int>> nextCycle;
  nextCycle.reserve(previousCycle.size() * possibleWeights.size());
  for (const auto& pw : possibleWeights) {
    for (auto cycle : previousCycle) {
      cycle.push_back(pw);
      nextCycle.push_back(cycle);
    }
  }
  return nextCycle;
}

bool trimCycle(std::vector<std::vector<int>>& nextCycle,
               std::vector<std::vector<int>>& matchesSpectrum,
               const std::map<int, int>& spectrumToWeights,
               const std::vector<int> specVec) {
  int pos = (int)nextCycle.size();
  std::vector<int> positions;
  for (const auto& cycle : iter::reversed(nextCycle)) {
    --pos;
    auto currentFragments = getAllVecLinearFragments(cycle);
    std::map<int, int> fragmentWeightCounts;
    for (const auto& fragment : currentFragments) {
      auto currentSum = vectorSum(fragment);
      ++fragmentWeightCounts[currentSum];
    }
    bool faulted = false;
    for (const auto& weight : fragmentWeightCounts) {
      if (spectrumToWeights.find(weight.first) == spectrumToWeights.end() ||
          spectrumToWeights.at(weight.first) < weight.second) {
        nextCycle.erase(nextCycle.begin() + pos);
        positions.push_back(pos);
        faulted = true;
        break;
      }
    }
    if (faulted) {
      continue;
    }
    std::vector<int> theoSpectrum;
    theoSpectrum.reserve(currentFragments.size());
    auto currentCycloFragments = getAllVecCycloFragments(cycle);
    for (const auto& fragment : currentCycloFragments) {
      auto currentSum = vectorSum(fragment);
      theoSpectrum.push_back(currentSum);
    }
    std::sort(theoSpectrum.begin(), theoSpectrum.end());
    if (theoSpectrum == specVec) {
      matchesSpectrum.push_back(cycle);
      nextCycle.erase(nextCycle.begin() + pos);
    }
  }
  return !nextCycle.empty();
}

std::vector<std::vector<int>> getPossibleProteinsForSpectrum(
    const std::string& spectrum, bool useAllWeights) {
  std::vector<int> specVec = stringToVector<int>(spectrum);
  std::map<int, int> spectrumToWeights;
  for (const auto& i : specVec) {
    ++spectrumToWeights[i];
  }
  std::vector<int> possibleWeights;
  if (useAllWeights) {
    possibleWeights = getVectorOfMapKeys(aminoAcidInfo::infos::weightIntToAminoAcid);
    // for(const auto & weight : )
  } else {
    possibleWeights = getRealPossibleWeights(specVec);
  }

  std::vector<std::vector<int>> initialCycle;
  for (auto pw : possibleWeights) {
    initialCycle.emplace_back(std::vector<int>(1, pw));
  }
  bool keepGrowing = true;
  std::vector<std::vector<int>> matchesSpectrum;
  int count = 0;
  while (keepGrowing) {
    initialCycle = growNextCycle(initialCycle, possibleWeights);
    keepGrowing =
        trimCycle(initialCycle, matchesSpectrum, spectrumToWeights, specVec);
    std::sort(initialCycle.begin(), initialCycle.end());
    ++count;
  }
  return matchesSpectrum;
}

int scoreSpectrumAgreement(const std::map<int, int>& spectrum,
                           const std::map<int, int>& currentSpectrum) {
  int ans = 0;
  for (const auto& weight : currentSpectrum) {
    // std::cout<<"Weight: "<<weight.first<<" count:
    // "<<weight.second<<std::endl;
    if (spectrum.find(weight.first) == spectrum.end()) {

    } else if (spectrum.at(weight.first) <= weight.second) {
      ans += spectrum.at(weight.first);
      // std::cout<<"ansFirstOption: "<<ans<<std::endl;
    } else if (spectrum.at(weight.first) > weight.second) {
      ans += weight.second;
      // std::cout<<"ansSecondOption: "<<ans<<std::endl;
    } else {
      // this shouldn't happen
      std::cout << "ERROR!!" << std::endl;
    }
  }
  return ans;
}

std::multimap<int, std::vector<int>, std::greater<int>> growNextCycleScore(
    const std::multimap<int, std::vector<int>, std::greater<int>>&
        previousCycle,
    const std::vector<int>& possibleWeights,
    const std::map<int, int>& spectrumCounts, int parentMass, bool linear) {
  std::multimap<int, std::vector<int>, std::greater<int>> nextCycle;
  // nextCycle.reserve(previousCycle.size()*possibleWeights.size());
  for (const auto& pw : possibleWeights) {
    for (auto cycle : previousCycle) {
      cycle.second.push_back(pw);
      std::vector<std::vector<int>> currentFragments;
      if (linear) {
        currentFragments = getAllVecLinearFragments(cycle.second);
      } else {
        currentFragments = getAllVecCycloFragments(cycle.second);
      }
      std::map<int, int> fragmentWeightCounts;
      bool tooBig = false;
      // std::cout<<std::endl;
      for (const auto& fragment : currentFragments) {
        // std::cout<<"fragment: ";printVector(fragment);
        auto currentSum = vectorSum(fragment);
        // std::cout<<"Current sum: "<<currentSum<<std::endl;
        if (currentSum > parentMass) {
          tooBig = true;
        }
        ++fragmentWeightCounts[currentSum];
      }
      if (!tooBig) {
        int currentScore =
            scoreSpectrumAgreement(spectrumCounts, fragmentWeightCounts);
        // std::cout<<"currentScore: "<<currentScore<<std::endl;
        nextCycle.insert({currentScore, cycle.second});
      }
      // exit(1);
    }
  }
  return nextCycle;
}

std::multimap<int, std::vector<int>, std::greater<int>> trimCycleScore(
    std::multimap<int, std::vector<int>, std::greater<int>>& nextCycle,
    std::multimap<int, std::vector<int>, std::greater<int>>& matchesSpectrum,
    int parentMass, int leaderBoardNumber, int& currentLeader) {
  std::multimap<int, std::vector<int>, std::greater<int>> ans;
  int pos = 0;
  int lastScore;
  for (const auto& cycle : nextCycle) {
    if (pos == 0) {
      lastScore = cycle.first;
    }
    ++pos;
    if (pos > leaderBoardNumber && cycle.first != lastScore) {

    } else {
      int mass = vectorSum(cycle.second);
      if (mass == parentMass && cycle.first == currentLeader) {
        matchesSpectrum.insert(cycle);
      } else if (mass == parentMass && cycle.first >= currentLeader) {
        matchesSpectrum.clear();
        matchesSpectrum.insert(cycle);
        currentLeader = cycle.first;
      }
      lastScore = cycle.first;
      ans.insert(cycle);
    }
  }
  return ans;
}
std::multimap<int, std::vector<int>, std::greater<int>>
getPossibleProteinsForSpectrum(const std::string& spectrum,
                               int leaderBoardNumber, bool verbose,
                               bool useAllWeights, bool convolution,
                               int convolutionCutOff, bool linear) {
  std::vector<int> specVec = stringToVector<int>(spectrum);
  std::sort(specVec.begin(), specVec.end());
  int parentMass = specVec.back();
  std::map<int, int> spectrumToWeights;
  for (const auto& i : specVec) {
    ++spectrumToWeights[i];
  }
  /*
   for(const auto & specW : spectrumToWeights){
   std::cout<<specW.first<< ":"<<specW.second<<std::endl;
   }*/
  std::vector<int> possibleWeights;
  if (convolution) {
    possibleWeights = topConvolutionWeights(specVec, convolutionCutOff);
    std::cout << "convolutionCutOff: " << convolutionCutOff << std::endl;
    std::cout << "convolutionPossibleWeights: ";
    printVector(possibleWeights);
  } else if (useAllWeights) {
    possibleWeights = std::vector<int>(144);
    std::iota(possibleWeights.begin(), possibleWeights.end(), 57);
  } else {
    possibleWeights = getVectorOfMapKeys(aminoAcidInfo::infos::weightIntToAminoAcid);
  }

  std::multimap<int, std::vector<int>, std::greater<int>> initialCycle;
  for (auto pw : possibleWeights) {
    initialCycle.insert({0, std::vector<int>(1, pw)});
  }
  bool keepGrowing = true;
  std::multimap<int, std::vector<int>, std::greater<int>> matchesSpectrum;
  int currentLeader = 0;
  int count = 0;
  while (keepGrowing) {
    std::cout << "on cycle " << count << std::endl;
    initialCycle = growNextCycleScore(initialCycle, possibleWeights,
                                      spectrumToWeights, parentMass, linear);
    if (verbose) {
      for (const auto& cycle : initialCycle) {
        std::cout << cycle.first << " " << vectorToString(cycle.second)
                  << std::endl;
      }
    }
    initialCycle = trimCycleScore(initialCycle, matchesSpectrum, parentMass,
                                  leaderBoardNumber, currentLeader);
    if (initialCycle.empty()) {
      keepGrowing = false;
    }
    /*if (initialCycle.empty() || !matchesSpectrum.empty()) {
      keepGrowing=false;
    }*/
    ++count;
  }
  return matchesSpectrum;
}
std::map<int, int> getConvolutionWeights(std::vector<int> experimentalSpectrum,
                                         int multiplicityCutOff, int lowerBound,
                                         int upperBound) {
  std::sort(experimentalSpectrum.begin(), experimentalSpectrum.end());
  std::map<int, int, std::greater<int>> countsOfDifferences;
  for (const auto& i : iter::range(experimentalSpectrum.size())) {
    for (const auto& j : iter::range(i + 1, experimentalSpectrum.size())) {
      int currentDifference = experimentalSpectrum[j] - experimentalSpectrum[i];
      if (currentDifference >= lowerBound && currentDifference <= upperBound) {
        ++countsOfDifferences[currentDifference];
      }
    }
  }
  std::map<int, int> ans;
  for (const auto& count : countsOfDifferences) {
    if (count.second >= multiplicityCutOff) {
      ans.insert(count);
    }
  }
  return ans;
}
std::vector<int> convolutionWeights(std::vector<int> experimentalSpectrum,
                                    int multiplicityCutOff, int lowerBound,
                                    int upperBound) {
  std::map<int, int> counts = getConvolutionWeights(
      experimentalSpectrum, multiplicityCutOff, lowerBound, upperBound);
  return getVectorOfMapKeys(counts);
}
std::vector<int> topConvolutionWeights(std::vector<int> experimentalSpectrum,
                                       int mItems, int lowerBound,
                                       int upperBound) {
  std::map<int, int> allCounts =
      getConvolutionWeights(experimentalSpectrum, 1, lowerBound, upperBound);
  std::multimap<int, int, std::greater<int>> byCounts;
  for (const auto& count : allCounts) {
    byCounts.insert({count.second, count.first});
  }
  int counting = 0;
  int lastCount = byCounts.begin()->first;
  // int testCount=0;
  /*for(const auto & count : byCounts){
    testCount++;
    std::cout<<"testCount: "<<testCount<<" currentCount: "<<count.first<<"
  weight: "<<count.second<<std::endl;
  }*/
  for (const auto& count : byCounts) {
    ++counting;
    if (counting > mItems && count.first != lastCount) {
      break;
    }
    lastCount = count.first;
  }
  // std::cout<<"lastcount: "<<lastCount<<std::endl;
  return convolutionWeights(experimentalSpectrum, lastCount, lowerBound,
                            upperBound);
}
int64_t getMinCoins(int64_t change, const std::vector<int64_t>& coins,
                    std::unordered_map<int64_t, int64_t>& cache) {
  if (0 == change) {
    return 0;
  }
  int64_t count = change;
  for (auto coin : coins) {
    if (change - coin < 0) {

    } else {
      int64_t q = change - coin;
      if (cache.find(q) == cache.end()) {
        auto t = getMinCoins(q, coins, cache) + 1;
        cache[q] = t;
      }
      if (cache[q] < count) {
        count = cache[q];
      }
    }
  }
  return count;
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
    double currentMedian = vectorMedianCopy(getWindowQuals(qual, qWindowSize, i));
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
double getChangeInHydro(const char& firstAA, const char& secondAA) {
  return std::abs(aminoAcidInfo::infos::allInfo.at(firstAA).acidHydrophobicity_ -
                  aminoAcidInfo::infos::allInfo.at(secondAA).acidHydrophobicity_);
}

/*
std::vector<double> getHydroChanges(const std::string& originalCodon,
                                    const VecStr& mutantCodons,
                                    const std::map<std::string, char>& code) {
  std::vector<double> ans;
  if (originalCodon.size() != 3) {
    std::cout << "codon needs to be 3 bases" << std::endl;
    std::cout << originalCodon << std::endl;
    std::cout << originalCodon.size() << std::endl;
    return ans;
  }
  char originalAA = code.at(originalCodon);
  if (originalAA == '*') {
    return ans;
  }
  for (const auto& codon : mutantCodons) {
    char mutantAA = code.at(codon);
    if (mutantAA == '*') {
      continue;
    }
    ans.emplace_back(getChangeInHydro(originalAA, mutantAA));
  }
  return ans;
}*/


VecStr createDegenStrs(const std::string & str){
	VecStr ans;

	char first = *str.begin();
	std::cout << first << std::endl;
	std::cout << std::endl;
	if(degenBaseExapnd.find(first) != degenBaseExapnd.end()){
		for(const auto & nextChar : degenBaseExapnd.at(first)){
			ans.emplace_back(std::string(1,nextChar));
		}
	}else{
		ans.emplace_back(std::string(1,first));
	}

	printVector(ans);
	uint32_t charCount = 0;
	for(const auto & c : str){
		++charCount;
		if(charCount == 1){
			continue;
		}
		if(degenBaseExapnd.find(c) != degenBaseExapnd.end()){
			VecStr adding;
			for(auto & current : ans){
				uint32_t count = 0;
				std::string copy = current;
				for(const auto & nextChar : degenBaseExapnd.at(c)){
					if(count == 0){
						current.push_back(nextChar);
					}else{
						adding.emplace_back(copy);
						adding.back().push_back(nextChar);
					}
					++count;
				}
			}
			addOtherVec(ans,adding);
		}else{
			for(auto & current : ans){
				current.push_back(c);
			}
		}
	}
	return ans;
}




table getSeqPosTab(const std::string & str){
	table out(VecStr{"char", "pos"});
	for(const auto & pos : iter::range(str.size())){
		out.content_.emplace_back(toVecStr(str[pos], pos));
	}
	return out;
}



}  // namespace bib
