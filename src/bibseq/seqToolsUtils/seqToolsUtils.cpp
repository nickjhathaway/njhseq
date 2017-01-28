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
		std::string dirName = bfs::basename(dir);
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
				bfs::path targetDir = bib::files::makeDir(topDir,
						bib::files::MkdirPar(targetDirs.first, false));
				for (auto & sampDirs : targetDirs.second) {
					if (verbose) {
						std::cout << bib::bashCT::bold << "Making Samp Dir: "
								<< bib::bashCT::green << bib::files::make_path(targetDir,  sampDirs.first)
								<< bib::bashCT::reset << std::endl;
					}

					bfs::path sampDir = bib::files::makeDir(targetDir,
							bib::files::MkdirPar(sampDirs.first, false));
					for (auto & rep : sampDirs.second) {
						if (verbose) {
							std::cout << bib::bashCT::bold << "Making Rep Dir: "
									<< bib::bashCT::blue <<  bib::files::make_path(sampDir, rep.first)
									<< bib::bashCT::reset << std::endl;
						}

						bfs::path repDir = bib::files::makeDir(sampDir,
								bib::files::MkdirPar(rep.first, false));
						rep.second = bib::files::make_path(cwd, repDir).string();
					}
				}
			}
		} else {
			for (auto & targetDirs : sampleDirWithSubDirs) {
				for (auto & sampDirs : targetDirs.second) {
					std::string sampDir = bib::files::join(topDir, sampDirs.first).string();
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
								bib::files::MkdirPar(rep.first, false)).string();
						rep.second = bib::files::join(cwd, repDir).string();
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
	bfs::path indexDir = bib::files::makeDir(mainDirectoryName,
			bib::files::MkdirPar("locationByIndex", false));
	for (const auto & targetDirs : sampleDirWithSubDirs) {
		std::ofstream indexFile;
		openTextFile(indexFile, bib::files::make_path(indexDir, targetDirs.first).string(), ".tab.txt", false,
				false);
		for (const auto & sampDirs : targetDirs.second) {
			for (auto & rep : sampDirs.second) {
				indexFile << rep.first << "\t" << rep.second << std::endl;
			}
		}
	}
}









std::vector<char> getAAs(bool noStopCodon){
	std::vector<char> ans = getVectorOfMapKeys(aminoAcidInfo::infos::allInfo);
	if(noStopCodon){
		removeElement(ans, '*');
	}
	return ans;
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
