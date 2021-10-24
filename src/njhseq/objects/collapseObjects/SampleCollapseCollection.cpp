/*
 * SampleCollapseCollection.cpp
 *
 *  Created on: Jul 24, 2016
 *      Author: nick
 */

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


#include "SampleCollapseCollection.hpp"
#include "njhseq/objects/seqObjects/Clusters/clusterUtils.hpp"


namespace njhseq {

namespace collapse {


Json::Value SampleCollapseCollection::PreFilteringCutOffs::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["clusterSizeCutOff"] = njh::json::toJson(clusterSizeCutOff );
	ret["sampleMinReadCount"] = njh::json::toJson(sampleMinReadCount );
	ret["replicateMinReadCount"] = njh::json::toJson(replicateMinReadCount );

	return ret;
}


std::unordered_map<std::string, double> SampleCollapseCollection::processCustomCutOffs(
		const bfs::path &customCutOffsFnp, const VecStr &allSamples,
		double defaultFracCutOff) {
	std::unordered_map<std::string, double> ret;
	if ("" != customCutOffsFnp.string()) {
		table customCutOffsTab(customCutOffsFnp.string(), "\t", true);
		customCutOffsTab.checkForColumnsThrow(VecStr { "sample", "cutOff" },
				__PRETTY_FUNCTION__);
		for (const auto rowPos : iter::range(customCutOffsTab.content_.size())) {
			ret[customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos("sample")]] =
					njh::lexical_cast<double>(
							customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos(
									"cutOff")]);
		}
	}
	for (const auto &samp : allSamples) {
		if (!njh::in(samp, ret)) {
			ret[samp] = defaultFracCutOff;
		}
	}
	return ret;
}



SampleCollapseCollection::PreFilteringCutOffs::PreFilteringCutOffs() {
}
SampleCollapseCollection::PreFilteringCutOffs::PreFilteringCutOffs(
		const Json::Value & coreJson) {
	njh::json::MemberChecker checker(coreJson);

	checker.failMemberCheckThrow(VecStr {"clusterSizeCutOff", "sampleMinReadCount",
			"replicateMinReadCount"}, __PRETTY_FUNCTION__);
	clusterSizeCutOff = coreJson["clusterSizeCutOff"].asUInt();
	sampleMinReadCount = coreJson["sampleMinReadCount"].asUInt();
	replicateMinReadCount = coreJson["replicateMinReadCount"].asUInt();

}

SampleCollapseCollection::SampleCollapseCollection(SeqIOOptions inputOptions,
		const bfs::path & inputDir,
		const bfs::path & outputDir,
		const PopNamesInfo & popNames,
		PreFilteringCutOffs preFiltCutOffs) :
		inputOptions_(inputOptions),
		masterInputDir_(inputDir),
		masterOutputDir_(outputDir),
		popNames_(popNames),
		preFiltCutOffs_(preFiltCutOffs){
	samplesOutputDir_ = njh::files::make_path(masterOutputDir_, "samplesOutput");
	populationOutputDir_ = njh::files::make_path(masterOutputDir_, "population");
	njh::files::makeDir(njh::files::MkdirPar(samplesOutputDir_.string(), true));

}

SampleCollapseCollection::SampleCollapseCollection(const Json::Value & coreJson){

	njh::json::MemberChecker checker(coreJson);

	checker.failMemberCheckThrow(VecStr {"preFiltCutOffs_", "inputOptions_",
			"masterInputDir_", "masterOutputDir_", "popNames_",
			"populationOutputDir_", "samplesOutputDir_", "lowRepCntSamples_" }, __PRETTY_FUNCTION__);

	inputOptions_ = SeqIOOptions(coreJson["inputOptions_"].toStyledString());
	masterInputDir_ = coreJson["masterInputDir_"].asString();
	masterOutputDir_ = coreJson["masterOutputDir_"].asString();
	std::set<std::string> controlSamples;//this is for an update so older directories don't error out
	if(coreJson["popNames_"].isMember("controlSamples_")){
		controlSamples = njh::json::jsonArrayToSet<std::string>(
				coreJson["popNames_"]["controlSamples_"], [](const Json::Value &val) {
					return val.asString();
				});
	}
	popNames_ = PopNamesInfo(coreJson["popNames_"]["populationName_"].asString(),
			njh::json::jsonArrayToSet<std::string>(coreJson["popNames_"]["samples_"],
					[](const Json::Value &val) {
						return val.asString();
					})
			,
			controlSamples);
	preFiltCutOffs_ = PreFilteringCutOffs(coreJson["preFiltCutOffs_"]);

	populationOutputDir_ = coreJson["populationOutputDir_"].asString();
	samplesOutputDir_ = coreJson["samplesOutputDir_"].asString();
	passingSamples_ = njh::json::jsonArrayToVec<std::string>(coreJson["passingSamples_"],
			[](const Json::Value & val){ return val.asString();});
	lowRepCntSamples_ = njh::json::jsonArrayToVec<std::string>(coreJson["lowRepCntSamples_"],
				[](const Json::Value & val){ return val.asString();});

	//load in group meta data if it is available
	if(bfs::exists(njh::files::make_path(masterOutputDir_, "groups", "groupMetaData.json"))){
		groupMetaData_ = std::make_unique<MultipleGroupMetaData>(
				MultipleGroupMetaData::fromJson(
						njh::json::parseFile(
								njh::files::make_path(masterOutputDir_, "groups",
										"groupMetaData.json").string())));

		groupDataPaths_ = std::make_unique<AllGroupDataPaths>(
				njh::files::make_path(njh::files::normalize(masterOutputDir_), "groups"),
				groupMetaData_);
	}
}


void SampleCollapseCollection::addGroupMetaData(
		const bfs::path & groupingsFile) {

	groupMetaData_ = std::make_unique<MultipleGroupMetaData>(njh::files::normalize(groupingsFile), popNames_.samples_);
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << njh::conToStr(popNames_.samples_) << std::endl;
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << njh::conToStr(groupMetaData_->samples_) << std::endl;
	//std::cout << njh::conToStr(groupMetaData_->missingSamples_) << std::endl;


	/**@todo add some output to output directory that states the samples missing and missing meta*/
	groupDataPaths_ = std::make_unique<AllGroupDataPaths>(
			njh::files::make_path(njh::files::normalize(masterOutputDir_), "groups"),
			groupMetaData_);
}

uint32_t SampleCollapseCollection::numOfSamples() const {
	return popNames_.samples_.size();
}

bool SampleCollapseCollection::hasSample(const std::string & sampleName) const {
	return popNames_.hasSample(sampleName);
}

void SampleCollapseCollection::checkForSampleThrow(const std::string & funcName,
		const std::string & sampleName) const {
	if (!hasSample(sampleName)) {
		std::stringstream ss;
		ss << funcName << ": error, sample " << sampleName << " isn't in samples_" << "\n";
		ss << "Options are: " << njh::conToStr(popNames_.samples_, ",") << "\n";
		throw std::runtime_error { ss.str() };
	}
}

void SampleCollapseCollection::setUpSample(const std::string & sampleName,
		aligner & alignerObj, const collapser & collapserObj,
		const ChimeraOpts & chiOpts) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	auto sampleDir = njh::files::make_path(masterInputDir_, sampleName);
	auto analysisFiles = njh::files::listAllFiles(sampleDir.string(), true, {
			std::regex { "^" + inputOptions_.firstName_.string() + "$" } }, 2);
	std::vector<RepFile> repFiles;
	for (const auto & af : analysisFiles) {
		auto fileToks = njh::tokenizeString(
				bfs::relative(af.first, sampleDir).string(), "/");
		if (2 != fileToks.size()) {
			std::stringstream ss;
			ss << sampleDir << std::endl;
			ss << "File path should be two levels deep, not " << fileToks.size()
					<< " for " << bfs::relative(af.first, sampleDir).string()
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		repFiles.emplace_back(fileToks[0], af.first);
	}
	return setUpSample(sampleName, repFiles, alignerObj, collapserObj, chiOpts);
}


void SampleCollapseCollection::setUpSample(const std::string & sampleName,
		const std::unordered_map<std::string, RepSeqs> & seqsByRep,
		aligner & alignerObj,
		const collapser & collapserObj, const ChimeraOpts & chiOpts){
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);

	std::vector<std::vector<cluster>> inputClusters;
	std::vector<std::vector<cluster>> lowCntInputClusters;

	for (const auto & rep : seqsByRep) {

		std::vector<cluster> clusters = baseCluster::convertVectorToClusterVector<cluster>(rep.second.repSeqs_);
		readVecSorter::sortReadVector(clusters, "totalCount");
		// consider adding the sample name in the name as well
		if(rep.second.reNameInput_){
			renameReadNamesNewClusters(clusters, rep.second.repName_, true, true, false);
		}
		if (chiOpts.checkChimeras_) {
			for(auto & clus : clusters){
				clus.seqBase_.unmarkAsChimeric();
			}
			collapserObj.markChimeras(clusters, alignerObj, chiOpts);
		}

		clusterVec::allSetFractionClusters(clusters);
		if(rep.second.reNameInput_){
			readVec::allUpdateName(clusters);
		}
		double totalCount = readVec::getTotalReadCount(clusters);
		if (preFiltCutOffs_.replicateMinReadCount > 0
				&& totalCount < preFiltCutOffs_.replicateMinReadCount) {
			lowCntInputClusters.emplace_back(clusters);
		} else {
			inputClusters.emplace_back(clusters);
		}
	}
	if(!inputClusters.empty()){
		std::lock_guard<std::mutex> lock(mut_);
		if (njh::in(sampleName, sampleCollapses_)) {
			sampleCollapses_[sampleName] = std::make_shared<sampleCollapse>(
					inputClusters, sampleName, preFiltCutOffs_.clusterSizeCutOff);
		} else {
			sampleCollapses_.emplace(sampleName,
					std::make_shared<sampleCollapse>(inputClusters, sampleName,
							preFiltCutOffs_.clusterSizeCutOff));
		}
	} else {
		std::lock_guard<std::mutex> lock(mut_);
		if (njh::in(sampleName, sampleCollapses_)) {
			sampleCollapses_[sampleName] = std::make_shared<sampleCollapse>(
					lowCntInputClusters, sampleName, preFiltCutOffs_.clusterSizeCutOff);
		} else {
			sampleCollapses_.emplace(sampleName,
					std::make_shared<sampleCollapse>(lowCntInputClusters, sampleName,
							preFiltCutOffs_.clusterSizeCutOff));
		}
		lowRepCntSamples_.emplace_back(sampleName);
	}
}

void SampleCollapseCollection::setUpSample(const std::string & sampleName,
		const std::vector<RepFile> & analysisFiles,
		aligner & alignerObj,
		const collapser & collapserObj,
		const ChimeraOpts & chiOpts){
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	std::unordered_map<std::string, RepSeqs> allRepSeqs;
	for (const auto & repf : analysisFiles) {
		auto sampleOpts = inputOptions_;
		sampleOpts.firstName_ = repf.repFnp_.string();
		SeqInput reader(sampleOpts);

		if(njh::in(repf.repName_, allRepSeqs)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have replicate name: " << repf.repName_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		allRepSeqs.emplace(repf.repName_,RepSeqs(repf.repName_, reader.readAllReads<seqInfo>()));
	}
	return setUpSample(sampleName, allRepSeqs, alignerObj, collapserObj, chiOpts);
}

bfs::path SampleCollapseCollection::getSampleFinalHapsPath(
		const std::string & sampleName) const {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	return njh::files::make_path(samplesOutputDir_, sampleName, "final",
			sampleName + inputOptions_.getOutExtension());
}

bfs::path SampleCollapseCollection::getPopFinalHapsPath() const {
	return njh::files::make_path(populationOutputDir_,
			"PopSeqs" + inputOptions_.getOutExtension());

}


bfs::path SampleCollapseCollection::getPopInfoPath() const {
	return njh::files::make_path(populationOutputDir_,
			"populationCluster.tab.txt.gz");
}


bfs::path SampleCollapseCollection::getFinalSampHapsPath(const std::string & sample) const {
	checkForSampleThrow(sample, __PRETTY_FUNCTION__);
	return njh::files::make_path(masterOutputDir_, "final", sample + inputOptions_.getOutExtension());
}

bfs::path SampleCollapseCollection::getSampInfoPath() const {
	/**@todo need to make this standard, currently set by SeekDeep processClusters*/
	return njh::files::make_path(masterOutputDir_, "selectedClustersInfo.tab.txt.gz");
}

bfs::path SampleCollapseCollection::getHapIdTabPath() const {
	return njh::files::make_path(masterOutputDir_,
			"hapIdTable.tab.txt.gz");
}



void SampleCollapseCollection::setUpSampleFromPrevious(
		const std::string & sampleName) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	auto sampPtr = std::make_shared<sampleCollapse>(sampleName);
	auto & samp = *sampPtr;
	auto sampleDir = njh::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName);

	std::map<std::string, double> allSampCounts;
	std::map<std::string, sampInfo> allSampInfos;

	std::map<std::string, double> finalSampCounts;
	std::map<std::string, sampInfo> finalSampInfos;

//	auto finalClustersFiles = njh::files::filesInFolder(
//			njh::files::make_path(sampleDir, "final/clusters"));
//	auto excludedClustersFiles = njh::files::filesInFolder(
//			njh::files::make_path(sampleDir, "excluded/clusters"));
	seqInfo seq;
	auto finalClustersFileOpts = SeqIOOptions(njh::files::make_path(sampleDir, "final/clusters/inputClusters" + inputOptions_.getOutExtension()),
			inputOptions_.inFormat_, true);
	std::unordered_map<std::string, std::vector<seqInfo>> finalClusterInputClusters;
	if(finalClustersFileOpts.inExists()){
		SeqInput reader(finalClustersFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			allSampCounts[seq.getOwnSampName()] += seq.cnt_;
			finalSampCounts[seq.getOwnSampName()] += seq.cnt_;
			MetaDataInName seqMeta(seq.name_);
			auto topClusterName = seqMeta.getMeta("TopClusterName");
			seqMeta.removeMeta("TopClusterName");
			if(seqMeta.meta_.size() == 0){
				MetaDataInName::removeMetaDataInName(seq.name_);
			}else{
				seqMeta.resetMetaInName(seq.name_);
			}
			finalClusterInputClusters[topClusterName].emplace_back(seq);
		}
	}
	MapStrStr refCompInfosCollapsed;
	MapStrStr refCompInfosExcluded;
	if (bfs::exists(
			njh::files::join(sampleDir.string(), "refCompInfos.tab.txt"))) {
		table refCompTab(
				njh::files::join(sampleDir.string(), "refCompInfos.tab.txt"), "\t",
				true);
		//std::cout << __PRETTY_FUNCTION__ << " " << sampleName << std::endl;
		for (const auto & row : refCompTab.content_) {
			//std::cout << "\t" << row[refCompTab.getColPos("read")] << " :  " <<  row[refCompTab.getColPos("bestExpected")] << std::endl;
			refCompInfosCollapsed[row[refCompTab.getColPos("read")]] =
					row[refCompTab.getColPos("bestExpected")];
		}
	}

	if (bfs::exists(
			njh::files::join(sampleDir.string(), "refCompInfosExcluded.tab.txt"))) {
		table refCompTab(
				njh::files::join(sampleDir.string(), "refCompInfosExcluded.tab.txt"), "\t",
				true);
		//std::cout << __PRETTY_FUNCTION__ << " " << sampleName << std::endl;
		for (const auto & row : refCompTab.content_) {
			//std::cout << "\t" << row[refCompTab.getColPos("read")] << " :  " <<  row[refCompTab.getColPos("bestExpected")] << std::endl;
			refCompInfosExcluded[row[refCompTab.getColPos("read")]] =
					row[refCompTab.getColPos("bestExpected")];
		}
	}
	std::unordered_map<std::string, std::vector<seqInfo>> excludedInputClusters;
	auto excludedClustersFileOpts = SeqIOOptions(
			njh::files::make_path(sampleDir,
					"excluded/clusters/inputClusters" + inputOptions_.getOutExtension()),
			inputOptions_.inFormat_, true);
	if (excludedClustersFileOpts.inExists()) {
		SeqInput reader(excludedClustersFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			allSampCounts[seq.getOwnSampName()] += seq.cnt_;
			MetaDataInName seqMeta(seq.name_);
			auto topClusterName = seqMeta.getMeta("TopClusterName");
			seqMeta.removeMeta("TopClusterName");
			if(seqMeta.meta_.size() == 0){
//				std::cout << seq.name_ << std::endl;
				MetaDataInName::removeMetaDataInName(seq.name_);
//				std::cout << seq.name_ << std::endl;
//				exit(1);

			}else{
				seqMeta.resetMetaInName(seq.name_);
			}
			excludedInputClusters[topClusterName].emplace_back(seq);
		}
	}



	for (const auto & sCount : allSampCounts) {
		allSampInfos[sCount.first] = sampInfo(sCount.first, sCount.second);
	}

	for (const auto & sCount : finalSampCounts) {
		finalSampInfos[sCount.first] = sampInfo(sCount.first, sCount.second);
	}
	auto finalFileOpts = SeqIOOptions(
			njh::files::make_path(sampleDir, "final/",
					sampleName + inputOptions_.getOutExtension()).string(),
			inputOptions_.inFormat_, true);

	auto excludedFileOpts = SeqIOOptions(
			njh::files::make_path(sampleDir, "excluded/",
					sampleName + inputOptions_.getOutExtension()).string(),
			inputOptions_.inFormat_, true);

	seqInfo subInfo;
	if(finalFileOpts.inExists()){
		SeqInput reader(finalFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
//			auto finalSubFileOpts = SeqIOOptions(
//					njh::files::make_path(sampleDir, "final/clusters",
//							seq.name_ + inputOptions_.getOutExtension()).string(),
//					inputOptions_.inFormat_, true);
//			SeqInput subReader(finalSubFileOpts);
//			auto subReads = subReader.readAllReads<seqInfo>();
			samp.collapsed_.clusters_.emplace_back(
					sampleCluster(seq, finalClusterInputClusters.at(seq.name_), finalSampInfos));
			samp.collapsed_.clusters_.back().expectsString =
					refCompInfosCollapsed[samp.collapsed_.clusters_.back().seqBase_.name_];
			for (const auto & subRead : finalClusterInputClusters.at(seq.name_)) {
				samp.input_.clusters_.emplace_back(
						sampleCluster(subRead, allSampInfos));
			}
		}
	}
	if(excludedFileOpts.inExists()){
		SeqInput reader(excludedFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
//			auto excludedSubFileOpts = SeqIOOptions(
//					njh::files::make_path(sampleDir, "excluded/clusters",
//							seq.name_ + inputOptions_.getOutExtension()).string(),
//					inputOptions_.inFormat_, true);
//			SeqInput subReader(excludedSubFileOpts);
//			auto subReads = subReader.readAllReads<seqInfo>();
			samp.excluded_.clusters_.emplace_back(
					sampleCluster(seq, excludedInputClusters.at(seq.name_), allSampInfos));
			samp.excluded_.clusters_.back().expectsString =
					refCompInfosExcluded[samp.excluded_.clusters_.back().seqBase_.name_];
			for (const auto & subRead : excludedInputClusters.at(seq.name_)) {
				samp.input_.clusters_.emplace_back(
						sampleCluster(subRead, allSampInfos));
			}
		}
	}
	samp.updateExclusionInfos();
	samp.updateInitialInfos();
	samp.updateCollapsedInfos();

	sampleCollapses_[sampleName] = sampPtr;


}

void SampleCollapseCollection::clusterSample(const std::string & sampleName,
		aligner & alignerObj, const collapser & collapserObj,
		const CollapseIterations & colIters) {
	std::string sortBy = "fraction";
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);

	auto samp = sampleCollapses_.at(sampleName);
	samp->cluster(collapserObj, colIters, sortBy, alignerObj);

	samp->updateCollapsedInfos();
	samp->updateExclusionInfos();
	samp->renameClusters(sortBy);
}

void SampleCollapseCollection::collapseLowFreqOneOffsSample(
		const std::string & sampleName, aligner & alignerObj,
		const collapser & collapserObj, double lowFreqMultiplier) {

	std::string sortBy = "fraction";
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);

	auto samp = sampleCollapses_.at(sampleName);
	samp->collapseLowFreqOneOffs(lowFreqMultiplier, alignerObj, collapserObj);

	samp->updateCollapsedInfos();
	samp->updateExclusionInfos();
	samp->renameClusters(sortBy);
}

void SampleCollapseCollection::dumpSample(const std::string & sampleName) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	auto samp = sampleCollapses_.at(sampleName);
	auto sampDir = njh::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName);
	auto finalDir = njh::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName, "final");
	auto finalClustersDir = njh::files::make_path(masterOutputDir_,
			"samplesOutput", sampleName, "final", "clusters");
	auto excludedDir = njh::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName, "excluded");
	auto excludedClusteredDir = njh::files::make_path(masterOutputDir_,
			"samplesOutput", sampleName, "excluded", "clusters");

	njh::files::makeDir(njh::files::MkdirPar(sampDir.string(), true));
	njh::files::makeDirP(njh::files::MkdirPar(finalClustersDir.string()));
	njh::files::makeDirP(njh::files::MkdirPar(excludedClusteredDir.string()));
	//final
	if(!samp->collapsed_.clusters_.empty()){
		SeqIOOptions finalOutOpts(
				njh::files::join(finalDir.string(), sampleName).string()
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		SeqOutput::write(samp->collapsed_.clusters_, finalOutOpts);
		SeqIOOptions clusOutOpts(
				njh::files::make_path(finalClustersDir.string(), "inputClusters").string()
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		clusOutOpts.out_.overWriteFile_ = true;
		SeqOutput inputClustersWriter(clusOutOpts);
		inputClustersWriter.openOut();
		for (const auto & clus : samp->collapsed_.clusters_) {
			for(auto & subClus : clus.reads_){
				MetaDataInName subClusMeta;
				bool hadMeta = false;
				if(MetaDataInName::nameHasMetaData(subClus->seqBase_.name_)){
					subClusMeta.processNameForMeta(subClus->seqBase_.name_, true);
					hadMeta = true;
				}
				auto writeCluster = subClus->seqBase_;
				subClusMeta.addMeta("TopClusterName", clus.seqBase_.name_);
				if(hadMeta){
					subClusMeta.resetMetaInName(writeCluster.name_);
				}else{
					subClusMeta.resetMetaInName(writeCluster.name_, writeCluster.name_.rfind("_t"));
				}
				inputClustersWriter.write(writeCluster);
			}
		}
	}

	//excluded
	if(!samp->excluded_.clusters_.empty()){
		SeqIOOptions excluddOutOpts(
				njh::files::join(excludedDir.string(), sampleName).string()
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		SeqOutput::write(samp->excluded_.clusters_, excluddOutOpts);

		SeqIOOptions clusOutOpts(
				njh::files::make_path(excludedClusteredDir.string(), "inputClusters").string()
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		clusOutOpts.out_.overWriteFile_ = true;
		SeqOutput inputClustersWriter(clusOutOpts);
		inputClustersWriter.openOut();
		for (const auto & clus : samp->excluded_.clusters_) {
			for(auto & subClus : clus.reads_){
				MetaDataInName subClusMeta;
				bool hadMeta = false;
				if(MetaDataInName::nameHasMetaData(subClus->seqBase_.name_)){
					subClusMeta.processNameForMeta(subClus->seqBase_.name_, true);
					hadMeta = true;
				}
				auto writeCluster = subClus->seqBase_;
				subClusMeta.addMeta("TopClusterName", clus.seqBase_.name_);
				if(hadMeta){
					subClusMeta.resetMetaInName(writeCluster.name_);
				}else{
					subClusMeta.resetMetaInName(writeCluster.name_, writeCluster.name_.rfind("_t"));
				}
				inputClustersWriter.write(writeCluster);
			}
		}
	}
	if (!samp->collapsed_.clusters_.empty()
			&& "" != samp->collapsed_.clusters_.front().expectsString) {
		std::map<std::string, std::string> refCompInfos;
		for (const auto & clus : samp->collapsed_.clusters_) {
			refCompInfos[clus.seqBase_.name_] = clus.expectsString;
		}


		table refCompTab(refCompInfos, VecStr { "read", "bestExpected" });
		refCompTab.outPutContents(
				TableIOOpts(
						OutOptions(
								njh::files::join(sampDir.string(), "refCompInfos.tab.txt")),
						"\t", true));

		std::map<std::string, std::string> refCompInfosExcluded;
		for (const auto & clus : samp->excluded_.clusters_) {
			refCompInfosExcluded[clus.seqBase_.name_] = clus.expectsString;
		}
		table refCompExcludedTab(refCompInfosExcluded, VecStr { "read", "bestExpected" });
		refCompExcludedTab.outPutContents(
				TableIOOpts(
						OutOptions(
								njh::files::join(sampDir.string(), "refCompInfosExcluded.tab.txt")),
						"\t", true));
	}
	if(!keepSampleInfoInMemory_){
		clearSample(sampleName);
	}
}

void SampleCollapseCollection::clearSample(const std::string & sampleName) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	sampleCollapses_.at(sampleName) = nullptr;
}

void SampleCollapseCollection::clearSamples(const VecStr & sampleNames) {
	for (const auto & sampleName : sampleNames) {
		clearSample(sampleName);
	}
}


void SampleCollapseCollection::setPassingSamples(){
	passingSamples_.clear();
	for (const auto & sampleName : popNames_.samples_) {
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampleName);
		}
		if(!njh::in(sampleName, lowRepCntSamples_) && sampleCollapses_[sampleName]->collapsed_.info_.totalReadCount_ >= preFiltCutOffs_.sampleMinReadCount){
			passingSamples_.emplace_back(sampleName);
		}
		if(!keepSampleInfoInMemory_){
			clearSample(sampleName);
		}
	}
}

std::vector<sampleCluster> SampleCollapseCollection::createPopInput() {
	passingSamples_.clear();
	oututSampClusToOldNameKey_.clear();
	std::vector<sampleCluster> output;
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for (const auto & sampleName : popNames_.samples_) {
		if(!keepSampleInfoInMemory_){
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			setUpSampleFromPrevious(sampleName);
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "sampleName: " <<  sampleName << std::endl;
//		std::cout << "njh::in(sampleName, lowRepCntSamples_): " << njh::colorBool(njh::in(sampleName, lowRepCntSamples_))<< std::endl;
//		std::cout << "njh::in(sampleName, sampleCollapses_): " << njh::colorBool(njh::in(sampleName, sampleCollapses_))<< std::endl;

		if (njh::in(sampleName, lowRepCntSamples_)
				|| sampleCollapses_[sampleName]->collapsed_.info_.totalReadCount_
						< preFiltCutOffs_.sampleMinReadCount) {
			continue;
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		passingSamples_.emplace_back(sampleName);
		double totalReadCnt = 0;
		for (const auto &out : sampleCollapses_[sampleName]->collapsed_.clusters_) {
			totalReadCnt += out.seqBase_.cnt_;
		}
		std::map<std::string, sampInfo> outSampInfos { { sampleName, sampInfo(sampleName, totalReadCnt) } };
		for (const auto &out : sampleCollapses_[sampleName]->collapsed_.clusters_) {
			output.emplace_back(sampleCluster(out.seqBase_, outSampInfos));
			//output.emplace_back(sampleCluster(out.createRead(), outSampInfos));
			std::string oldName = output.back().seqBase_.name_;
			output.back().updateName();
			oututSampClusToOldNameKey_[sampleName][oldName] = output.back().seqBase_.name_;
			output.back().reads_.front()->seqBase_.name_ = output.back().seqBase_.name_;
		}
		if(!keepSampleInfoInMemory_){
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			clearSample(sampleName);
		}
	}



//	std::cout << "size of sample Cluster" << sizeof(sampleCluster) << std::endl;

//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return output;
}

void SampleCollapseCollection::doPopulationClustering(
		const std::vector<sampleCluster> & input, aligner & alignerObj,
		const collapser & collapserObj, const CollapseIterations & popColIters) {

	popCollapse_ = std::make_unique<collapse::populationCollapse>(input, popNames_.populationName_);

	popCollapse_->popCluster(collapserObj, popColIters, "fraction", alignerObj);

	std::unordered_map<std::string, uint32_t> clusterTotals;
	for (const auto & clus : popCollapse_->collapsed_.clusters_) {
		for (const auto & subClus : clus.reads_) {
			++clusterTotals[subClus->getOwnSampName()];
		}
	}
	popCollapse_->collapsed_.info_.cois_.clear();
	for(const auto & sampleCounts :clusterTotals){
		popCollapse_->collapsed_.info_.updateCoi(sampleCounts.second);
	}
}

bool SampleCollapseCollection::excludeCommonlyLowFreqHaps(double lowFreqCutOff){
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	std::vector<uint32_t> lowFreqHaps;
	for(const auto clusPos : iter::range(popCollapse_->collapsed_.clusters_.size())){
		const auto & clus = popCollapse_->collapsed_.clusters_[clusPos];
		auto avgFrac = clus.getCumulativeFrac()/clus.sampInfos().size();
		//std::cout << clus.seqBase_.name_ << " avgFrac: " << avgFrac << std::endl;
		if(avgFrac < lowFreqCutOff){
			bool save = false;
			for(const auto & subClus : clus.reads_){
				if(subClus->seqBase_.frac_ > 2 * lowFreqCutOff){
					save = true;
				}
			}
			if(!save){
				lowFreqHaps.emplace_back(clusPos);
			}
		}
	}
	if(lowFreqHaps.empty()){
		return false;
	}
	std::unordered_map<std::string, std::vector<std::string>> hapsToBeRemovePerSample;
	for(const auto & lowFreqHap : lowFreqHaps){
		for(const auto & clus : popCollapse_->collapsed_.clusters_[lowFreqHap].reads_){
			hapsToBeRemovePerSample[clus->getOwnSampName()].emplace_back(oututSampClusToOldNameKey_[clus->getOwnSampName()][clus->seqBase_.name_]);
		}
	}
	for(const auto & sampClusters : hapsToBeRemovePerSample){
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampClusters.first);
		}
		std::vector<uint32_t> toBeExcluded;
		for(const auto clusPos : iter::range(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.size())){
			if(njh::in(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[clusPos].seqBase_.name_, sampClusters.second)){
				toBeExcluded.push_back(clusPos);
			}
		}
		std::sort(toBeExcluded.rbegin(), toBeExcluded.rend());
		for(const auto exclude : toBeExcluded){
			MetaDataInName filteredMeta;
			if(MetaDataInName::nameHasMetaData(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.name_)){
				filteredMeta = MetaDataInName(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.name_);
			}
			filteredMeta.addMeta("ExcludeCommonlyLowFreq", "TRUE");
			if (sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.isChimeric()) {
				filteredMeta.addMeta("ExcludeIsChimeric", "TRUE", true);
			}
			sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.resetMetaInName(filteredMeta);

			sampleCollapses_.at(sampClusters.first)->excluded_.clusters_.emplace_back(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude]);
			sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.erase(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.begin() + exclude);
		}
		sampleCollapses_.at(sampClusters.first)->updateAfterExclustion();
		sampleCollapses_.at(sampClusters.first)->renameClusters("fraction");

		if(!keepSampleInfoInMemory_){
			dumpSample(sampClusters.first);
		}
	}
	return true;
}


bool SampleCollapseCollection::excludeOneSampOnlyHaps(double fracCutOff){
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	std::vector<uint32_t> lowFreqHaps;
	std::unordered_map<std::string, VecStr> samplesWithUniqHaps;
	for(const auto & popClus : popCollapse_->collapsed_.clusters_){
		if(1 == popClus.sampleClusters().size()){
			//auto avgFrac = popClus.getCumulativeFrac()/popClus.sampInfos().size();
			auto avgFrac = popClus.getCumulativeFrac();
			if(avgFrac < fracCutOff){
				//this should just be one anyways
				for(const auto & clus : popClus.reads_){
					samplesWithUniqHaps[clus->getOwnSampName()].emplace_back(oututSampClusToOldNameKey_[clus->getOwnSampName()][clus->seqBase_.name_]);
				}
			}
		}
	}
	if(samplesWithUniqHaps.empty()){
		return false;
	}

	for(const auto & sampClusters : samplesWithUniqHaps){
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampClusters.first);
		}
		std::vector<uint32_t> toBeExcluded;
		for(const auto clusPos : iter::range(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.size())){
			if(njh::in(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[clusPos].seqBase_.name_, sampClusters.second)){
				toBeExcluded.push_back(clusPos);
			}
		}
		std::sort(toBeExcluded.rbegin(), toBeExcluded.rend());
		for(const auto exclude : toBeExcluded){
			MetaDataInName filteredMeta;
			if(MetaDataInName::nameHasMetaData(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.name_)){
				filteredMeta = MetaDataInName(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.name_);
			}
			filteredMeta.addMeta("ExcludeOneSampOnly", "TRUE");
			if (sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.isChimeric()) {
				filteredMeta.addMeta("ExcludeIsChimeric", "TRUE", true);
			}
			sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.resetMetaInName(filteredMeta);

			sampleCollapses_.at(sampClusters.first)->excluded_.clusters_.emplace_back(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude]);
			sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.erase(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.begin() + exclude);
		}
		sampleCollapses_.at(sampClusters.first)->updateAfterExclustion();
		sampleCollapses_.at(sampClusters.first)->renameClusters("fraction");
		if(!keepSampleInfoInMemory_){
			dumpSample(sampClusters.first);
		}
	}
	return true;
}

bool SampleCollapseCollection::excludeOneSampOnlyOneOffHaps(double fracCutOff, aligner & alignerObj){
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	std::vector<uint32_t> lowFreqHaps;
	std::unordered_map<std::string, VecStr> samplesWithUniqHaps;
	for(const auto & popClus : popCollapse_->collapsed_.clusters_){
		if(1 == popClus.sampleClusters().size()){
			//auto avgFrac = popClus.getCumulativeFrac()/popClus.sampInfos().size();
			auto avgFrac = popClus.getCumulativeFrac();
			if(avgFrac < fracCutOff){
				//this should just be one anyways
				for(const auto & clus : popClus.reads_){
					samplesWithUniqHaps[clus->getOwnSampName()].emplace_back(oututSampClusToOldNameKey_[clus->getOwnSampName()][clus->seqBase_.name_]);
				}
			}
		}
	}
	if(samplesWithUniqHaps.empty()){
		return false;
	}

	for(const auto & sampClusters : samplesWithUniqHaps){
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampClusters.first);
		}
		std::vector<uint32_t> toBeExcluded;
		for(const auto clusPos : iter::range(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.size())){
			if(njh::in(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[clusPos].seqBase_.name_, sampClusters.second)){
				bool oneOffAnother = false;
				for(const auto & clus : sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_){
					alignerObj.alignCacheGlobal(clus, sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[clusPos]);
					alignerObj.profileAlignment(clus, sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[clusPos],
							false, false, false);
					if(alignerObj.comp_.distances_.getNumOfEvents(true) <=1){
						oneOffAnother = true;
						break;
					}
				}
				if(oneOffAnother){
					toBeExcluded.push_back(clusPos);
				}
			}
		}
		std::sort(toBeExcluded.rbegin(), toBeExcluded.rend());
		for(const auto exclude : toBeExcluded){
			MetaDataInName filteredMeta;
			if(MetaDataInName::nameHasMetaData(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.name_)){
				filteredMeta = MetaDataInName(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.name_);
			}
			filteredMeta.addMeta("ExcludeOneSampOnlyOneOffAnother", "TRUE");
			if (sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.isChimeric()) {
				filteredMeta.addMeta("ExcludeIsChimeric", "TRUE", true);
			}
			sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude].seqBase_.resetMetaInName(filteredMeta);

			sampleCollapses_.at(sampClusters.first)->excluded_.clusters_.emplace_back(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_[exclude]);
			sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.erase(sampleCollapses_.at(sampClusters.first)->collapsed_.clusters_.begin() + exclude);
		}
		sampleCollapses_.at(sampClusters.first)->updateAfterExclustion();
		sampleCollapses_.at(sampClusters.first)->renameClusters("fraction");
		if(!keepSampleInfoInMemory_){
			dumpSample(sampClusters.first);
		}
	}
	return true;
}

void SampleCollapseCollection::dumpPopulation(const bfs::path & popDir, bool dumpTable) {
	if (popCollapse_) {
		auto popSubClusDir = njh::files::make_path(popDir, "clusters");
		//clear previous population dir or create if for the first time
		njh::files::makeDir(njh::files::MkdirPar(popDir.string(), true));
		njh::files::makeDir(njh::files::MkdirPar(popSubClusDir.string(), true));
		SeqIOOptions popOutOpts(
				njh::files::make_path(popDir,
						"PopSeqs" + inputOptions_.getOutExtension()),
				inputOptions_.outFormat_);

		SeqOutput::write(popCollapse_->collapsed_.clusters_, popOutOpts);
		if (!popCollapse_->collapsed_.clusters_.empty()
				&& "" != popCollapse_->collapsed_.clusters_.front().expectsString) {
			MapStrStr refCompInfos;
			for (const auto & clus : popCollapse_->collapsed_.clusters_) {
				refCompInfos[clus.seqBase_.name_] = clus.expectsString;
			}
			table refCompTab(refCompInfos, VecStr { "read", "bestExpected" });
			refCompTab.outPutContents(
					TableIOOpts(
							OutOptions(
									njh::files::join(popDir.string(), "refCompInfos.tab.txt")),
							"\t", true));

			MapStrStr aaTypedInfos;
			for (const auto & clus : popCollapse_->collapsed_.clusters_) {
				std::string h_aaTyped = "";
				if(clus.meta_.containsMeta("h_AATyped")){
					h_aaTyped = clus.meta_.getMeta("h_AATyped");
				}
				aaTypedInfos[clus.seqBase_.seq_] = h_aaTyped;
			}
			table aaTypedTab(aaTypedInfos, VecStr { "seq", "h_AATyped" });
			aaTypedTab.outPutContents(
					TableIOOpts(
							OutOptions(
									njh::files::join(popDir.string(), "aaTyped.tab.txt")),
							"\t", true));


		}
		std::unordered_map<std::string, double> readTotals;
		SeqIOOptions inputPopClustersOpts(
									njh::files::join(popSubClusDir.string(), "inputClusters").string()
											+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		inputPopClustersOpts.out_.overWriteFile_ = true;
		SeqOutput inputPopClusWriter(inputPopClustersOpts);

		inputPopClusWriter.openOut();
		for (const auto & clus : popCollapse_->collapsed_.clusters_) {
			for (const auto & subClus : clus.reads_) {
				readTotals[subClus->getOwnSampName()] += subClus->seqBase_.cnt_;

				MetaDataInName subClusMeta;
				bool hadMeta = false;
				if(MetaDataInName::nameHasMetaData(subClus->seqBase_.name_)){
					subClusMeta.processNameForMeta(subClus->seqBase_.name_, true);
					hadMeta = true;
				}
				auto writeCluster = subClus->seqBase_;
				subClusMeta.addMeta("TopClusterName", clus.seqBase_.name_);
				if(hadMeta){
					subClusMeta.resetMetaInName(writeCluster.name_);
				}else{
					subClusMeta.resetMetaInName(writeCluster.name_, writeCluster.name_.rfind("_f"));
				}
				inputPopClusWriter.write(writeCluster);
			}
		}
		table readNumsTab(readTotals, VecStr { "sample", "readTotal" });
		readNumsTab.sortTable("sample", false);
		readNumsTab.outPutContents(
				TableIOOpts(
						OutOptions(njh::files::join(popDir.string(), "readTotals.tab.txt")),
						"\t", true));
		if(dumpTable){
			printPopulationCollapseInfo(njh::files::make_path(popDir, "populationCluster.tab.txt.gz"));
		}
		popCollapse_ = nullptr;
	}
}

std::vector<seqInfo> SampleCollapseCollection::genOutPopSeqsPerSample() const{
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	std::vector<seqInfo> outseqs;
	for(const auto & seq : popCollapse_->collapsed_.clusters_){
		for(const auto & subSeq : seq.reads_){
			auto topSeqCopy = seq.seqBase_;
			topSeqCopy.name_ = subSeq->seqBase_.name_;
			//auto subSeqCopy = subSeq->seqBase_;
			auto sample = topSeqCopy.getOwnSampName();
			MetaDataInName subseqMeta;
			subseqMeta.addMeta("PopUID", seq.getStubName(true));
			subseqMeta.addMeta("sample", sample);
			subseqMeta.addMeta("readCount", subSeq->seqBase_.cnt_);
			subseqMeta.addMeta("p_name", popNames_.populationName_);
			if(nullptr != groupMetaData_){
				auto sampMeta = groupMetaData_->getMetaForSample(sample, getVectorOfMapKeys(groupMetaData_->groupData_));
				subseqMeta.addMeta(sampMeta, false);
			}
			subseqMeta.resetMetaInName(topSeqCopy.name_, topSeqCopy.name_.rfind("_f"));
			outseqs.emplace_back(topSeqCopy);
		}
	}
	return outseqs;
}

void SampleCollapseCollection::dumpPopulation(bool dumpTable) {
	dumpPopulation(populationOutputDir_, dumpTable);
}

void SampleCollapseCollection::loadInPreviousPop() {
	loadInPreviousPop(popNames_.samples_, populationOutputDir_);
}

void SampleCollapseCollection::loadInPreviousPop(const std::set<std::string> & samples){
	loadInPreviousPop(samples, populationOutputDir_);
}

void SampleCollapseCollection::loadInPreviousPop(const std::set<std::string> & samples, const bfs::path & popDir) {


	popCollapse_ = std::make_unique<collapse::populationCollapse>(
			popNames_.populationName_);

	table readNumsTab(
			njh::files::make_path(popDir,
					"readTotals.tab.txt").string(), "\t", true);
	std::unordered_map<std::string, sampInfo> sampInfos;
	for (const auto & row : readNumsTab.content_) {
		auto sampName = row[readNumsTab.getColPos("sample")];
		sampInfos[sampName] = sampInfo(sampName,
				njh::lexical_cast<double>(row[readNumsTab.getColPos("readTotal")]));
	}

	MapStrStr aaTyped;
	if (bfs::exists(njh::files::join(popDir.string(), "aaTyped.tab.txt"))) {
		table aaTypedTab(njh::files::join(popDir.string(), "aaTyped.tab.txt"),
				"\t", true);
		for (const auto & row : aaTypedTab.content_) {
			aaTyped[row[aaTypedTab.getColPos("seq")]] =
					row[aaTypedTab.getColPos("h_AATyped")];
		}
	}

	MapStrStr refCompInfos;
	if (bfs::exists(njh::files::join(popDir.string(), "refCompInfos.tab.txt"))) {
		table refCompTab(njh::files::join(popDir.string(), "refCompInfos.tab.txt"),
				"\t", true);
		for (const auto & row : refCompTab.content_) {
			refCompInfos[row[refCompTab.getColPos("read")]] =
					row[refCompTab.getColPos("bestExpected")];
		}
	}


	auto popInputClustersOpts = SeqIOOptions(
			njh::files::make_path(popDir,
					"clusters/inputClusters" + inputOptions_.getOutExtension()),
			inputOptions_.inFormat_, true);
	SeqInput inputReader(popInputClustersOpts);
	inputReader.openIn();
	std::unordered_map<std::string, std::vector<seqInfo>> inputPopClusters;
	seqInfo inputSeq;
	while(inputReader.readNextRead(inputSeq)){
		MetaDataInName inputMeta(inputSeq.name_);
		auto topClusterName = inputMeta.getMeta("TopClusterName");
		inputMeta.removeMeta("TopClusterName");
		if(inputMeta.meta_.size() == 0){
			MetaDataInName::removeMetaDataInName(inputSeq.name_);
		}else{
			inputMeta.resetMetaInName(inputSeq.name_);
		}
		inputPopClusters[topClusterName].emplace_back(inputSeq);
	}

	SeqIOOptions popInOpts = SeqIOOptions(getPopFinalHapsPath().string(),
			inputOptions_.inFormat_, true);
	SeqInput reader(popInOpts);
	reader.openIn();
	seqInfo seq;
	std::vector<sampleCluster> input;
	while (reader.readNextRead(seq)) {
		std::vector<seqInfo> subSeqs;
		std::set<std::string> subSamples;
//		auto subClustersPath = njh::files::make_path(popSubClusDir,
//				seq.name_ + inputOptions_.getOutExtension());
//		SeqIOOptions subClusInOpts = SeqIOOptions(subClustersPath.string(),
//				inputOptions_.inFormat_, true);
//		SeqInput subReader(subClusInOpts);
//		subReader.openIn();
//		seqInfo subSeq;
//		while (subReader.readNextRead(subSeq)) {
		for(auto & subSeq : inputPopClusters.at(seq.name_)){
			//a little hacky but has to be done for certain pipeline;
			std::string inputSampleName = subSeq.getOwnSampName();
			if(MetaDataInName::nameHasMetaData(subSeq.name_)){
				MetaDataInName meta(subSeq.name_);
				if(meta.containsMeta("samp")){
					inputSampleName = meta.getMeta("samp");
				}
			}
			if (njh::in(inputSampleName, samples)) {
				subSeq.cnt_ = subSeq.frac_
						* sampInfos.at(subSeq.getOwnSampName()).runReadCnt_;
				subSeqs.emplace_back(subSeq);
				subSamples.insert(subSeq.getOwnSampName());
				std::map<std::string, sampInfo> currentInfo { { subSeq.getOwnSampName(),
						sampInfos.at(subSeq.getOwnSampName()) } };
				input.emplace_back(sampleCluster(subSeq, currentInfo));
			}
		}
		if (!subSeqs.empty()) {
			std::map<std::string, sampInfo> currentInfos;
			for (const auto & samp : subSamples) {
				currentInfos[samp] = sampInfos[samp];
			}
			popCollapse_->collapsed_.clusters_.emplace_back(seq, subSeqs,
					currentInfos);
			popCollapse_->collapsed_.clusters_.back().updateSampInfosFracs();
			popCollapse_->collapsed_.clusters_.back().expectsString =
					refCompInfos[popCollapse_->collapsed_.clusters_.back().seqBase_.name_];
		}
	}
	popCollapse_->addInput(input);
	popCollapse_->updateCollapsedInfos();
	std::unordered_map<std::string, uint32_t> clusterTotals;
	for (const auto & clus : popCollapse_->collapsed_.clusters_) {
		for (const auto & subClus : clus.reads_) {
			++clusterTotals[subClus->getOwnSampName()];
		}
	}
	popCollapse_->collapsed_.info_.cois_.clear();
	for(const auto & sampleCounts :clusterTotals){
		popCollapse_->collapsed_.info_.updateCoi(sampleCounts.second);
	}

	for(auto & clus : popCollapse_->collapsed_.clusters_){
		clus.addMeta("h_AATyped", aaTyped[clus.seqBase_.seq_], true);
	}
}

void SampleCollapseCollection::renamePopWithSeqs(
		const std::vector<readObject> & otherPopSeqs, comparison allowableErrors) {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	popCollapse_->renameToOtherPopNames(otherPopSeqs, allowableErrors);
}

void SampleCollapseCollection::addRefMetaToName(
		const std::vector<readObject> & otherPopSeqs, comparison allowableErrors) {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	popCollapse_->addRefMetaToName(otherPopSeqs, allowableErrors);
}

void SampleCollapseCollection::comparePopToRefSeqs(
		const std::vector<readObject> & expectedSeqs, aligner & alignerObj) {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	popCollapse_->collapsed_.checkAgainstExpected(expectedSeqs, alignerObj, false);
}


void SampleCollapseCollection::checkForPopCollapseThrow(
		const std::string & funcName) const {
	if (nullptr == popCollapse_) {
		std::stringstream ss;
		ss << funcName << ", error popCollapse_ not loaded\n";
		throw std::runtime_error { ss.str() };
	}
}


void SampleCollapseCollection::printPopulationCollapseInfo(const bfs::path& fileName)const {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	auto popTab = genPopulationCollapseInfo();
	bool zip = njh::endsWith(fileName.string(), ".gz");
	popTab.outPutContents(TableIOOpts(OutOptions(fileName.string(), zip ? ".gz" : ".txt"), "\t", true));
}

table SampleCollapseCollection::genPopulationCollapseInfo() const {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	VecStr header = toVecStr("h_popUID",
			populationCollapse::getPopInfoHeaderVec());
	//add gorup meta data if loaded
	if (groupMetaLoaded()) {
		for (const auto & group : groupMetaData_->groupData_) {
			addOtherVec(header,
					toVecStr("g_" + group.first + "_subGroupCounts",
							"g_" + group.first + "_subGroupFracs"));
		}
	}
	addOtherVec(header, toVecStr(
			sampleCluster::getPopHapInfoHeaderVec(), "bestExpected"));
	table ret(header);
	for (const auto& clus : popCollapse_->collapsed_.clusters_) {
		auto row = toVecStr(clus.getStubName(false),
				popCollapse_->getPopInfoVec());
		//add gorup meta data if loaded
		if(nullptr != groupMetaData_){
			std::vector<std::string> currentSamples;
			for(const auto & samp : clus.sampInfos()){;
				if(samp.second.numberOfClusters_ > 0){
					currentSamples.emplace_back(samp.first);
				}
			}
			auto groupPopInfos = groupMetaData_->getGroupPopInfos(currentSamples);
			for(const auto & info : groupPopInfos){
				addOtherVec(row, toVecStr(info.groupCountsStr(), info.groupFracsStr()));
			}
		}
		addOtherVec(row,toVecStr(clus.getPopHapInfoVec(popCollapse_->collapsed_.info_.totalReadCount_,
						popCollapse_->collapsed_.numOfReps()),
				clus.expectsString) );
		ret.addRow(row);
	}
	return ret;
}

void SampleCollapseCollection::printSampleCollapseInfo(const bfs::path& fileName){
	auto sampTab = genSampleCollapseInfo(popNames_.samples_);
	bool zip = njh::endsWith(fileName.string(), ".gz");
	sampTab.outPutContents(TableIOOpts(OutOptions(fileName.string(), zip ? ".gz" : ".txt"), "\t", true));
}




void SampleCollapseCollection::printAllSubClusterInfo(const OutOptions& outOpts, bool skipExcludeReadCntCutOff) {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	OutputStream allSubClusterInfo(outOpts);
	uint32_t maxRunCount = 0;
	std::set<std::string> allMetaFields;
	for (const auto & sampName : popNames_.samples_) {
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampName);
		}
		if (njh::in(sampName, lowRepCntSamples_)
				|| sampleCollapses_[sampName]->collapsed_.info_.totalReadCount_
						< preFiltCutOffs_.sampleMinReadCount) {
			continue;
		}
		for(const auto & excluded : sampleCollapses_.at(sampName)->excluded_.clusters_){
			if(MetaDataInName::nameHasMetaData(excluded.seqBase_.name_)){
				MetaDataInName excmeta(excluded.seqBase_.name_);
				for(const auto & mf : excmeta.meta_){
					allMetaFields.emplace(mf.first);
				}
			}
		}
		for(const auto & collapsed : sampleCollapses_.at(sampName)->collapsed_.clusters_){
			if(MetaDataInName::nameHasMetaData(collapsed.seqBase_.name_)){
				MetaDataInName excmeta(collapsed.seqBase_.name_);
				for(const auto & mf : excmeta.meta_){
					allMetaFields.emplace(mf.first);
				}
			}
		}
		if (maxRunCount
				< sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size()) {
			maxRunCount =
					sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size();
		}
		if(!keepSampleInfoInMemory_){
			clearSample(sampName);
		}
	}
	std::string delim = "\t";
	allSubClusterInfo << "s_name";
	allSubClusterInfo << delim << "h_PopUID";
	allSubClusterInfo << delim << "c_Name";
	allSubClusterInfo << delim << "c_ReadCnt";
	allSubClusterInfo << delim << "c_RunCnt";
	allSubClusterInfo << delim << "c_seq";
	allSubClusterInfo << delim << "c_ClusterCnt";
	allSubClusterInfo << delim << "c_IncludedInFinalAnalysis";
	allSubClusterInfo << delim << "c_bestExpected";

	for(const auto & mf : allMetaFields){
		allSubClusterInfo << delim << "c_" << mf;
	}
	for(const auto runNum : iter::range<uint32_t>(1, maxRunCount + 1)){
		allSubClusterInfo << delim << "R" << runNum << "_Name";
		allSubClusterInfo << delim << "R" << runNum << "_ReadCnt";
		allSubClusterInfo << delim << "R" << runNum << "_RepTotalCnt";
	}

	allSubClusterInfo << std::endl;
	for (const auto& sampName : popNames_.samples_) {
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampName);
		}
		if (njh::in(sampName, lowRepCntSamples_)
				|| sampleCollapses_[sampName]->collapsed_.info_.totalReadCount_
						< preFiltCutOffs_.sampleMinReadCount) {
			continue;
		}
		auto sampPtr = sampleCollapses_.at(sampName);
		for (const auto clusPos : iter::range(sampPtr->collapsed_.clusters_.size())) {
			const auto & clus = sampPtr->collapsed_.clusters_[clusPos];
			allSubClusterInfo  << sampName
					<< delim << popCollapse_->collapsed_.clusters_[popCollapse_->collapsed_.subClustersPositions_.at(
							clus.getStubName(true))].seqBase_.name_;
			allSubClusterInfo
						<< delim << clus.seqBase_.name_
						<< delim << clus.seqBase_.cnt_
						<< delim << clus.numberOfRuns()
						<< delim << clus.seqBase_.seq_
						<< delim << clus.reads_.size();
			allSubClusterInfo << delim << "TRUE";
			allSubClusterInfo << delim << clus.expectsString;
			MetaDataInName seqMeta;
			if(MetaDataInName::nameHasMetaData(clus.seqBase_.name_)){
				seqMeta = MetaDataInName(clus.seqBase_.name_);
			}
			for(const auto & mField : allMetaFields){
				allSubClusterInfo << delim;
				if(seqMeta.containsMeta(mField)){
					allSubClusterInfo << seqMeta.getMeta(mField);
				}else{
					allSubClusterInfo << "NA";
				}
			}
			for (const auto &info : sampPtr->input_.info_.infos_) {
				allSubClusterInfo << delim;
				auto search = clus.sampInfos().find(info.first);
				if (search == clus.sampInfos().end() || search->second.readCnt_ == 0) {
					allSubClusterInfo << info.first;
					allSubClusterInfo << delim << 0;
					allSubClusterInfo << delim << info.second.runReadCnt_;
				} else {
					allSubClusterInfo << info.first;
					allSubClusterInfo << delim << clus.sampInfos().at(info.first).readCnt_;
					allSubClusterInfo << delim << info.second.runReadCnt_;
				}
			}
			if(clus.sampInfos().size() < maxRunCount){
				for(uint32_t i = 0; i < maxRunCount - clus.sampInfos().size(); ++i){
					allSubClusterInfo << delim << delim << delim;
				}
			}
			allSubClusterInfo << std::endl;
		}
		for (const auto clusPos : iter::range(sampPtr->excluded_.clusters_.size())) {
			const auto & clus = sampPtr->excluded_.clusters_[clusPos];

			{
				MetaDataInName seqMeta;
				if(MetaDataInName::nameHasMetaData(clus.seqBase_.name_)){
					seqMeta = MetaDataInName(clus.seqBase_.name_);
				}
				if(skipExcludeReadCntCutOff && seqMeta.containsMeta("ExcludeFailedReadCutOff")){
					continue;
				}
			}
			allSubClusterInfo  << sampName
					;
			std::string matchingPopName;
			for(const auto & pop : popCollapse_->collapsed_.clusters_){
				if(pop.seqBase_.seq_ == clus.seqBase_.seq_){
					matchingPopName = pop.seqBase_.name_;
				}
			}
			allSubClusterInfo << delim << matchingPopName;

			allSubClusterInfo
						<< delim << clus.seqBase_.name_
						<< delim << clus.seqBase_.cnt_
						<< delim << clus.numberOfRuns()
						<< delim << clus.seqBase_.seq_
						<< delim << clus.reads_.size();
			allSubClusterInfo << delim << "FALSE";
			allSubClusterInfo << delim << clus.expectsString;
			MetaDataInName seqMeta;
			if(MetaDataInName::nameHasMetaData(clus.seqBase_.name_)){
				seqMeta = MetaDataInName(clus.seqBase_.name_);
			}
			for(const auto & mField : allMetaFields){
				allSubClusterInfo << delim;
				if(seqMeta.containsMeta(mField)){
					allSubClusterInfo << seqMeta.getMeta(mField);
				}else{
					allSubClusterInfo << "NA";
				}
			}
			for (const auto &info : sampPtr->input_.info_.infos_) {
				allSubClusterInfo << delim;
				auto search = clus.sampInfos().find(info.first);
				if (search == clus.sampInfos().end() || search->second.readCnt_ == 0) {
					allSubClusterInfo << info.first;
					allSubClusterInfo << delim << 0;
					allSubClusterInfo << delim << info.second.runReadCnt_;
				} else {
					allSubClusterInfo << info.first;
					allSubClusterInfo << delim << clus.sampInfos().at(info.first).readCnt_;
					allSubClusterInfo << delim << info.second.runReadCnt_;
				}
			}
			if(clus.sampInfos().size() < maxRunCount){
				for(uint32_t i = 0; i < maxRunCount - clus.sampInfos().size(); ++i){
					allSubClusterInfo << delim << delim << delim;
				}
			}
			allSubClusterInfo << std::endl;
		}
		if(!keepSampleInfoInMemory_){
			clearSample(sampName);
		}
	}
}

table SampleCollapseCollection::genSampleCollapseInfo(
		const std::set<std::string> & samples) {
	/**@todo currently a hot mess below, still needs to be redone, just left over code from the old days
	 *
	 */
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	std::string delim = "\t";
	bool checkingExpected = true;
	uint32_t maxRunCount = 0;
	std::vector<VecStr> rows;
	for (const auto& sampName : samples) {
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampName);
		}
		if (njh::in(sampName, lowRepCntSamples_)
				|| sampleCollapses_[sampName]->collapsed_.info_.totalReadCount_
						< preFiltCutOffs_.sampleMinReadCount) {
			continue;
		}
		if (maxRunCount
				< sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size()) {
			maxRunCount =
					sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size();
		}
		if(!keepSampleInfoInMemory_){
			clearSample(sampName);
		}
	}
	for (const auto& sampName : samples) {
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampName);
		}
		if (njh::in(sampName, lowRepCntSamples_)
				|| sampleCollapses_[sampName]->collapsed_.info_.totalReadCount_
						< preFiltCutOffs_.sampleMinReadCount) {
			continue;
		}
		for (const auto clusPos : iter::range(
				sampleCollapses_.at(sampName)->collapsed_.clusters_.size())) {
			const auto & clus =
					sampleCollapses_.at(sampName)->collapsed_.clusters_[clusPos];
			std::stringstream rowStream;
			rowStream  << sampName
					<< delim << popCollapse_->populationName_
					<< delim << popCollapse_->collapsed_.clusters_[popCollapse_->collapsed_.subClustersPositions_.at(
							clus.getStubName(true))].getPopInfo(
							popCollapse_->collapsed_.info_.totalReadCount_,
							popCollapse_->collapsed_.info_.numberOfClusters_,
							popCollapse_->collapsed_.info_.infos_.size(), delim)
							<< delim << sampleCollapses_.at(sampName)->getSimpleSampInfo(delim);
			if(nullptr != groupMetaData_){
				for(const auto & group : groupMetaData_->groupData_){
					rowStream << delim << group.second->getGroupForSample(sampName);
				}
			}
			rowStream << delim << clusPos
							<< delim << clus.seqBase_.name_
							<< delim << clus.getClusterInfo(delim)
							<< delim << clus.getRepsInfo(
							sampleCollapses_.at(sampName)->input_.info_.infos_,
							sampleCollapses_.at(sampName)->excluded_.info_.infos_,
							sampleCollapses_.at(sampName)->collapsed_.info_.infos_,
							maxRunCount,
							checkingExpected, delim);
			rows.emplace_back(tokenizeString(rowStream.str(), delim, true));
		}
		if(!keepSampleInfoInMemory_){
			clearSample(sampName);
		}
	}
	std::stringstream headerStream;
	headerStream << "s_Sample"
			<< delim << "p_name"
			<< delim << sampleCluster::getPopInfoHeader(delim)
			<< delim << collapse::sampleCollapse::getSimpleSampInfoHeader(delim);
	if(nullptr != groupMetaData_){
		for(const auto & group : groupMetaData_->groupData_){
			headerStream << delim << "s_" << group.first;
		}
	}
	headerStream << delim << "c_clusterID"
			<< delim << "c_name"
			<< delim << sampleCluster::getClusterInfoHeader(delim) << delim
			<< sampleCluster::getRepsInfoHeader(maxRunCount, checkingExpected, delim);
	table ret(tokenizeString(headerStream.str(), delim));
	//std::cout << headerStream.str() << std::endl;
	ret.addRows(rows);
	return ret;
}

void SampleCollapseCollection::symlinkInSampleFinals() const {
	auto finalDir = njh::files::make_path(masterOutputDir_, "final");
	njh::files::makeDir(njh::files::MkdirPar(finalDir.string(), true));
	for (const auto & samp : popNames_.samples_) {
		auto sampFinalHapFile = getSampleFinalHapsPath(samp);
		if(bfs::exists(sampFinalHapFile)){
			bfs::create_symlink(bfs::relative(sampFinalHapFile, finalDir), njh::files::make_path(finalDir, sampFinalHapFile.filename()));
		}
	}
}

bool SampleCollapseCollection::groupMetaLoaded() const{
	return nullptr != groupMetaData_;
}

void SampleCollapseCollection::createGroupInfoFiles(){
	if(nullptr != groupMetaData_){
		auto groupsTopDir = njh::files::make_path(masterOutputDir_, "groups");
		njh::files::makeDir(njh::files::MkdirPar(groupsTopDir.string(), true));
		std::ofstream groupMetaJsonFile;
		openTextFile(groupMetaJsonFile,
				njh::files::make_path(groupsTopDir, "groupMetaData.json").string(),
				".json", false, true);
		groupMetaJsonFile << groupMetaData_->toJson() << std::endl;
		for(const auto & group : groupMetaData_->groupData_){
			auto mainGroupDir = njh::files::make_path(groupsTopDir, group.first);
			njh::files::makeDir(njh::files::MkdirPar(mainGroupDir.string(), true));
			std::unordered_map<std::string, table> popTabs;
			std::unordered_map<std::string, table> sampTabs;
			std::unordered_map<std::string, table> hapIdTabs;


			for(const auto & subGroup : group.second->subGroupToSamples_){
				loadInPreviousPop(subGroup.second);
				popTabs[subGroup.first] = genPopulationCollapseInfo();
				sampTabs[subGroup.first] = genSampleCollapseInfo(subGroup.second);
				hapIdTabs[subGroup.first] = genHapIdTable(subGroup.second);
				popCollapse_ = nullptr;
			}

			std::unordered_map<std::string, VecStr> popUids;
			for(const auto & pop : popTabs){
				popUids[pop.first] = pop.second.getColumnLevels("h_popUID");
			}
			std::unordered_map<std::string, VecStr> uniquePopUids;
			for(const auto & popUid : popUids){
				for(const auto & uid : popUid.second){
					bool unique = true;
					for(const auto & otherUid : popUids){
						if(popUid.first != otherUid.first){
							if(njh::in(uid, otherUid.second)){
								unique = false;
								break;
							}
						}
					}
					if(unique){
						uniquePopUids[popUid.first].emplace_back(uid);
					}
				}
			}
			for(auto & popTab : popTabs){
				popTab.second.addColumn({group.first + ":" + popTab.first}, "g_GroupName");
				popTab.second.addColumn({njh::conToStr(uniquePopUids[popTab.first])}, "g_hapsFoundOnlyInThisGroup");
				popTab.second.addColumn({estd::to_string(uniquePopUids[popTab.first].size())}, "p_TotalUniqueHaplotypes");
			}
			for(auto & sampTab : sampTabs){
				auto currentHeader = sampTab.second.columnNames_;
				sampTab.second.addColumn({group.first + ":" + sampTab.first}, "g_GroupName");
				removeElement(currentHeader, std::string("s_Sample"));
				VecStr outHeader{"s_Sample", "g_GroupName"};
				addOtherVec(outHeader, currentHeader);
				sampTab.second = sampTab.second.getColumns(outHeader);
			}
			VecStr groupInfoColNames { "g_GroupName", "p_TotalInputReadCnt",
					"p_TotalInputClusterCnt", "p_TotalPopulationSampCnt",
					"p_TotalHaplotypes", "p_TotalUniqueHaplotypes", "p_ExpectedHeterozygosity",
					"p_meanCoi", "p_medianCoi", "p_minCoi",
					"p_maxCoi", "g_hapsFoundOnlyInThisGroup"};
			//info on all the sub groups
			table outTab(groupInfoColNames);

			for(const auto & popTab : popTabs){
				auto subGroupDir = njh::files::make_path(mainGroupDir, popTab.first);
				njh::files::makeDir(njh::files::MkdirPar(subGroupDir.string()));
				popTab.second.outPutContents(TableIOOpts(OutOptions(njh::files::make_path(subGroupDir,"popFile.tab.txt.gz")), "\t", true));
				sampTabs.at(popTab.first).outPutContents(TableIOOpts(OutOptions(njh::files::make_path(subGroupDir,"sampFile.tab.txt.gz")), "\t", true));
				hapIdTabs.at(popTab.first).outPutContents(TableIOOpts(OutOptions(njh::files::make_path(subGroupDir,"hapIdTable.tab.txt.gz")), "\t", true));
				std::ofstream subGroupMetaJsonFile;
				openTextFile(subGroupMetaJsonFile,
						njh::files::make_path(subGroupDir, "subGroupNamesData.json").string(),
						".json", false, true);
				auto popNames = popTab.second.getColumn("h_popUID");
				auto sampNames = sampTabs.at(popTab.first).getColumnLevels("s_Name");
				Json::Value nameMetaData;
				nameMetaData["popUIDs"] = njh::json::toJson(popNames);
				nameMetaData["sampNames"] = njh::json::toJson(sampNames);
				subGroupMetaJsonFile << nameMetaData << std::endl;
				outTab.rbind(popTab.second.getColumns(groupInfoColNames), false);
			}
			outTab = outTab.getUniqueRows();
			outTab.sortTable("g_GroupName", false);
			outTab.outPutContents(TableIOOpts(OutOptions(njh::files::make_path(mainGroupDir,"groupInfo.tab.txt.gz")), "\t", true));
		}

		bool verbose = false;
		if(verbose){
			std::cout << "Missing meta data for the following samples:" << std::endl;
			std::cout << njh::conToStr(groupMetaData_->missingMetaForSamples_, "\n") << std::endl;
			std::cout << "Missing clustering data for the following samples:" << std::endl;
			std::cout << njh::conToStr(groupMetaData_->missingSamples_, "\n") << std::endl;
		}
	}
}

void SampleCollapseCollection::outputRepAgreementInfo() {
	bfs::path sampRepInfoDir = njh::files::makeDir(masterOutputDir_.string(),
			njh::files::MkdirPar("sampRepAgreementInfo", true));
	table rmseTab(VecStr { "sampleName", "RMSE" });
	table repInfoTab(
			VecStr { "sampleName", "clusterName", "repName", "fraction" });
	for (const auto & sampleName : popNames_.samples_) {
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(sampleName);
		}
		if (njh::in(sampleName, lowRepCntSamples_)
				|| sampleCollapses_[sampleName]->collapsed_.info_.totalReadCount_
						< preFiltCutOffs_.sampleMinReadCount) {
			continue;
		}
		auto & samp = *(sampleCollapses_[sampleName]);
		auto repInfoforSample = samp.collapsed_.getReplicateInfo();
		for (const auto & row : repInfoforSample.content_) {
			VecStr addingRow { sampleName };
			addOtherVec(addingRow, row);
			repInfoTab.content_.emplace_back(addingRow);
		}
		if (!samp.collapsed_.clusters_.empty()) {
			rmseTab.content_.emplace_back(
					toVecStr(sampleName, samp.collapsed_.getRMSE()));
		}
		if(!keepSampleInfoInMemory_){
			clearSample(sampleName);
		}
	}
	repInfoTab.outPutContents(
			TableIOOpts(
					OutOptions(sampRepInfoDir.string() + "replicatesFractions", ".tab.txt"), "\t",
					repInfoTab.hasHeader_));
	rmseTab.sortTable("RMSE", true);
	rmseTab.outPutContents(
			TableIOOpts(OutOptions(sampRepInfoDir.string() + "RMSE.tab.txt", ".tab.txt"), "\t",
					rmseTab.hasHeader_));
}


double SampleCollapseCollection::estimateNumberOfStrains(const std::set<std::string> & samples){
	/**@todo check to see if previous pop has already been loaded with these samples, can avoid loading again */

	loadInPreviousPop(samples);
	uint32_t numOfSingletons = 0;
	uint32_t numOfDoublets = 0;
	double totalNum = popCollapse_->collapsed_.clusters_.size();
	for(const auto & clus : popCollapse_->collapsed_.clusters_){
		if(1 == clus.numberOfRuns()){
			++numOfSingletons;
		}else if (2 == clus.numberOfRuns()){
			++numOfDoublets;
		}
	}
	if(0 == numOfSingletons){
		numOfSingletons = 1;
	}
	if(0 == numOfDoublets){
		numOfDoublets = 1;
	}
	//double variance = numOfDoublets * (std::pow(numOfSingletons/(numOfDoublets * 4.0), 4) + std::pow(numOfSingletons/(numOfDoublets), 3),std::pow(numOfSingletons/(numOfDoublets * 2.0), 2) );
	return totalNum + std::pow(numOfSingletons,2.0)/(2.0 * numOfDoublets);
}

double SampleCollapseCollection::calculatePIE(const std::set<std::string> & samples){
	/**@todo check to see if previous pop has already been loaded with these samples, can avoid loading again */
	loadInPreviousPop(samples);
	double sum = 0;
	for(const auto & clus : popCollapse_->collapsed_.clusters_){
		sum += clus.numberOfRuns() * (clus.numberOfRuns() - 1);
	}
	//original measure of dominace
	//double numOfSamples = len(samples);
	double dominance = sum/(len(samples) * (len(samples ) - 1));
	return 1 - dominance;
}


void SampleCollapseCollection::createCoreJsonFile() const{
	auto coreJson = toJsonCore();
	std::ofstream outCoreJsonFile;
	OutOptions outOpts(njh::files::make_path(masterOutputDir_, "coreInfo.json"));
	outOpts.overWriteFile_ = true;
	openTextFile(outCoreJsonFile, outOpts);
	outCoreJsonFile << coreJson << std::endl;
}

Json::Value SampleCollapseCollection::toJsonCore() const{
	Json::Value ret;
	ret["inputOptions_"] = njh::json::toJson(inputOptions_);
	ret["masterInputDir_"] = njh::json::toJson(njh::files::normalize(masterInputDir_));
	ret["masterOutputDir_"] = njh::json::toJson(njh::files::normalize(masterOutputDir_));
	ret["samplesOutputDir_"] = njh::json::toJson(njh::files::normalize(samplesOutputDir_));
	ret["populationOutputDir_"] = njh::json::toJson(njh::files::normalize(populationOutputDir_));
	ret["popNames_"] = njh::json::toJson(popNames_);
	ret["passingSamples_"] = njh::json::toJson(passingSamples_);
	ret["lowRepCntSamples_"] = njh::json::toJson(lowRepCntSamples_);
	ret["preFiltCutOffs_"] = njh::json::toJson(preFiltCutOffs_);
	return ret;
}


table SampleCollapseCollection::genHapIdTable(){
	return genHapIdTable(popNames_.samples_);
}

table SampleCollapseCollection::genHapIdTable(const std::set<std::string> & samples){

	if(nullptr == popCollapse_){
		loadInPreviousPop(samples);
	}
	VecStr popUidNames;
	for(const auto & clus : popCollapse_->collapsed_.clusters_){
		popUidNames.emplace_back(clus.getStubName(true));
	}
	njh::sort(popUidNames);
	std::vector<std::vector<std::string>> inputContent = std::vector<std::vector<std::string>>{popUidNames.size(), std::vector<std::string>{}};
	for(const auto popPos : iter::range(popUidNames.size())){
		inputContent[popPos].emplace_back(popUidNames[popPos]);
	}
	table ret(inputContent, VecStr{"#PopUID"});
	for (const auto & samp : samples) {
		if(!keepSampleInfoInMemory_){
			setUpSampleFromPrevious(samp);
		}
		if(njh::in(samp, lowRepCntSamples_) || sampleCollapses_[samp]->collapsed_.info_.totalReadCount_ < preFiltCutOffs_.sampleMinReadCount
				){
			continue;
		}
		std::unordered_map<std::string, double> popUidToFrac;
		for (const auto clusPos : iter::range(
				sampleCollapses_.at(samp)->collapsed_.clusters_.size())) {
			const auto & clus =
					sampleCollapses_.at(samp)->collapsed_.clusters_[clusPos];
			auto popUid = popCollapse_->collapsed_.clusters_[popCollapse_->collapsed_.subClustersPositions_.at(
										clus.getStubName(true))].getStubName(true);
			popUidToFrac[popUid] = clus.seqBase_.frac_;
		}
		VecStr addingCol;
		for(const auto & pop : popUidNames){
			auto search = popUidToFrac.find(pop);
			if(popUidToFrac.end() != search){
				addingCol.emplace_back(estd::to_string(search->second));
			}else{
				addingCol.emplace_back("0");
			}
		}
		ret.addColumn(addingCol, samp);
		if(!keepSampleInfoInMemory_){
			clearSample(samp);
		}
	}
	return ret;
}

void SampleCollapseCollection::excludeOnFrac(const std::string & sampleName,
		const std::unordered_map<std::string, double> & customCutOffsMap,
		bool fracExcludeOnlyInFinalAverageFrac) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	if (fracExcludeOnlyInFinalAverageFrac) {
		sampleCollapses_.at(sampleName)->excludeFraction(
				customCutOffsMap.at(sampleName), true);
	} else {
		sampleCollapses_.at(sampleName)->excludeFractionAnyRep(
				customCutOffsMap.at(sampleName), true);
	}
}

}  // namespace collapse
}  // namespace njhseq

