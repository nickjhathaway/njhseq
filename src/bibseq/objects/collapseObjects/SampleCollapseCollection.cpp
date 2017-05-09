/*
 * SampleCollapseCollection.cpp
 *
 *  Created on: Jul 24, 2016
 *      Author: nick
 */

#include "SampleCollapseCollection.hpp"
#include "bibseq/objects/seqObjects/Clusters/clusterUtils.hpp"


namespace bibseq {

namespace collapse {



SampleCollapseCollection::SampleCollapseCollection(SeqIOOptions inputOptions,
		const bfs::path & inputDir,
		const bfs::path & outputDir,
		const PopNamesInfo & popNames,
		uint32_t clusterSizeCutOff) :
		inputOptions_(inputOptions),
		masterInputDir_(inputDir),
		masterOutputDir_(outputDir),
		popNames_(popNames),
		clusterSizeCutOff_(clusterSizeCutOff) {
	samplesOutputDir_ = bib::files::make_path(masterOutputDir_, "samplesOutput");
	populationOutputDir_ = bib::files::make_path(masterOutputDir_, "population");
	bib::files::makeDir(bib::files::MkdirPar(samplesOutputDir_.string(), true));

}

SampleCollapseCollection::SampleCollapseCollection(const Json::Value & coreJson){

	bib::json::MemberChecker checker(coreJson);

	checker.failMemberCheckThrow(VecStr { "clusterSizeCutOff_", "inputOptions_",
			"masterInputDir_", "masterOutputDir_", "popNames_",
			"populationOutputDir_", "samplesOutputDir_" }, __PRETTY_FUNCTION__);

	inputOptions_ = SeqIOOptions(coreJson["inputOptions_"].toStyledString());
	masterInputDir_ = coreJson["masterInputDir_"].asString();
	masterOutputDir_ = coreJson["masterOutputDir_"].asString();
	popNames_ = PopNamesInfo(coreJson["popNames_"]["populationName_"].asString(),
			bib::json::jsonArrayToSet<std::string>(coreJson["popNames_"]["samples_"],
					[](const Json::Value & val){ return val.asString();}));
	clusterSizeCutOff_ = coreJson["clusterSizeCutOff_"].asUInt();

	populationOutputDir_ = coreJson["populationOutputDir_"].asString();
	samplesOutputDir_ = coreJson["samplesOutputDir_"].asString();
	//load in group meta data if it is available
	if(bfs::exists(bib::files::make_path(masterOutputDir_, "groups", "groupMetaData.json"))){
		groupMetaData_ = std::make_unique<MultipleGroupMetaData>(
				MultipleGroupMetaData::fromJson(
						bib::json::parseFile(
								bib::files::make_path(masterOutputDir_, "groups",
										"groupMetaData.json").string())));

		groupDataPaths_ = std::make_unique<AllGroupDataPaths>(
				bib::files::make_path(bib::files::normalize(masterOutputDir_), "groups"),
				groupMetaData_);
	}
}


void SampleCollapseCollection::addGroupMetaData(
		const bfs::path & groupingsFile) {

	groupMetaData_ = std::make_unique<MultipleGroupMetaData>(
			bib::files::normalize(groupingsFile), popNames_.samples_);
	/**@todo add some output to output directory that states the samples missing and missing meta*/
	groupDataPaths_ = std::make_unique<AllGroupDataPaths>(
			bib::files::make_path(bib::files::normalize(masterOutputDir_), "groups"),
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
		ss << "Options are: " << bib::conToStr(popNames_.samples_, ",") << "\n";
		throw std::runtime_error { ss.str() };
	}
}

void SampleCollapseCollection::setUpSample(const std::string & sampleName,
		aligner & alignerObj, const collapser & collapserObj,
		const ChimeraOpts & chiOpts) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	auto sampleDir = bib::files::make_path(masterInputDir_, sampleName);
	auto analysisFiles = bib::files::listAllFiles(sampleDir.string(), true, {
			std::regex { "^" + inputOptions_.firstName_.string() + "$" } }, 2);

	std::vector<RepFile> repFiles;
	for (const auto & af : analysisFiles) {
		auto fileToks = bib::tokenizeString(
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
		const std::vector<RepFile> & analysisFiles,
		aligner & alignerObj,
		const collapser & collapserObj,
		const ChimeraOpts & chiOpts){
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);

	std::vector<std::vector<cluster>> inputClusters;
	for (const auto & repf : analysisFiles) {
		auto sampleOpts = inputOptions_;
		sampleOpts.firstName_ = repf.repFnp_.string();
		SeqInput reader(sampleOpts);

		std::vector<cluster> clusters = baseCluster::convertVectorToClusterVector<
				cluster>(reader.readAllReads<seqInfo>());
		readVecSorter::sortReadVector(clusters, "totalCount");
		// consider adding the sample name in the name as well
		if(repf.reNameInput_){
			renameReadNamesNewClusters(clusters, repf.repName_, true, true, false);
		}

		if (chiOpts.checkChimeras_) {
			collapserObj.markChimeras(clusters, alignerObj, chiOpts);
		}
		clusterVec::allSetFractionClusters(clusters);
		if(repf.reNameInput_){
			readVec::allUpdateName(clusters);
		}
		inputClusters.emplace_back(clusters);
	}
	{
		std::lock_guard<std::mutex> lock(mut_);
		if (bib::in(sampleName, sampleCollapses_)) {
			sampleCollapses_[sampleName] = std::make_shared<sampleCollapse>(
					inputClusters, sampleName, clusterSizeCutOff_);
		} else {
			sampleCollapses_.emplace(sampleName,
					std::make_shared<sampleCollapse>(inputClusters, sampleName,
							clusterSizeCutOff_));
		}
	}
}

bfs::path SampleCollapseCollection::getSampleFinalHapsPath(
		const std::string & sampleName) const {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	return bib::files::make_path(samplesOutputDir_, sampleName, "final",
			sampleName + inputOptions_.getOutExtension());
}

bfs::path SampleCollapseCollection::getPopFinalHapsPath() const {
	return bib::files::make_path(populationOutputDir_,
			"PopSeqs" + inputOptions_.getOutExtension());

}


bfs::path SampleCollapseCollection::getPopInfoPath() const {
	return bib::files::make_path(populationOutputDir_,
			"populationCluster.tab.txt");
}

bfs::path SampleCollapseCollection::getSampInfoPath() const {
	/**@todo need to make this standard, currently set by SeekDeep processClusters*/
	return bib::files::make_path(masterOutputDir_, "selectedClustersInfo.tab.txt");
}
bfs::path SampleCollapseCollection::getHapIdTabPath() const {
	return bib::files::make_path(masterOutputDir_,
			"hapIdTable.tab.txt");
}



void SampleCollapseCollection::setUpSampleFromPrevious(
		const std::string & sampleName) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	auto sampPtr = std::make_shared<sampleCollapse>(sampleName);
	auto & samp = *sampPtr;
	auto sampleDir = bib::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName);

	std::map<std::string, double> allSampCounts;
	std::map<std::string, sampInfo> allSampInfos;

	std::map<std::string, double> finalSampCounts;
	std::map<std::string, sampInfo> finalSampInfos;

	auto finalClustersFiles = bib::files::filesInFolder(
			bib::files::make_path(sampleDir, "final/clusters"));
	auto excludedClustersFiles = bib::files::filesInFolder(
			bib::files::make_path(sampleDir, "excluded/clusters"));
	seqInfo seq;
	for (const auto & finalClustersFile : finalClustersFiles) {
		auto finalClustersFileOpts = SeqIOOptions(finalClustersFile.string(),
				inputOptions_.inFormat_, true);
		if(finalClustersFileOpts.inExists()){
			SeqInput reader(finalClustersFileOpts);
			reader.openIn();
			while (reader.readNextRead(seq)) {
				allSampCounts[seq.getOwnSampName()] += seq.cnt_;
				finalSampCounts[seq.getOwnSampName()] += seq.cnt_;
			}
		}
	}

	MapStrStr refCompInfosCollapsed;
	MapStrStr refCompInfosExcluded;
	if (bfs::exists(
			bib::files::join(sampleDir.string(), "refCompInfos.tab.txt"))) {
		table refCompTab(
				bib::files::join(sampleDir.string(), "refCompInfos.tab.txt"), "\t",
				true);
		//std::cout << __PRETTY_FUNCTION__ << " " << sampleName << std::endl;
		for (const auto & row : refCompTab.content_) {
			//std::cout << "\t" << row[refCompTab.getColPos("read")] << " :  " <<  row[refCompTab.getColPos("bestExpected")] << std::endl;
			refCompInfosCollapsed[row[refCompTab.getColPos("read")]] =
					row[refCompTab.getColPos("bestExpected")];
		}
	}

	if (bfs::exists(
			bib::files::join(sampleDir.string(), "refCompInfosExcluded.tab.txt"))) {
		table refCompTab(
				bib::files::join(sampleDir.string(), "refCompInfosExcluded.tab.txt"), "\t",
				true);
		//std::cout << __PRETTY_FUNCTION__ << " " << sampleName << std::endl;
		for (const auto & row : refCompTab.content_) {
			//std::cout << "\t" << row[refCompTab.getColPos("read")] << " :  " <<  row[refCompTab.getColPos("bestExpected")] << std::endl;
			refCompInfosExcluded[row[refCompTab.getColPos("read")]] =
					row[refCompTab.getColPos("bestExpected")];
		}
	}

	for (const auto & excludedClustersFile : excludedClustersFiles) {
		auto excludedClustersFileOpts = SeqIOOptions(excludedClustersFile.string(),
				inputOptions_.inFormat_, true);
		SeqInput reader(excludedClustersFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			allSampCounts[seq.getOwnSampName()] += seq.cnt_;
		}
	}

	for (const auto & sCount : allSampCounts) {
		allSampInfos[sCount.first] = sampInfo(sCount.first, sCount.second);
	}

	for (const auto & sCount : finalSampCounts) {
		finalSampInfos[sCount.first] = sampInfo(sCount.first, sCount.second);
	}
	auto finalFileOpts = SeqIOOptions(
			bib::files::make_path(sampleDir, "final/",
					sampleName + inputOptions_.getOutExtension()).string(),
			inputOptions_.inFormat_, true);

	auto excludedFileOpts = SeqIOOptions(
			bib::files::make_path(sampleDir, "excluded/",
					sampleName + inputOptions_.getOutExtension()).string(),
			inputOptions_.inFormat_, true);

	seqInfo subInfo;
	if(finalFileOpts.inExists()){
		SeqInput reader(finalFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			auto finalSubFileOpts = SeqIOOptions(
					bib::files::make_path(sampleDir, "final/clusters",
							seq.name_ + inputOptions_.getOutExtension()).string(),
					inputOptions_.inFormat_, true);
			SeqInput subReader(finalSubFileOpts);
			auto subReads = subReader.readAllReads<seqInfo>();
			samp.collapsed_.clusters_.emplace_back(
					sampleCluster(seq, subReads, finalSampInfos));
			samp.collapsed_.clusters_.back().expectsString =
					refCompInfosCollapsed[samp.collapsed_.clusters_.back().seqBase_.name_];
			for (const auto & subRead : subReads) {
				samp.input_.clusters_.emplace_back(
						sampleCluster(subRead, allSampInfos));
			}
		}
	}
	if(excludedFileOpts.inExists()){
		SeqInput reader(excludedFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			auto excludedSubFileOpts = SeqIOOptions(
					bib::files::make_path(sampleDir, "excluded/clusters",
							seq.name_ + inputOptions_.getOutExtension()).string(),
					inputOptions_.inFormat_, true);
			SeqInput subReader(excludedSubFileOpts);
			auto subReads = subReader.readAllReads<seqInfo>();
			samp.excluded_.clusters_.emplace_back(
					sampleCluster(seq, subReads, allSampInfos));
			samp.excluded_.clusters_.back().expectsString =
					refCompInfosExcluded[samp.excluded_.clusters_.back().seqBase_.name_];
			for (const auto & subRead : subReads) {
				samp.input_.clusters_.emplace_back(
						sampleCluster(subRead, allSampInfos));
			}
		}
	}
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

void SampleCollapseCollection::dumpSample(const std::string & sampleName) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	auto samp = sampleCollapses_.at(sampleName);
	auto sampDir = bib::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName);
	auto finalDir = bib::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName, "final");
	auto finalClustersDir = bib::files::make_path(masterOutputDir_,
			"samplesOutput", sampleName, "final", "clusters");
	auto excludedDir = bib::files::make_path(masterOutputDir_, "samplesOutput",
			sampleName, "excluded");
	auto excludedClusteredDir = bib::files::make_path(masterOutputDir_,
			"samplesOutput", sampleName, "excluded", "clusters");

	bib::files::makeDir(bib::files::MkdirPar(sampDir.string(), true));

	bib::files::makeDirP(bib::files::MkdirPar(finalClustersDir.string()));

	bib::files::makeDirP(bib::files::MkdirPar(excludedClusteredDir.string()));
	//final
	if(!samp->collapsed_.clusters_.empty()){
		SeqIOOptions finalOutOpts(
				bib::files::join(finalDir.string(), sampleName).string()
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);

		SeqOutput::write(samp->collapsed_.clusters_, finalOutOpts);
		for (const auto & clus : samp->collapsed_.clusters_) {
			SeqIOOptions clusOutOpts(
					bib::files::join(finalClustersDir.string(), clus.seqBase_.name_).string()
							+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
			clus.writeClusters(clusOutOpts);
		}
	}

	//excluded
	if(!samp->excluded_.clusters_.empty()){
		SeqIOOptions excluddOutOpts(
				bib::files::join(excludedDir.string(), sampleName).string()
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		SeqOutput::write(samp->excluded_.clusters_, excluddOutOpts);
		for (const auto & clus : samp->excluded_.clusters_) {
			SeqIOOptions clusOutOpts(
					bib::files::join(excludedClusteredDir.string(), clus.seqBase_.name_).string()
							+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
			clus.writeClusters(clusOutOpts);
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
								bib::files::join(sampDir.string(), "refCompInfos.tab.txt")),
						"\t", true));

		std::map<std::string, std::string> refCompInfosExcluded;
		for (const auto & clus : samp->excluded_.clusters_) {
			refCompInfosExcluded[clus.seqBase_.name_] = clus.expectsString;
		}
		table refCompExcludedTab(refCompInfosExcluded, VecStr { "read", "bestExpected" });
		refCompExcludedTab.outPutContents(
				TableIOOpts(
						OutOptions(
								bib::files::join(sampDir.string(), "refCompInfosExcluded.tab.txt")),
						"\t", true));
	}
	clearSample(sampleName);
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

void SampleCollapseCollection::investigateChimeras(double chiCutOff,
		aligner & alignerObj) {
	comparison onlyHpErrors;
	onlyHpErrors.oneBaseIndel_ = 1;
	onlyHpErrors.twoBaseIndel_ = 0;
	onlyHpErrors.largeBaseIndel_ = .99;

	table chiInfoTab(VecStr { "sample", "numClustersSaved",
			"totalClustersChecked", "clusterSavedNames" });
	for (const auto & sampleName : popNames_.samples_) {
		VecStr clustersSavedFromChi;
		setUpSampleFromPrevious(sampleName);
		uint32_t clustersNotSaved = 0;
		//first mark the suspicious clusters as being chimeric
		for (auto &clus : sampleCollapses_[sampleName]->collapsed_.clusters_) {
			if (clus.isClusterAtLeastChimericCutOff(chiCutOff)) {
				clus.seqBase_.markAsChimeric();
			}
		}
		//now check to see if it is ever the top two variants of a sample and if it is unmark it
		for (auto & clus : sampleCollapses_[sampleName]->collapsed_.clusters_) {
			if (clus.seqBase_.isChimeric()) {
				bool saved = false;
				for (const auto & otherSampleName : popNames_.samples_) {
					if (otherSampleName != sampleName
							&& !bib::containsSubString(stringToLowerReturn(otherSampleName),
									"control")
							&& !bib::containsSubString(stringToLowerReturn(sampleName),
									"control")) {
						auto otherSampOpts = SeqIOOptions(
								bib::files::make_path(masterOutputDir_, "samplesOutput",
										otherSampleName, "final",
										otherSampleName + inputOptions_.getOutExtension()).string(),
								inputOptions_.inFormat_, true);
						SeqInput otherReader(otherSampOpts);
						otherReader.openIn();
						uint32_t seqCount = 0;
						seqInfo otherSeq;
						while (otherReader.readNextRead(otherSeq) && seqCount < 2) {
							++seqCount;
							if (otherSeq.frac_ < 0.01) {
								continue;
							}
							alignerObj.alignCacheGlobal(otherSeq, clus.seqBase_);
							alignerObj.profilePrimerAlignment(otherSeq, clus.seqBase_);

							if (onlyHpErrors.passErrorProfile(alignerObj.comp_)) {
								clus.seqBase_.unmarkAsChimeric();
								for (auto & subRead : clus.reads_) {
									subRead->seqBase_.unmarkAsChimeric();
								}
								clus.resetInfos();
								clustersSavedFromChi.emplace_back(clus.seqBase_.name_);
								saved = true;
								break;
							}

						}
					}
				}
				if (!saved) {
					++clustersNotSaved;
				}
			}
			sampleCollapses_[sampleName]->collapsed_.setSetInfo();
		}
		chiInfoTab.content_.emplace_back(
				toVecStr(sampleName,
						getPercentageString(clustersSavedFromChi.size(),
								clustersSavedFromChi.size() + clustersNotSaved),
						clustersSavedFromChi.size() + clustersNotSaved,
						vectorToString(clustersSavedFromChi, ",")));
		dumpSample(sampleName);
	}
	TableIOOpts chiOutOptions(
			OutOptions(masterOutputDir_.string() + "chiInfo", ".tab.txt"), "\t",
			chiInfoTab.hasHeader_);
	chiInfoTab.outPutContents(chiOutOptions);
}

std::vector<sampleCluster> SampleCollapseCollection::createPopInput() {
	std::vector<sampleCluster> output;
	for (const auto & sampleName : popNames_.samples_) {
		setUpSampleFromPrevious(sampleName);
		auto & samp = *(sampleCollapses_[sampleName]);
		double totalReadCnt_ = 0;
		for (const auto &out : samp.collapsed_.clusters_) {
			totalReadCnt_ += out.seqBase_.cnt_;
		}
		std::map<std::string, sampInfo> outSampInfos { { sampleName, sampInfo(
				sampleName, totalReadCnt_) } };
		for (const auto &out : samp.collapsed_.clusters_) {
			output.emplace_back(sampleCluster(out.createRead(), outSampInfos));
			output.back().updateName();
			output.back().reads_.front()->seqBase_.name_ =
					output.back().seqBase_.name_;
		}
		clearSample(sampleName);
	}
	return output;
}

void SampleCollapseCollection::doPopulationClustering(
		const std::vector<sampleCluster> & input, aligner & alignerObj,
		const collapser & collapserObj, const CollapseIterations & popColIters) {
	popCollapse_ = std::make_unique<collapse::populationCollapse>(input,
			popNames_.populationName_);
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

void SampleCollapseCollection::dumpPopulation(const bfs::path & popDir, bool dumpTable) {
	if (popCollapse_) {
		auto popSubClusDir = bib::files::make_path(popDir, "clusters");
		//clear previous population dir or create if for the first time
		bib::files::makeDir(bib::files::MkdirPar(popDir.string(), true));
		bib::files::makeDir(bib::files::MkdirPar(popSubClusDir.string(), true));
		SeqIOOptions popOutOpts(
				bib::files::join(popDir.string(), "PopSeqs").string()
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
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
									bib::files::join(popDir.string(), "refCompInfos.tab.txt")),
							"\t", true));

		}
		std::unordered_map<std::string, double> readTotals;
		for (const auto & clus : popCollapse_->collapsed_.clusters_) {
			clus.writeClusters(
					SeqIOOptions(
							bib::files::join(popSubClusDir.string(), clus.seqBase_.name_).string()
									+ inputOptions_.getOutExtension(), inputOptions_.outFormat_));
			for (const auto & subClus : clus.reads_) {
				readTotals[subClus->getOwnSampName()] += subClus->seqBase_.cnt_;
			}
		}
		table readNumsTab(readTotals, VecStr { "sample", "readTotal" });
		readNumsTab.sortTable("sample", false);
		readNumsTab.outPutContents(
				TableIOOpts(
						OutOptions(bib::files::join(popDir.string(), "readTotals.tab.txt")),
						"\t", true));
		if(dumpTable){
			printPopulationCollapseInfo(bib::files::make_path(popDir, "populationCluster.tab.txt"));
		}
		popCollapse_ = nullptr;
	}
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
	auto popSubClusDir = bib::files::make_path(popDir, "clusters");
	popCollapse_ = std::make_unique<collapse::populationCollapse>(
			popNames_.populationName_);

	table readNumsTab(
			bib::files::make_path(popDir,
					"readTotals.tab.txt").string(), "\t", true);
	std::unordered_map<std::string, sampInfo> sampInfos;
	for (const auto & row : readNumsTab.content_) {
		auto sampName = row[readNumsTab.getColPos("sample")];
		sampInfos[sampName] = sampInfo(sampName,
				bib::lexical_cast<double>(row[readNumsTab.getColPos("readTotal")]));
	}

	MapStrStr refCompInfos;
	if (bfs::exists(bib::files::join(popDir.string(), "refCompInfos.tab.txt"))) {
		table refCompTab(bib::files::join(popDir.string(), "refCompInfos.tab.txt"),
				"\t", true);
		for (const auto & row : refCompTab.content_) {
			refCompInfos[row[refCompTab.getColPos("read")]] =
					row[refCompTab.getColPos("bestExpected")];
		}
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
		auto subClustersPath = bib::files::make_path(popSubClusDir,
				seq.name_ + inputOptions_.getOutExtension());
		SeqIOOptions subClusInOpts = SeqIOOptions(subClustersPath.string(),
				inputOptions_.inFormat_, true);
		SeqInput subReader(subClusInOpts);
		subReader.openIn();
		seqInfo subSeq;
		while (subReader.readNextRead(subSeq)) {
			//a little hacky but has to be done for certain pipeline;
			std::string inputSampleName = subSeq.getOwnSampName();
			if(MetaDataInName::nameHasMetaData(subSeq.name_)){
				MetaDataInName meta(subSeq.name_);
				if(meta.containsMeta("samp")){
					inputSampleName = meta.getMeta("samp");
				}
			}
			if (bib::in(inputSampleName, samples)) {
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
	popTab.outPutContents(TableIOOpts(OutOptions(fileName.string(), ".txt"), "\t", true));
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
	sampTab.outPutContents(TableIOOpts(OutOptions(fileName.string(), ".txt"), "\t", true));
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
		setUpSampleFromPrevious(sampName);
		if (maxRunCount
				< sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size()) {
			maxRunCount =
					sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size();
		}
		clearSample(sampName);
	}
	for (const auto& sampName : samples) {
		setUpSampleFromPrevious(sampName);
		for (const auto clusPos : iter::range(
				sampleCollapses_.at(sampName)->collapsed_.clusters_.size())) {
			const auto & clus =
					sampleCollapses_.at(sampName)->collapsed_.clusters_[clusPos];
			std::stringstream rowStream;
			rowStream  << sampName
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
							<< delim << clus.getClusterInfo(delim)
							<< delim << clus.getRepsInfo(
							sampleCollapses_.at(sampName)->input_.info_.infos_,
							sampleCollapses_.at(sampName)->excluded_.info_.infos_,
							sampleCollapses_.at(sampName)->collapsed_.info_.infos_,
							maxRunCount,
							checkingExpected, delim);
			rows.emplace_back(tokenizeString(rowStream.str(), delim, true));
		}
		clearSample(sampName);
	}
	std::stringstream headerStream;
	headerStream << "s_Sample"
			<< delim << sampleCluster::getPopInfoHeader(delim)
			<< delim << collapse::sampleCollapse::getSimpleSampInfoHeader(delim);
	if(nullptr != groupMetaData_){
		for(const auto & group : groupMetaData_->groupData_){
			headerStream << delim << "s_" << group.first;
		}
	}
	headerStream << delim << "c_clusterID"
			<< delim << sampleCluster::getClusterInfoHeader(delim) << delim
			<< sampleCluster::getRepsInfoHeader(maxRunCount, checkingExpected, delim);
	table ret(tokenizeString(headerStream.str(), delim));
	ret.addRows(rows);
	return ret;
}

void SampleCollapseCollection::symlinkInSampleFinals() const {
	auto finalDir = bib::files::make_path(masterOutputDir_, "final");
	bib::files::makeDir(bib::files::MkdirPar(finalDir.string(), true));
	for (const auto & samp : popNames_.samples_) {
		auto sampFinalHapFile = getSampleFinalHapsPath(samp);
		if(bfs::exists(sampFinalHapFile)){
			bfs::create_symlink(bfs::relative(sampFinalHapFile, finalDir), bib::files::make_path(finalDir, sampFinalHapFile.filename()));
		}
	}
}

bool SampleCollapseCollection::groupMetaLoaded() const{
	return nullptr != groupMetaData_;
}

void SampleCollapseCollection::createGroupInfoFiles(){
	if(nullptr != groupMetaData_){
		auto groupsTopDir = bib::files::make_path(masterOutputDir_, "groups");
		bib::files::makeDir(bib::files::MkdirPar(groupsTopDir.string(), true));
		std::ofstream groupMetaJsonFile;
		openTextFile(groupMetaJsonFile,
				bib::files::make_path(groupsTopDir, "groupMetaData.json").string(),
				".json", false, true);
		groupMetaJsonFile << groupMetaData_->toJson() << std::endl;
		for(const auto & group : groupMetaData_->groupData_){
			auto mainGroupDir = bib::files::make_path(groupsTopDir, group.first);
			bib::files::makeDir(bib::files::MkdirPar(mainGroupDir.string(), true));
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
							if(bib::in(uid, otherUid.second)){
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
				popTab.second.addColumn({bib::conToStr(uniquePopUids[popTab.first])}, "g_hapsFoundOnlyInThisGroup");
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
					"p_TotalHaplotypes", "p_TotalUniqueHaplotypes", "p_meanCoi", "p_medianCoi", "p_minCoi",
					"p_maxCoi", "g_hapsFoundOnlyInThisGroup"};
			//info on all the sub groups
			table outTab(groupInfoColNames);

			for(const auto & popTab : popTabs){
				auto subGroupDir = bib::files::make_path(mainGroupDir, popTab.first);
				bib::files::makeDir(bib::files::MkdirPar(subGroupDir.string()));
				popTab.second.outPutContents(TableIOOpts(OutOptions(bib::files::make_path(subGroupDir,"popFile.tab.txt")), "\t", true));
				sampTabs.at(popTab.first).outPutContents(TableIOOpts(OutOptions(bib::files::make_path(subGroupDir,"sampFile.tab.txt")), "\t", true));
				hapIdTabs.at(popTab.first).outPutContents(TableIOOpts(OutOptions(bib::files::make_path(subGroupDir,"hapIdTable.tab.txt")), "\t", true));
				std::ofstream subGroupMetaJsonFile;
				openTextFile(subGroupMetaJsonFile,
						bib::files::make_path(subGroupDir, "subGroupNamesData.json").string(),
						".json", false, true);
				auto popNames = popTab.second.getColumn("h_popUID");
				auto sampNames = sampTabs.at(popTab.first).getColumnLevels("s_Name");
				Json::Value nameMetaData;
				nameMetaData["popUIDs"] = bib::json::toJson(popNames);
				nameMetaData["sampNames"] = bib::json::toJson(sampNames);
				subGroupMetaJsonFile << nameMetaData << std::endl;
				outTab.rbind(popTab.second.getColumns(groupInfoColNames), false);
			}
			outTab = outTab.getUniqueRows();
			outTab.sortTable("g_GroupName", false);
			outTab.outPutContents(TableIOOpts(OutOptions(bib::files::make_path(mainGroupDir,"groupInfo.tab.txt")), "\t", true));
		}

		bool verbose = false;
		if(verbose){
			std::cout << "Missing meta data for the following samples:" << std::endl;
			std::cout << bib::conToStr(groupMetaData_->missingMetaForSamples_, "\n") << std::endl;
			std::cout << "Missing clustering data for the following samples:" << std::endl;
			std::cout << bib::conToStr(groupMetaData_->missingSamples_, "\n") << std::endl;
		}
	}
}

void SampleCollapseCollection::outputRepAgreementInfo() {
	bfs::path sampRepInfoDir = bib::files::makeDir(masterOutputDir_.string(),
			bib::files::MkdirPar("sampRepAgreementInfo", true));
	table rmseTab(VecStr { "sampleName", "RMSE" });
	table repInfoTab(
			VecStr { "sampleName", "clusterName", "repName", "fraction" });
	for (const auto & sampleName : popNames_.samples_) {
		setUpSampleFromPrevious(sampleName);
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
		clearSample(sampleName);
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
	OutOptions outOpts(bib::files::make_path(masterOutputDir_, "coreInfo.json"));
	outOpts.overWriteFile_ = true;
	openTextFile(outCoreJsonFile, outOpts);
	outCoreJsonFile << coreJson << std::endl;
}

Json::Value SampleCollapseCollection::toJsonCore() const{
	Json::Value ret;
	ret["inputOptions_"] = bib::json::toJson(inputOptions_);
	ret["masterInputDir_"] = bib::json::toJson(bfs::absolute(masterInputDir_));
	ret["masterOutputDir_"] = bib::json::toJson(bfs::absolute(masterOutputDir_));
	ret["samplesOutputDir_"] = bib::json::toJson(bfs::absolute(samplesOutputDir_));
	ret["populationOutputDir_"] = bib::json::toJson(bfs::absolute(populationOutputDir_));
	ret["popNames_"] = bib::json::toJson(popNames_);
	ret["clusterSizeCutOff_"] = bib::json::toJson(clusterSizeCutOff_);
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
	bib::sort(popUidNames);
	std::vector<std::vector<std::string>> inputContent = std::vector<std::vector<std::string>>{popUidNames.size(), std::vector<std::string>{}};
	for(const auto & popPos : iter::range(popUidNames.size())){
		inputContent[popPos].emplace_back(popUidNames[popPos]);
	}
	table ret(inputContent, VecStr{"#PopUID"});
	for (const auto & samp : samples) {
		setUpSampleFromPrevious(samp);
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
		clearSample(samp);
	}
	return ret;
}


}  // namespace collapse
}  // namespace bibseq

