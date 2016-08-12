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
		const bfs::path & inputDir, const bfs::path & outputDir,
		const PopNamesInfo & popNames, uint32_t clusterSizeCutOff) :
		inputOptions_(inputOptions), masterInputDir_(inputDir), masterOutputDir_(
				outputDir), popNames_(popNames), clusterSizeCutOff_(clusterSizeCutOff) {
	samplesOutputDir_ = bib::files::make_path(masterOutputDir_, "samplesOutput");
	populationOutputDir_ = bib::files::make_path(masterOutputDir_, "population");
	bib::files::makeDir(bib::files::MkdirPar(samplesOutputDir_.string(), true));

}



void SampleCollapseCollection::addGroupMetaData(
		const bfs::path & groupingsFile) {
	groupMetaData_ = std::make_unique<MultipleGroupMetaData>(groupingsFile,
			popNames_.samples_);

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
		ss << funcName << ": error, sample isn't in samples_" << "\n";
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
			std::regex { "^" + inputOptions_.firstName_ + "$" } }, 2);
	std::vector<std::vector<cluster>> inputClusters;
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
		auto sampleOpts = inputOptions_;
		sampleOpts.firstName_ = af.first.string();
		SeqInput reader(sampleOpts);

		std::vector<cluster> clusters = baseCluster::convertVectorToClusterVector<
				cluster>(reader.readAllReads<seqInfo>());
		readVecSorter::sortReadVector(clusters, "totalCount");
		// consider adding the sample name in the name as well
		renameReadNamesNewClusters(clusters, fileToks[0], true, true, false);
		if (chiOpts.checkChimeras_) {
			collapserObj.markChimeras(clusters, alignerObj, chiOpts);
		}
		clusterVec::allSetFractionClusters(clusters);
		readVec::allUpdateName(clusters);
		inputClusters.emplace_back(clusters);
	}
	if (bib::in(sampleName, sampleCollapses_)) {
		sampleCollapses_[sampleName] = std::make_shared<sampleCollapse>(
				inputClusters, sampleName, clusterSizeCutOff_);
	} else {
		sampleCollapses_.emplace(sampleName,
				std::make_shared<sampleCollapse>(inputClusters, sampleName,
						clusterSizeCutOff_));
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
		SeqInput reader(finalClustersFileOpts);
		reader.openIn();
		while (reader.readNextRead(seq)) {
			allSampCounts[seq.getOwnSampName()] += seq.cnt_;
			finalSampCounts[seq.getOwnSampName()] += seq.cnt_;
		}
	}

	MapStrStr refCompInfos;
	if (bfs::exists(
			bib::files::join(sampleDir.string(), "refCompInfos.tab.txt"))) {
		table refCompTab(
				bib::files::join(sampleDir.string(), "refCompInfos.tab.txt"), "\t",
				true);
		for (const auto & row : refCompTab.content_) {
			refCompInfos[row[refCompTab.getColPos("read")]] =
					row[refCompTab.getColPos("bestExpected")];
		}
	}

	for (const auto & excludedClustersFile : excludedClustersFiles) {
		auto finalClustersFileOpts = SeqIOOptions(excludedClustersFile.string(),
				inputOptions_.inFormat_, true);
		SeqInput reader(finalClustersFileOpts);
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
	{
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
					refCompInfos[samp.collapsed_.clusters_.back().seqBase_.name_];
			for (const auto & subRead : subReads) {
				samp.input_.clusters_.emplace_back(
						sampleCluster(subRead, allSampInfos));
			}
		}
	}
	{
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
					refCompInfos[samp.excluded_.clusters_.back().seqBase_.name_];
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
	auto & samp = *(sampleCollapses_[sampleName]);
	samp.cluster(collapserObj, colIters, sortBy, alignerObj);

	samp.updateCollapsedInfos();
	samp.updateExclusionInfos();
	samp.renameClusters(sortBy);
}

void SampleCollapseCollection::dumpSample(const std::string & sampleName) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	auto & samp = *(sampleCollapses_[sampleName]);
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
	SeqIOOptions finalOutOpts(
			bib::files::join(finalDir.string(), sampleName)
					+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
	SeqOutput::write(samp.collapsed_.clusters_, finalOutOpts);
	for (const auto & clus : samp.collapsed_.clusters_) {
		SeqIOOptions clusOutOpts(
				bib::files::join(finalClustersDir.string(), clus.seqBase_.name_)
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		clus.writeClusters(clusOutOpts);
	}
	SeqIOOptions excluddOutOpts(
			bib::files::join(excludedDir.string(), sampleName)
					+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
	SeqOutput::write(samp.excluded_.clusters_, excluddOutOpts);
	for (const auto & clus : samp.excluded_.clusters_) {
		SeqIOOptions clusOutOpts(
				bib::files::join(excludedClusteredDir.string(), clus.seqBase_.name_)
						+ inputOptions_.getOutExtension(), inputOptions_.outFormat_);
		clus.writeClusters(clusOutOpts);
	}

	if (!samp.collapsed_.clusters_.empty()
			&& "" != samp.collapsed_.clusters_.front().expectsString) {

		std::map<std::string, std::string> refCompInfos;
		for (const auto & clus : samp.collapsed_.clusters_) {
			refCompInfos[clus.seqBase_.name_] = clus.expectsString;
		}
		for (const auto & clus : samp.excluded_.clusters_) {
			refCompInfos[clus.seqBase_.name_] = clus.expectsString;
		}
		table refCompTab(refCompInfos, VecStr { "read", "bestExpected" });
		refCompTab.outPutContents(
				TableIOOpts(
						OutOptions(
								bib::files::join(sampDir.string(), "refCompInfos.tab.txt")),
						"\t", true));
	}
	clearSample(sampleName);
}

void SampleCollapseCollection::clearSample(const std::string & sampleName) {
	checkForSampleThrow(__PRETTY_FUNCTION__, sampleName);
	sampleCollapses_[sampleName] = nullptr;
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
	popCollapse_ = std::make_shared<collapse::populationCollapse>(input,
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

void SampleCollapseCollection::dumpPopulation(const bfs::path & popDir) {
	if (popCollapse_) {
		auto popSubClusDir = bib::files::make_path(popDir, "clusters");
		//clear previous population dir or create if for the first time
		bib::files::makeDir(bib::files::MkdirPar(popDir.string(), true));
		bib::files::makeDir(bib::files::MkdirPar(popSubClusDir.string(), true));
		SeqIOOptions popOutOpts(
				bib::files::join(popDir.string(), "PopSeqs")
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
							bib::files::join(popSubClusDir.string(), clus.seqBase_.name_)
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
		printPopulationCollapseInfo(bib::files::make_path(popDir, "populationCluster.tab.txt"));
		popCollapse_ = nullptr;
	}
}

void SampleCollapseCollection::dumpPopulation() {
	dumpPopulation(populationOutputDir_);
}

void SampleCollapseCollection::loadInPreviousPop() {
	loadInPreviousPop(popNames_.samples_, populationOutputDir_);
}

void SampleCollapseCollection::loadInPreviousPop(const std::set<std::string> & samples){
	loadInPreviousPop(samples, populationOutputDir_);
}

void SampleCollapseCollection::loadInPreviousPop(const std::set<std::string> & samples, const bfs::path & popDir) {
	auto popSubClusDir = bib::files::make_path(popDir, "clusters");
	popCollapse_ = std::make_shared<collapse::populationCollapse>(
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
			if (bib::in(subSeq.getOwnSampName(), samples)) {
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
		const std::vector<readObject> & otherPopSeqs) {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	popCollapse_->renameToOtherPopNames(otherPopSeqs);
}

void SampleCollapseCollection::comparePopToRefSeqs(
		const std::vector<readObject> & expectedSeqs, aligner & alignerObj) {
	checkForPopCollapseThrow(__PRETTY_FUNCTION__);
	popCollapse_->collapsed_.checkAgainstExpected(expectedSeqs, alignerObj, false,
			false);
}


void SampleCollapseCollection::checkForPopCollapseThrow(
		const std::string & funcName) const {
	if (nullptr == popCollapse_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error popCollapse_ not loaded\n";
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
	if (nullptr != groupMetaData_) {
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
		if(nullptr != groupMetaData_){
			for(const auto & group : groupMetaData_->groupData_){
				std::map<std::string, uint32_t> subGroupCounts;;
				for(const auto & samp : clus.sampInfos()){;
					if(samp.second.numberOfClusters_ > 0){
						++subGroupCounts[group.second->getGroupForSample(samp.first)];
					}
				}
				std::string groupCountsStr;
				std::string groupFracsStr;
				for(const auto & subGroupCount : subGroupCounts){
					groupCountsStr += bib::pasteAsStr(subGroupCount.first, ":", subGroupCount.second, ";");
					groupFracsStr +=  bib::pasteAsStr(subGroupCount.first, ":", subGroupCount.second/static_cast<double>(clus.numberOfRuns()), ";");
				}
				addOtherVec(row, toVecStr(groupCountsStr, groupFracsStr));
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
		auto currentInfos = sampleCollapses_.at(sampName)->getAllInfoMap(
				checkingExpected, delim,
				sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size());
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
							sampleCollapses_.at(sampName)->collapsed_.numOfReps(),
							checkingExpected, delim);
			rows.emplace_back(tokenizeString(rowStream.str(), delim, true));
		}
		if (maxRunCount
				< sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size()) {
			maxRunCount =
					sampleCollapses_.at(sampName)->collapsed_.info_.infos_.size();
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

void SampleCollapseCollection::createGroupInfoFiles(){
	if(nullptr != groupMetaData_){
		auto groupsTopDir = bib::files::make_path(masterOutputDir_, "groups");
		bib::files::makeDir(bib::files::MkdirPar(groupsTopDir.string(), true));
		for(const auto & group : groupMetaData_->groupData_){
			auto mainGroupDir = bib::files::make_path(groupsTopDir, group.first);
			bib::files::makeDir(bib::files::MkdirPar(mainGroupDir.string(), true));
			std::unordered_map<std::string, table> popTabs;
			std::unordered_map<std::string, table> sampTabs;
			for(const auto & subGroup : group.second->subGroupToSamples_){
				loadInPreviousPop(subGroup.second);
				popTabs[subGroup.first] = genPopulationCollapseInfo();
				sampTabs[subGroup.first] = genSampleCollapseInfo(subGroup.second);
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
			for(const auto & popTab : popTabs){
				auto subGroupDir = bib::files::make_path(mainGroupDir, popTab.first);
				bib::files::makeDir(bib::files::MkdirPar(subGroupDir.string()));
				popTab.second.outPutContents(TableIOOpts(OutOptions(bib::files::make_path(subGroupDir,"popFile.tab.txt").string()), "\t", true));
				sampTabs.at(popTab.first).outPutContents(TableIOOpts(OutOptions(bib::files::make_path(subGroupDir,"sampFile.tab.txt").string()), "\t", true));
			}
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
	std::string sampRepInfoDir = bib::files::makeDir(masterOutputDir_.string(),
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
					OutOptions(sampRepInfoDir + "replicatesFractions", ".tab.txt"), "\t",
					repInfoTab.hasHeader_));
	rmseTab.sortTable("RMSE", true);
	rmseTab.outPutContents(
			TableIOOpts(OutOptions(sampRepInfoDir + "RMSE.tab.txt", ".tab.txt"), "\t",
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

}  // namespace collapse
}  // namespace bibseq





























