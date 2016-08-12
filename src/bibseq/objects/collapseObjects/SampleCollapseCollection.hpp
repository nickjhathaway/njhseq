#pragma once
/*
 * SampleCollapseCollection.hpp
 *
 *  Created on: Jul 24, 2016
 *      Author: nick
 */

#include "bibseq/objects/collapseObjects/sampleCollapse.hpp"
#include "bibseq/objects/collapseObjects/populationCollapse.hpp"
#include "bibseq/objects/collapseObjects/SampleCollapseCollectionUtils.h"



namespace bibseq {
namespace collapse {



class SampleCollapseCollection {

public:
	SampleCollapseCollection(SeqIOOptions inputOptions,
			const bfs::path & inputDir,
			const bfs::path & outputDir,
			const PopNamesInfo & popNames,
			uint32_t clusterSizeCutOff);

	const SeqIOOptions inputOptions_;
	const bfs::path masterInputDir_;
	const bfs::path masterOutputDir_;
private:
	bfs::path samplesOutputDir_;
	bfs::path populationOutputDir_;
public:
	PopNamesInfo popNames_;
	uint32_t clusterSizeCutOff_;
	std::map<std::string, std::shared_ptr<collapse::sampleCollapse>> sampleCollapses_;
	std::shared_ptr<populationCollapse> popCollapse_;
	std::unique_ptr<MultipleGroupMetaData> groupMetaData_;

	void addGroupMetaData(const bfs::path & groupingsFile);

	void setUpSample(const std::string & sampleName, aligner & alignerObj,
			const collapser & collapserObj, const ChimeraOpts & chiOpts);

	void setUpSampleFromPrevious(const std::string & sampleName);

	void clusterSample(const std::string & sampleName, aligner & alignerObj,
			const collapser & collapserObj, const CollapseIterations & colIters);

	void dumpSample(const std::string & sampleName);

	void clearSample(const std::string & sampleName);
	void clearSamples(const VecStr & sampleNames);

	bfs::path getSampleFinalHapsPath(const std::string & sampleName) const;
	bfs::path getPopFinalHapsPath() const;

	uint32_t numOfSamples() const;

	bool hasSample(const std::string & sampleName) const;

	void checkForSampleThrow(const std::string & funcName, const std::string & sampleName) const;

	void investigateChimeras(double chiCutOff, aligner & alignerObj);

	std::vector<sampleCluster> createPopInput();

	void doPopulationClustering(const std::vector<sampleCluster> & input,
			aligner & alignerObj, const collapser & collapserObj,
			const CollapseIterations & popColIters);

	void dumpPopulation();
	void dumpPopulation(const bfs::path& outputPopDir);

	void loadInPreviousPop();
	void loadInPreviousPop(const std::set<std::string> & samples);
	void loadInPreviousPop(const std::set<std::string> & samples, const bfs::path& outputPopDir);

	void renamePopWithSeqs(const std::vector<readObject> & otherPopSeqs);

	void comparePopToRefSeqs(const std::vector<readObject> & expectedSeqs,
			aligner & alignerObj);


	void checkForPopCollapseThrow(const std::string & funcName) const;

	void printPopulationCollapseInfo(const bfs::path& fileName) const;
	table genPopulationCollapseInfo() const;

	void printSampleCollapseInfo(const bfs::path& fileName);
	table genSampleCollapseInfo(const std::set<std::string> & samples);


	void symlinkInSampleFinals() const;

	void createGroupInfoFiles();

	void outputRepAgreementInfo();

	double estimateNumberOfStrains(const std::set<std::string> & samples);

	double calculatePIE(const std::set<std::string> & samples);

};

}  // namespace collapse
}  // namespace bibseq

