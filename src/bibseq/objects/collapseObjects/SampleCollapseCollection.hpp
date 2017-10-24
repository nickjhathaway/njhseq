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

	class RepFile {
	public:
		RepFile(const std::string & repName, const bfs::path & repFnp) :
				repName_(repName), repFnp_(repFnp) {
		}
		std::string repName_;
		bfs::path repFnp_;
		bool reNameInput_ = true;
	};

	SampleCollapseCollection(SeqIOOptions inputOptions,
			const bfs::path & inputDir,
			const bfs::path & outputDir,
			const PopNamesInfo & popNames,
			uint32_t clusterSizeCutOff,
			uint32_t sampleMinReadCount);

	SampleCollapseCollection(const Json::Value & coreJson);

	SeqIOOptions inputOptions_;
	bfs::path masterInputDir_;
	bfs::path masterOutputDir_;

private:
	bfs::path samplesOutputDir_;
	bfs::path populationOutputDir_;
	std::mutex mut_;
public:
	PopNamesInfo popNames_{"", VecStr{}};
	VecStr passingSamples_;
	uint32_t clusterSizeCutOff_;
	uint32_t sampleMinReadCount_;
	std::map<std::string, std::shared_ptr<collapse::sampleCollapse>> sampleCollapses_;
	std::unique_ptr<populationCollapse> popCollapse_;
	std::unique_ptr<MultipleGroupMetaData> groupMetaData_;
	std::unique_ptr<AllGroupDataPaths> groupDataPaths_;



	void addGroupMetaData(const bfs::path & groupingsFile);
	void createGroupInfoFiles();
	bool groupMetaLoaded() const;

	void setUpSample(const std::string & sampleName,
			const std::vector<RepFile> & analysisFiles,
			aligner & alignerObj,
			const collapser & collapserObj, const ChimeraOpts & chiOpts);

	void setUpSample(const std::string & sampleName,
			aligner & alignerObj,
			const collapser & collapserObj,
			const ChimeraOpts & chiOpts);

	void setUpSampleFromPrevious(const std::string & sampleName);

	void clusterSample(const std::string & sampleName, aligner & alignerObj,
			const collapser & collapserObj, const CollapseIterations & colIters);

	void dumpSample(const std::string & sampleName);

	void clearSample(const std::string & sampleName);
	void clearSamples(const VecStr & sampleNames);

	bfs::path getSampleFinalHapsPath(const std::string & sampleName) const;
	bfs::path getPopFinalHapsPath() const;

	bfs::path getPopInfoPath() const;
	bfs::path getSampInfoPath() const;
	bfs::path getHapIdTabPath() const;

	uint32_t numOfSamples() const;

	bool hasSample(const std::string & sampleName) const;

	void checkForSampleThrow(const std::string & funcName, const std::string & sampleName) const;

	void investigateChimeras(double chiCutOff, aligner & alignerObj);

	std::vector<sampleCluster> createPopInput();

	void doPopulationClustering(const std::vector<sampleCluster> & input,
			aligner & alignerObj, const collapser & collapserObj,
			const CollapseIterations & popColIters);

	void dumpPopulation(bool dumpTable = true);
	void dumpPopulation(const bfs::path& outputPopDir, bool dumpTable = true);

	void loadInPreviousPop();
	void loadInPreviousPop(const std::set<std::string> & samples);
	void loadInPreviousPop(const std::set<std::string> & samples, const bfs::path& outputPopDir);

	void renamePopWithSeqs(const std::vector<readObject> & otherPopSeqs, comparison allowableErrors = comparison());
	void addRefMetaToName(const std::vector<readObject> & otherPopSeqs, comparison allowableErrors = comparison());

	void comparePopToRefSeqs(const std::vector<readObject> & expectedSeqs,
			aligner & alignerObj);


	void checkForPopCollapseThrow(const std::string & funcName) const;

	void printPopulationCollapseInfo(const bfs::path& fileName) const;
	table genPopulationCollapseInfo() const;

	void printSampleCollapseInfo(const bfs::path& fileName);
	table genSampleCollapseInfo(const std::set<std::string> & samples);


	void symlinkInSampleFinals() const;

	table genHapIdTable();
	table genHapIdTable(const std::set<std::string> & samples);

	void outputRepAgreementInfo();

	double estimateNumberOfStrains(const std::set<std::string> & samples);

	double calculatePIE(const std::set<std::string> & samples);

	Json::Value toJsonCore() const;

	void createCoreJsonFile() const;


};

}  // namespace collapse
}  // namespace bibseq

