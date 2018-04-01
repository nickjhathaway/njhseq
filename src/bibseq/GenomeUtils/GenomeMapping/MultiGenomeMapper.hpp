#pragma once
/*
 * MultiGenomeMapper.hpp
 *
 *  Created on: Mar 31, 2017
 *      Author: nick
 */

#include "bibseq/utils.h"
#include "bibseq/GenomeUtils/GenomeMapping/GenomicRegionCounter.hpp"
#include "bibseq/system.h"

namespace bibseq {

class MultiGenomeMapper {
public:

	struct inputParameters{
		bfs::path genomeDir_;
		std::string primaryGenome_;
		std::set<std::string> selectedGenomes_;
		uint32_t numThreads_ = 1;


		bib::files::MkdirPar workingDirectory_{""};

		bool keepTempFiles_ = false;
		bool verbose_ = false;

		Json::Value toJson() const;
	};

	MultiGenomeMapper(const inputParameters & pars);
	MultiGenomeMapper(const bfs::path & genomeDir,
			const std::string & primaryGenome);

	inputParameters pars_;

	Json::Value toJson() const;

	class Genome {
	public:
		Genome(const std::string & name, const bfs::path & fnp);
		std::string name_;
		bfs::path fnp_;
		bfs::path fnpTwoBit_;
		std::unordered_map<std::string, uint32_t> chromosomeLengths_;

		void createTwoBit();

		void buildBowtie2Index()const;

		Json::Value chromosomeLengths() const;

		Json::Value toJson() const;
	};

	std::unordered_map<std::string, std::unique_ptr<Genome>> genomes_;
	std::mutex mut_;

	void checkForGenomeThrow(const std::string & genome, const std::string & funcName) const;
	bool hasGenome(const std::string & genome) const;

	void loadInGenomes();
	void setUpGenomes();

	void bioIndexAllGenomes();

	void init();

	void setSelectedGenomes(const std::set<std::string> & genomes);
	void setSelectedGenomes(const VecStr & genomes);
	void setSelectedGenomes(const std::string & gnomesStr);

	std::vector<bfs::path> getGenomeFnps() const;

	struct AlignCmdOutput {
		AlignCmdOutput(const bfs::path & alignedFnp) :
				alignedFnp_(alignedFnp) {
		}
		bfs::path alignedFnp_;
		bib::sys::RunOutput rOutput_;
		Json::Value toJson() const{
			Json::Value ret;
			ret["class"] = bib::getTypeName(*this);
			ret["alignedFnp_"] = bib::json::toJson(alignedFnp_);
			ret["rOutput_"] = bib::json::toJson(rOutput_);
			return ret;
		}

	};

	std::unordered_map<std::string, AlignCmdOutput> alignToGenomes(
			const SeqIOOptions & inputOpts, const bfs::path & outputPrefix) const;

	std::unordered_map<std::string, AlignCmdOutput> alignToGenomesLastz(
			const SeqIOOptions & inputOpts, const bfs::path & outputPrefix,
			const BioCmdsUtils::LastZPars & pars = BioCmdsUtils::LastZPars()) const;

	std::unordered_map<std::string, std::vector<GenomicRegion>> getRegionsFromBams(
			const std::unordered_map<std::string, bfs::path> & bamFnps) const;

	std::vector<seqInfo> extractRegions(
			const std::unordered_map<std::string, GenomicRegion> & regions) const;

	seqInfo extractGenomeRegion(const std::string & genome, const GenomicRegion & region) const;


	std::vector<seqInfo> getRefSeqsWithPrimaryGenome(const GenomicRegion & region,
			const bfs::path & alignmentsDir,
			const BioCmdsUtils::LastZPars & lzPars,
			bool keepBestOnly) const;


	std::vector<seqInfo> getRefSeqsWithPrimaryGenome(const GenomicRegion & region,
			const bfs::path & alignmentsDir,
			const BioCmdsUtils::LastZPars & lzPars,
			bool keepBestOnly,
			bool extendAndTrim,
			uint32_t extendAndTrimLen,
			aligner & alignerObj) const;

	std::unordered_map<std::string, std::vector<seqInfo>> getRefSeqsWithPrimaryGenomeAll(
			const GenomicRegion & region, const bfs::path & alignmentsDir,
			const BioCmdsUtils::LastZPars & lzPars) const;


	static std::unordered_map<std::string, bfs::path> getBamFnps(const std::unordered_map<std::string, MultiGenomeMapper::AlignCmdOutput> & alignOutpus);

};

}  // namespace bibseq
