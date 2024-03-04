#pragma once
/*
 * MultiGenomeMapper.hpp
 *
 *  Created on: Mar 31, 2017
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
#include <njhseq/objects/Gene/GeneFromGffs.hpp>

#include "njhseq/utils.h"
#include "njhseq/GenomeUtils/GenomeMapping/GenomicRegionCounter.hpp"
#include "njhseq/objects/BioDataObject/reading.hpp"
#include "njhseq/system.h"
#include "njhseq/common/stdAddition.hpp"


namespace njhseq {

class MultiGenomeMapper {
public:
	struct IntersectedProteinInfo {
		IntersectedProteinInfo(std::string id, uint32_t aaStart, uint32_t aaStop,
		                  std::string desc) : id_(std::move(id)), aaStart_(aaStart), aaStop_(aaStop),
		                                      description_(std::move(desc)) {
		}

		std::string id_; //!< transcript ID
		uint32_t aaStart_; //!< 1-based
		uint32_t aaStop_; //!< 1-based
		std::string description_; //!< detailed description of the gene of the transcript

		MetaDataInName allMeta_; //!< any other associated meta data for the protein
	};

	struct inputParameters{
		inputParameters();
		inputParameters(const bfs::path & genomeDir, const std::string & primaryGenome);
		bfs::path genomeDir_;
		std::string primaryGenome_;
		std::set<std::string> selectedGenomes_;
		uint32_t numThreads_ = 1;

		bfs::path gffDir_;

		njh::files::MkdirPar workingDirectory_{""};

		bool keepTempFiles_ = false;
		bool verbose_ = false;

		intersectBedLocsWtihGffRecordsPars gffIntersectPars_;

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
		bfs::path gffFnp_;
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
	void loadGffFnps();
	void loadGffFnps(const bfs::path & gffDir);
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
		njh::sys::RunOutput rOutput_;
		Json::Value toJson() const{
			Json::Value ret;
			ret["class"] = njh::getTypeName(*this);
			ret["alignedFnp_"] = njh::json::toJson(alignedFnp_);
			ret["rOutput_"] = njh::json::toJson(rOutput_);
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


	template <typename BEDREC>
	std::unordered_map<std::string, std::vector<IntersectedProteinInfo>> addIntersectingGeneInfosToLocs(
		std::vector<BEDREC>& regions,
		const std::string&genome,
		const std::string&gffExtraAttributesStr = "description");

	struct getRefSeqsWithPrimaryGenomePars{
		BioCmdsUtils::LastZPars lzPars;
		bool keepBestOnly = false;
		bool byScore = false;
		bool extendAndTrim = false;
		uint32_t extendAndTrimLen = 10;
		bool shortNames = false;
	};

	std::vector<seqInfo> getRefSeqsWithPrimaryGenome(const GenomicRegion & region,
			const bfs::path & alignmentsDir,
			const getRefSeqsWithPrimaryGenomePars & pars,
			aligner & alignerObj) const;



	std::unordered_map<std::string, std::vector<seqInfo>> getRefSeqsWithPrimaryGenomeAll(
			const GenomicRegion & region, const bfs::path & alignmentsDir,
			const BioCmdsUtils::LastZPars & lzPars) const;

	static std::unordered_map<std::string, bfs::path> getBamFnps(const std::unordered_map<std::string, MultiGenomeMapper::AlignCmdOutput> & alignOutpus);
};


template<typename BEDREC>
std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>>
MultiGenomeMapper::addIntersectingGeneInfosToLocs(
	std::vector<BEDREC>&regions,
	const std::string&genome,
	const std::string&gffExtraAttributesStr) {
	checkForGenomeThrow(genome, __PRETTY_FUNCTION__);
	std::unordered_map<std::string, std::vector<IntersectedProteinInfo>> ret;


	intersectBedLocsWtihGffRecordsPars pars(genomes_.at(genome)->gffFnp_);
	pars.selectFeatures_ = VecStr{"gene", "protein_coding_gene"};
	pars.extraAttributes_ = tokenizeString(gffExtraAttributesStr, ",");
	auto jsonValues = intersectBedLocsWtihGffRecords(regions, pars);
	auto rawGeneIds = jsonValues.getMemberNames();
	std::set<std::string> rawReneIdsSet(rawGeneIds.begin(), rawGeneIds.end());

	auto rawGenes = GeneFromGffs::getGenesFromGffForIds(genomes_.at(genome)->gffFnp_, rawReneIdsSet);
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;

	for (const auto&gene: rawGenes) {
		bool failFilter = false;
		for (const auto&transcript: gene.second->mRNAs_) {
			if (njh::in(njh::strToLowerRet(transcript->type_), VecStr{"rrna", "trna", "snorna", "snrna", "ncrna"})) {
				failFilter = true;
				break;
			}
		}
		if (!failFilter) {
			genes[gene.first] = gene.second;
		}
	}
	for (const auto&reg: regions) {
		if (getRef(reg).extraFields_.empty()) {
			continue;
		}
		auto geneToks = njh::tokenizeString(getRef(reg).extraFields_.back(), "],[");
		if (geneToks.size() > 1) {
			geneToks.front().append("]");
			for (const auto pos: iter::range<uint32_t>(1, geneToks.size() - 1)) {
				geneToks[pos] = "[" + geneToks[pos] + "]";
			}
			geneToks.back() = "[" + geneToks.back();
		}
		std::string replacement;
		for (const auto&geneTok: geneToks) {
			MetaDataInName geneMeta(geneTok);
			if (njh::in(geneMeta.getMeta("ID"), genes)) {
				TwoBit::TwoBitFile tReader(genomes_.at(genome)->fnpTwoBit_);
				auto infos = njh::mapAt(genes, geneMeta.getMeta("ID"))->generateGeneSeqInfo(tReader, false);
				auto detailedName = njh::mapAt(genes, geneMeta.getMeta("ID"))->getGeneDetailedName();
				for (const auto&info: infos) {
					auto posInfos = info.second->getInfosByGDNAPos();

					auto minPos = vectorMinimum(getVectorOfMapKeys(posInfos));
					auto maxPos = vectorMaximum(getVectorOfMapKeys(posInfos));
					auto startPos = std::max(minPos, getRef(reg).chromStart_);
					auto stopPos = std::min(maxPos, getRef(reg).chromEnd_);
					if (stopPos == getRef(reg).chromEnd_) {
						stopPos = getRef(reg).chromEnd_ - 1;
					}
					//stop position is inclusive
					uint32_t aaStartPos = std::numeric_limits<uint32_t>::max();
					uint32_t aaStopPos = 0;

					for (const auto posInGene: iter::range(startPos, stopPos + 1)) {
						if (std::numeric_limits<uint32_t>::max() != posInfos[posInGene].aaPos_) {
							if (posInfos[posInGene].aaPos_ < aaStartPos) {
								aaStartPos = posInfos[posInGene].aaPos_;
							}
							if (posInfos[posInGene].aaPos_ > aaStopPos) {
								aaStopPos = posInfos[posInGene].aaPos_;
							}
						}
					}
					//if the region is completely within the intron the final stop and stops will be max() and 0;
					if (std::numeric_limits<uint32_t>::max() == aaStartPos) {
						aaStopPos = std::numeric_limits<uint32_t>::max();
					}
					else {
						//make 1 based
						++aaStopPos;
						++aaStartPos;
					}
					geneMeta.addMeta(info.first + "-AAStart", aaStartPos);
					geneMeta.addMeta(info.first + "-AAStop", aaStopPos);
					std::string description;
					bool allNAs = true;
					for (const auto&field: pars.extraAttributes_) {
						if (!description.empty()) {
							description += "::";
						}
						if ("NA" != geneMeta.getMeta(field)) {
							allNAs = false;
						}
						description += geneMeta.getMeta(field);
					}
					if (allNAs) {
						description = detailedName[info.first];
						geneMeta.addMeta("detailedDescription", detailedName[info.first]);
					}
					ret[getRef(reg).name_].emplace_back(geneMeta.getMeta("ID"), aaStartPos, aaStopPos, description);
				}
			}
			if (!replacement.empty()) {
				replacement += ",";
			}
			replacement += geneMeta.createMetaName();
		}
		getRef(reg).extraFields_.back() = replacement;
	}
	return ret;
}
}  // namespace njhseq
