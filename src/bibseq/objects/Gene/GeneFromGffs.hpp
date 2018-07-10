#pragma once

/*
 * GeneFromGffs.hpp
 *
 *  Created on: Aug 7, 2017
 *      Author: nick
 */

#include "bibseq/objects/Gene/GeneSeqInfo.hpp"
#include <TwoBit.h>

namespace bibseq {

class GeneFromGffs {
public:
	GeneFromGffs(const std::vector<std::shared_ptr<GFFCore>> & geneRecords);

	std::shared_ptr<GFFCore> gene_ = nullptr;

	std::vector<std::shared_ptr<GFFCore>> mRNAs_;
	std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> CDS_;
	std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> exons_;
	std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> polypeptides_;

	std::unordered_map<std::string, std::vector<std::shared_ptr<GFFCore>>> others_;

	std::unordered_map<std::string, std::string> getGeneDetailedName() const;

	std::string getOneGeneDetailedName() const;


	uint32_t numOfTranscripts() const;

	std::unordered_map<std::string, table> getIntronExonTables() const;

	std::unordered_map<std::string, std::vector<Bed6RecordCore>> getIntronExonBedLocs() const;

	std::unordered_map<std::string,std::shared_ptr<GeneSeqInfo>> generateGeneSeqInfo(TwoBit::TwoBitFile & tReader,
			bool oneBased) const;

	void writeOutGeneInfo(TwoBit::TwoBitFile & tReader, const OutOptions & outPrefix) const;
	void writeGffRecords(std::ostream & out) const;

	static std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> getGenesFromGffForIds(const bfs::path & gffFnp, const std::set<std::string> & ids);

};


} /* namespace bibseq */

