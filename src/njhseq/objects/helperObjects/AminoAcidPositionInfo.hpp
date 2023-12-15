#pragma once
//
// Created by Nicholas Hathaway on 12/15/23.
//



#include "njhseq/utils.h"
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/objects/Meta/MetaDataInName.hpp"


namespace njhseq {

class AminoAcidPositionInfo {
private:
	void intialize();
public:
	AminoAcidPositionInfo(table  infoTab, bool zeroBased);
	AminoAcidPositionInfo(const bfs::path&inputInfoFnp, bool zeroBased);

	[[nodiscard]] bool byRange() const;

	const bfs::path inputInfoFnp_;
	table infoTab_;
	bool zeroBased_ = false;
	std::set<std::string> ids_;
	std::map<std::string, std::set<uint32_t>> aminoPositionsPerId_;
	std::map<std::string, std::unordered_map<uint32_t, MetaDataInName>> metaDataForAAPos_;

};//



} //namespace njhseq

