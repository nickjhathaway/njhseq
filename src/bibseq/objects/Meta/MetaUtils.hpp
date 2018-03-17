#pragma once
/*
 * MetaUtils.hpp
 *
 *  Created on: Mar 14, 2018
 *      Author: nick
 */

#include "bibseq/objects/Meta/MetaDataInName.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"


namespace bibseq {


template<typename T>
table seqsToMetaTable(const std::vector<T> & seqs){
	std::set<std::string> allMetaKeys;
	for(auto & seq : seqs){
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			for(const auto & meta : metaData.meta_){
				allMetaKeys.emplace(meta.first);
			}
		}
	}
	VecStr colNames{"name", "count"};
	addOtherVec(colNames, VecStr(allMetaKeys.begin(), allMetaKeys.end()));
	colNames.emplace_back("seq");
	table outTab(colNames);

	for (const auto & seq : seqs) {
		VecStr row;
		std::unique_ptr<MetaDataInName> metaData = nullptr;
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			metaData = std::make_unique<MetaDataInName>(getSeqBase(seq).name_);
		}
		for (const auto & col : outTab.columnNames_) {
			if ("name" == col) {
				row.emplace_back(getSeqBase(seq).name_);
			} else if ("seq" == col) {
				row.emplace_back(getSeqBase(seq).seq_);
			} else if ("count" == col) {
				row.emplace_back(estd::to_string(getSeqBase(seq).cnt_));
			} else {
				if (nullptr != metaData && metaData->containsMeta(col)) {
					row.emplace_back(metaData->getMeta(col));
				} else {
					row.emplace_back("NA");
				}
			}
		}
		outTab.content_.emplace_back(row);
	}
	return outTab;
}

template<typename T>
table seqsToMetaTable(const std::vector<T> & seqs, const std::set<std::string> & fields){
//	VecStr colNames{"name", "count"};
//	addOtherVec(colNames, VecStr(fields.begin(), fields.end()));
//	colNames.emplace_back("seq");
	VecStr colNames(fields.begin(), fields.end());
	table outTab(colNames);
	for (const auto & seq : seqs) {
		VecStr row;
		std::unique_ptr<MetaDataInName> metaData = nullptr;
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			metaData = std::make_unique<MetaDataInName>(getSeqBase(seq).name_);
		}
		for (const auto & col : outTab.columnNames_) {
			if ("name" == col) {
				row.emplace_back(getSeqBase(seq).name_);
			} else if ("seq" == col) {
				row.emplace_back(getSeqBase(seq).seq_);
			} else if ("count" == col) {
				row.emplace_back(estd::to_string(getSeqBase(seq).cnt_));
			} else {
				if (nullptr != metaData && metaData->containsMeta(col)) {
					row.emplace_back(metaData->getMeta(col));
				} else {
					row.emplace_back("NA");
				}
			}
		}
		outTab.content_.emplace_back(row);
	}
	return outTab;
}



}  // namespace bibseq



