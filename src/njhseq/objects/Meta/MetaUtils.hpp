#pragma once
/*
 * MetaUtils.hpp
 *
 *  Created on: Mar 14, 2018
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
#include "njhseq/objects/Meta/MetaDataInName.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"


namespace njhseq {


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



}  // namespace njhseq



