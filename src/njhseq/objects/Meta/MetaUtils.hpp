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



class MetaFieldSeqFilterer {
public:
	struct MetaFieldSeqFiltererPars {
		std::string filteringMetaField_;
		std::string groupingMetaField_;
		std::string tieBreakerPreferenceField_;
		bool keepCommonSeqsWhenFiltering_ = false;
	};

	MetaFieldSeqFilterer(const MetaFieldSeqFiltererPars & pars) :
			pars_(pars) {
	}

	struct MetaFieldSeqFiltererRes {
		std::vector<uint32_t> filteredAllSeqs;
		VecStr failedGroups;
	};

	template<typename T>
	MetaFieldSeqFiltererRes filterSeqs(const std::vector<T> & allSeqs){
		//key = K1: ExperimentSample, K2: sample, V2: vector of seq positions
		// second map is ordered map to make results the same every time rather than determined by random order of unordered map
		std::unordered_map<std::string, std::map<std::string, std::vector<uint32_t>>> seqsByFilteringField;
		for(const auto & seqPos : iter::range(allSeqs.size())){
			const auto & seq = getSeqBase(allSeqs[seqPos]);
			MetaDataInName seqMeta(seq.name_);
			seqsByFilteringField[seqMeta.getMeta(pars_.filteringMetaField_)][seqMeta.getMeta(pars_.groupingMetaField_)].emplace_back(seqPos);
		}
		return filterSeqs(allSeqs, seqsByFilteringField);
	}

	template<typename T>
	MetaFieldSeqFiltererRes filterSeqs(const std::vector<T> & allSeqs, const std::vector<uint32_t> & selPositions){
		//key = K1: ExperimentSample, K2: sample, V2: vector of seq positions
		// second map is ordered map to make results the same every time rather than determined by random order of unordered map
		std::unordered_map<std::string, std::map<std::string, std::vector<uint32_t>>> seqsByFilteringField;
		for(const auto & seqPos : selPositions){
			if(seqPos > allSeqs.size()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " seqPos: " << seqPos << " out of range of allSeqs, size: " << allSeqs.size() << "\n";
				throw std::runtime_error{ss.str()};
			}
			const auto & seq = getSeqBase(allSeqs[seqPos]);
			MetaDataInName seqMeta(seq.name_);
			seqsByFilteringField[seqMeta.getMeta(pars_.filteringMetaField_)][seqMeta.getMeta(pars_.groupingMetaField_)].emplace_back(seqPos);
		}
		return filterSeqs(allSeqs, seqsByFilteringField);
	}


	//key = K1: ExperimentSample, K2: sample, V2: vector of seq positions
	// second map is ordered map to make results the same every time rather than determined by random order of unordered map
	template<typename T>
	MetaFieldSeqFiltererRes filterSeqs(const std::vector<T> & allSeqs, const std::unordered_map<std::string, std::map<std::string, std::vector<uint32_t>>> &seqsByFilteringField){
		MetaFieldSeqFiltererRes ret;
		for(const auto & seqsForFilterField : seqsByFilteringField){
			if(1 == seqsForFilterField.second.size()){
				addOtherVec(ret.filteredAllSeqs, seqsForFilterField.second.begin()->second);
			}else{
				//key K: sequence, V: count
				std::unordered_map<std::string, uint32_t> countsPerSeqs;
				for(const auto & seqsForSample : seqsForFilterField.second){
					for(const auto & seqPos : seqsForSample.second){
						++countsPerSeqs[getSeqBase(allSeqs[seqPos]).seq_];
					}
				}
				bool anyFailed = false;
				for(const auto & countsPerSeq : countsPerSeqs){
					if(countsPerSeq.second != seqsForFilterField.second.size()){
						anyFailed = true;
					}
				}
				if(anyFailed){
					ret.failedGroups.emplace_back(seqsForFilterField.first);
					if(pars_.keepCommonSeqsWhenFiltering_){
						for(const auto & seqPosForExp : seqsForFilterField.second.begin()->second){
							if(countsPerSeqs[getSeqBase(allSeqs[seqPosForExp]).seq_] == seqsForFilterField.second.size()){
								ret.filteredAllSeqs.emplace_back(seqPosForExp);
							}
						}
					} else {
						if("" != pars_.tieBreakerPreferenceField_) {
							for(const auto & subSample : seqsForFilterField.second) {
								MetaDataInName frontMeta(getSeqBase(allSeqs[subSample.second.front()]).name_);
								if(frontMeta.containsMeta(pars_.tieBreakerPreferenceField_) && frontMeta.getMeta<bool>("PreferredSample") ) {
									addOtherVec(ret.filteredAllSeqs, subSample.second);
									break;
								}
							}
						}
					}
				}else{
					bool added = false;
					if("" != pars_.tieBreakerPreferenceField_) {
						for(const auto & subSample : seqsForFilterField.second) {
							MetaDataInName frontMeta(getSeqBase(allSeqs[subSample.second.front()]).name_);
							if(frontMeta.containsMeta(pars_.tieBreakerPreferenceField_) && frontMeta.getMeta<bool>("PreferredSample") ) {
								addOtherVec(ret.filteredAllSeqs, subSample.second);
								added = true;
								break;
							}
						}
					}
					if(!added){
						addOtherVec(ret.filteredAllSeqs, seqsForFilterField.second.begin()->second);
					}
				}
			}
		}
		return ret;
	}


	MetaFieldSeqFiltererPars pars_;
};




}  // namespace njhseq



