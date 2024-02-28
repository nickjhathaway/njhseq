#pragma once
/*
 * MultipleGroupMetaData.hpp
 *
 *  Created on: Apr 7, 2017
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
#include "njhseq/objects/Meta/GroupMetaData.hpp"
#include "njhseq/objects/Meta/MetaDataInName.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

namespace njhseq {


/**@brief a class hold multiple groups of meta data for samples
 *
 */
class MultipleGroupMetaData {
public:
	/**@brief construct with a path to the meta data file and the clustered samples
	 *
	 * @param groupingsFile the meta data file path name
	 * @param availableSamples the samples that were clustered
	 */
	MultipleGroupMetaData(const bfs::path & groupingsFile,
			const std::set<std::string> & availableSamples);

	/**@brief construct with just the groupings file locations, no info loaded
	 *
	 * @param groupingsFile the path to a file containing group meta data;
	 */
	MultipleGroupMetaData(const bfs::path & groupingsFile);

	/**@brief construct with a table with the group info, groupingsFile_ will just be blank
	 *
	 * @param info table with group info
	 */
	MultipleGroupMetaData(const table & info,
			const std::set<std::string> & availableSamples);



	const bfs::path groupingsFile_; /**< file path to the group meta data */
	std::set<std::string> samples_; /**< samples that are both in the group meta group and that were clustered */
	std::set<std::string> missingSamples_; /**< samples that are in the group meta group but that woeren't clustered */
	std::set<std::string> missingMetaForSamples_; /**< samples that were clustered but weren't the group meta data */
	std::map<std::string, std::unique_ptr<GroupMetaData>> groupData_; /**< the multiple group data */

	/**@brief Set the group and sample data with input table, will clear any info currently held
	 *
	 * @param info a table containing the info, need to have more than 1 column, and one column must be sample, Sample, samples, or Samples
	 * @param availableSamples a set with the names of the samples that were actually clustered, will only consider these sample when reading
	 */
	void setInfoWithTable(const table & info, const std::set<std::string> & availableSamples);

	/**@brief Set the group and sample data with input table
	 *
	 * @param info the table to set info with, must have a column, sample, Sample, samples, or Samples
	 */
	void setInfoWithTable(const table & info);


	/**@brief clear the info in samples_, missingSamples_, missingMetaFroSamples_, and groupData_
	 *
	 */
	void resetInfo();

	/**@brief convert to json
	 *
	 * @return
	 */
	Json::Value toJson() const;

	void writeOutMetaFile(const OutOptions & outOpts) const;
	void writeOutMetaFile(const OutOptions & outOpts, const std::set<std::string> & samples) const;

	/**@brief create a MultipleGroupMetaData object from json
	 *
	 * @param jsonValue a json value representing MultipleGroupMetaData, likely created from MultipleGroupMetaData::toJson()
	 * @return a MultipleGroupMetaData object
	 */
	static MultipleGroupMetaData fromJson(const Json::Value & jsonValue);

	class GroupPopInfo {
	public:
		GroupPopInfo(const std::string & groupName, uint32_t numOfSamples);

		std::string groupName_;
		uint32_t numOfSamples_;
		std::map<std::string, uint32_t> subGroupCounts_;

		void increaseSubGroupCount(const std::string & subGroupName);
		std::string groupCountsStr() const;
		std::string groupFracsStr() const;

	};

	void checkForFieldsThrow(const VecStr & fields)const;

	void transformSubFields(const std::string & field, std::function<std::string(const std::string & subField)> transformer);

	std::vector<GroupPopInfo> getGroupPopInfos(const VecStr & samples);


	bool hasMetaField(const std::string & metaField) const;
	bool hasSample(const std::string & sample) const;

	MetaDataInName getMetaForSample(const std::string & name, const VecStr & fields) const;
	MetaDataInName getMetaForSample(const std::string & name) const;

	MetaDataInName genNaMeta(const VecStr & fields) const;


	table leftJoinWithMeta(const table & sampleTable,
			const std::string & sampleColumn) const;


	table leftJoinWithMeta(const table & sampleTable) const;

	template<typename T>
	bool attemptToAddSeqMeta(T & seq) const{
		if(MetaDataInName::nameHasMetaData(getSeqBase(seq).name_)){
			MetaDataInName metaData(getSeqBase(seq).name_);
			if(metaData.containsMeta("sample") && njh::in(metaData.getMeta("sample"), samples_)){
				auto seqMeta = getMetaForSample(metaData.getMeta("sample"), getVectorOfMapKeys(groupData_));
				metaData.addMeta(seqMeta, true);
				metaData.resetMetaInName(getSeqBase(seq).name_);
				return true;
			}else if(hasSample(getSeqBase(seq).name_)){
				auto seqMeta = getMetaForSample(getSeqBase(seq).name_, getVectorOfMapKeys(groupData_));
				seqMeta.resetMetaInName(getSeqBase(seq).name_);
				return true;
			}
		}else{
			if(hasSample(getSeqBase(seq).name_)){
				auto seqMeta = getMetaForSample(getSeqBase(seq).name_, getVectorOfMapKeys(groupData_));
				seqMeta.resetMetaInName(getSeqBase(seq).name_);
				return true;
			}
		}
		return false;
	}
};

}  // namespace njhseq


