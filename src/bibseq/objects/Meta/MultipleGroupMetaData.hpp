#pragma once
/*
 * MultipleGroupMetaData.hpp
 *
 *  Created on: Apr 7, 2017
 *      Author: nick
 */

#include "bibseq/objects/Meta/GroupMetaData.hpp"
#include "bibseq/objects/Meta/MetaDataInName.hpp"

namespace bibseq {


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

	/**@brief Set the group and sample data wtih input table, will clear any info currently held
	 *
	 * @param info a table containing the info, need to have more than 1 column, and one column must be Samples or samples
	 * @param availableSamples a set with the anmes of the samples that were actually clustered
	 */
	void setInfoWithTable(const table & info, const std::set<std::string> & availableSamples);

	/**@brief clear the info in samples_, missingSamples_, missingMetaFroSamples_, and groupData_
	 *
	 */
	void resetInfo();

	/**@brief convert to json
	 *
	 * @return
	 */
	Json::Value toJson() const;


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
	MetaDataInName genNaMeta(const VecStr & fields) const;

};

}  // namespace bibseq


