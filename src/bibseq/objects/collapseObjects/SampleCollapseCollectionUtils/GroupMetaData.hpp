#pragma once
/*
 * GroupMetaData.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: nick
 */


#include "bibseq/utils.h"


namespace bibseq {
class GroupMetaData {
public:

	GroupMetaData(const std::string & name);

	std::string name_;
	std::set<std::string> subGroupsLevels_;
	std::unordered_map<std::string, std::set<std::string>> subGroupToSamples_;
	std::unordered_map<std::string, std::string> sampleToSubGroup_;

	void addSampGroup(const std::string & sample, const std::string & subGroup);
	void setSubGroupsLevels();

	std::string getGroupForSample(const std::string & samp) const;

};

class MultipleGroupMetaData {
public:
	MultipleGroupMetaData(const bfs::path & groupingsFile,
			const std::set<std::string> & availableSamples);
	const bfs::path groupingsFile_;

	std::set<std::string> samples_;
	std::set<std::string> missingSamples_;
	std::set<std::string> missingMetaForSamples_;
	std::map<std::string, std::unique_ptr<GroupMetaData>> groupData_;

};


}  // namespace bibseq



