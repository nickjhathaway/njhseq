/*
 * GroupMetaData.cpp
 *
 *  Created on: Jul 30, 2016
 *      Author: nick
 */

#include "GroupMetaData.hpp"
#include "bibseq/objects/dataContainers/tables/table.hpp"

namespace bibseq {

GroupMetaData::GroupMetaData(const std::string & name) :
		name_(name) {
}

void GroupMetaData::addSampGroup(const std::string & sample,
		const std::string & subGroup) {
	std::string insertingGroup = subGroup;
	if("" == subGroup){
		insertingGroup = "NA";
	}
	subGroupToSamples_[insertingGroup].insert(sample);
	sampleToSubGroup_[sample] = insertingGroup;
}

void GroupMetaData::setSubGroupsLevels() {
	auto groupLevels = getVectorOfMapKeys(subGroupToSamples_);
	subGroupsLevels_ = std::set<std::string>{groupLevels.begin(), groupLevels.end()};
}

std::string GroupMetaData::getGroupForSample(const std::string & samp) const {
	std::string ret = "NA";
	auto search = sampleToSubGroup_.find(samp);
	if (sampleToSubGroup_.end() != search) {
		ret = search->second;
	}
	return ret;
}

Json::Value GroupMetaData::toJson() const {
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["name_"] = bib::json::toJson(name_);
	ret["subGroupsLevels_"] = bib::json::toJson(subGroupsLevels_);
	ret["subGroupToSamples_"] = bib::json::toJson(subGroupToSamples_);
	return ret;
}


VecStr GroupMetaData::getSampleNames() const {
	return bib::getVecOfMapKeys(sampleToSubGroup_);
}

GroupMetaData GroupMetaData::fromJson(const Json::Value & jsonValue){
	bib::json::MemberChecker checker(jsonValue);
	checker.failMemberCheckThrow(VecStr { "name_", "subGroupToSamples_"},
			__PRETTY_FUNCTION__);

	GroupMetaData ret(jsonValue["name_"].asString());

	auto subGroupNames = jsonValue["subGroupToSamples_"].getMemberNames();
	uint32_t groupPos = 0;
	for(const auto & group : jsonValue["subGroupToSamples_"]){
		for(const auto & samp : group){
			ret.addSampGroup(samp.asString(), subGroupNames[groupPos]);
		}
		++groupPos;
	}

	ret.setSubGroupsLevels();

	return ret;
}


}  // namespace bibseq


