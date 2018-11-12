/*
 * GroupMetaData.cpp
 *
 *  Created on: Jul 30, 2016
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
#include "GroupMetaData.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"

namespace njhseq {

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
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["name_"] = njh::json::toJson(name_);
	ret["subGroupsLevels_"] = njh::json::toJson(subGroupsLevels_);
	ret["subGroupToSamples_"] = njh::json::toJson(subGroupToSamples_);
	return ret;
}


VecStr GroupMetaData::getSampleNames() const {
	return njh::getVecOfMapKeys(sampleToSubGroup_);
}

GroupMetaData GroupMetaData::fromJson(const Json::Value & jsonValue){
	njh::json::MemberChecker checker(jsonValue);
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


}  // namespace njhseq


