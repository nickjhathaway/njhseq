/*
 * AllGroupDataPaths.cpp
 *
 *  Created on: Sep 13, 2016
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//


#include "AllGroupDataPaths.hpp"

namespace bibseq {


AllGroupDataPaths::GroupDataPaths::SubGroupDataPaths::SubGroupDataPaths(
		const bfs::path & mainDir) :
		mainDir_(mainDir), popFileFnp_(
				bib::files::make_path(mainDir, "popFile.tab.txt")), sampFileFnp_(
				bib::files::make_path(mainDir, "sampFile.tab.txt")),
				hapIdTabFnp_(bib::files::make_path(mainDir, "hapIdTable.tab.txt")),
				subGroupNamesDataFnp_(
				bib::files::make_path(mainDir, "subGroupNamesData.json")) {
}



Json::Value AllGroupDataPaths::GroupDataPaths::SubGroupDataPaths::toJson() const{
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["mainDir_"] = bib::json::toJson(mainDir_);
	ret["popFileFnp_"] = bib::json::toJson(popFileFnp_);
	ret["sampFileFnp_"] = bib::json::toJson(sampFileFnp_);
	ret["subGroupNamesDataFnp_"] = bib::json::toJson(subGroupNamesDataFnp_);
	return ret;
}

VecStr AllGroupDataPaths::GroupDataPaths::SubGroupDataPaths::readInPopUIDs() const {
	auto meta = bib::json::parseFile(subGroupNamesDataFnp_.string());
	return bib::json::jsonArrayToVec<std::string>(meta["popUIDs"],
			[](const Json::Value & val) {return val.asString();});
}

VecStr AllGroupDataPaths::GroupDataPaths::SubGroupDataPaths::readInSampNames() const {
	auto meta = bib::json::parseFile(subGroupNamesDataFnp_.string());
	return bib::json::jsonArrayToVec<std::string>(meta["sampNames"],
			[](const Json::Value & val) {return val.asString();});
}


AllGroupDataPaths::GroupDataPaths::GroupDataPaths(const bfs::path & mainDir,
		const std::set<std::string> & subGroups) :
		mainDir_(mainDir), groupInfoFnp_(
				bib::files::make_path(mainDir, "groupInfo.tab.txt")) {
	for (const auto & group : subGroups) {
		groupPaths_.emplace(group, bib::files::make_path(mainDir_, group));
	}
}

Json::Value AllGroupDataPaths::GroupDataPaths::toJson() const{
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["mainDir_"] = bib::json::toJson(mainDir_);
	ret["groupInfoFnp_"] = bib::json::toJson(groupInfoFnp_);
	ret["groupPaths_"] = bib::json::toJson(groupPaths_);

	return ret;
}

AllGroupDataPaths::AllGroupDataPaths(const bfs::path & mainDir,
		const std::unique_ptr<MultipleGroupMetaData> & groupMetaData) :
		mainDir_(mainDir) {
	if (nullptr == groupMetaData) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " error, groupMetaData is nullptr";
		throw std::runtime_error { ss.str() };
	}

	for (const auto & group : groupMetaData->groupData_) {
		allGroupPaths_.emplace(group.first,
				GroupDataPaths(bib::files::make_path(mainDir, group.first),
						group.second->subGroupsLevels_));
	}
}

Json::Value AllGroupDataPaths::toJson() const{
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["mainDir_"] = bib::json::toJson(mainDir_);
	ret["allGroupPaths_"] = bib::json::toJson(allGroupPaths_);

	return ret;
}




}  // namespace bibseq
