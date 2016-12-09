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
	subGroupToSamples_[subGroup].insert(sample);
	sampleToSubGroup_[sample] = subGroup;
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


void MultipleGroupMetaData::setInfoWithTable(const table & groupsTab,
		const std::set<std::string> & availableSamples){
	resetInfo();
	if (!groupsTab.containsColumn("Sample")
			&& !groupsTab.containsColumn("sample")) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, input table"
		   << " must contain the column \"Sample\""
			 << " or \"sample\"" << "\n";
		ss << "Current columns are: " << vectorToString(groupsTab.columnNames_, ",")
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	if (groupsTab.nCol() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, input table"
		   << " must contain at least two columns, current size is : "
				<< groupsTab.nCol() << "\n";
		throw std::runtime_error { ss.str() };
	}

	VecStr samples;
	if(groupsTab.containsColumn("Sample")){
		samples = groupsTab.getColumn("Sample");
	}else{
		samples = groupsTab.getColumn("sample");
	}

	for (const auto & samp : samples) {
		if (bib::in(samp, samples_) || bib::in(samp, missingSamples_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, have sample: " << samp
					<< " entered more than once in Sample column in " << groupingsFile_
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		if (bib::in(samp, availableSamples)) {
			samples_.insert(samp);
		} else {
			missingSamples_.insert(samp);
		}
	}
	for (const auto & samp : availableSamples) {
		if (!bib::in(samp, samples)) {
			missingMetaForSamples_.insert(samp);
		}
	}
	for (const auto & col : groupsTab.columnNames_) {
		if ("Sample" != col && "sample" != col) {
			if (bib::in(col, groupData_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << " Error, have grouping: " << col
						<< " entered more than once in the column names in "
						<< groupingsFile_ << "\n";
				throw std::runtime_error { ss.str() };
			}
			groupData_.emplace(col, std::make_unique<GroupMetaData>(col));
			VecStr grouping = groupsTab.getColumn(col);
			for (const auto pos : iter::range(grouping.size())) {
				if (bib::in(samples[pos], samples_)) {
					groupData_[col]->addSampGroup(samples[pos], grouping[pos]);
				}
			}
			for(const auto & samp : missingMetaForSamples_){
				groupData_[col]->addSampGroup(samp, "NA");
			}
			groupData_[col]->setSubGroupsLevels();
		}
	}
}

MultipleGroupMetaData::MultipleGroupMetaData(const bfs::path & groupingsFile,
		const std::set<std::string> & availableSamples) :
		groupingsFile_(groupingsFile) {
	table groupsTab(groupingsFile_.string(), "\t", true);
	setInfoWithTable(groupsTab, availableSamples);
}


MultipleGroupMetaData::MultipleGroupMetaData(const bfs::path & groupingsFile) :
		groupingsFile_(groupingsFile) {

}


MultipleGroupMetaData::MultipleGroupMetaData(const table & info, const std::set<std::string> & availableSamples):groupingsFile_("") {
	setInfoWithTable(info, availableSamples);
}

void MultipleGroupMetaData::resetInfo(){
	samples_.clear();
	missingSamples_.clear();
	missingMetaForSamples_.clear();
	groupData_.clear();
}

Json::Value MultipleGroupMetaData::toJson() const {
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));
	ret["groupingsFile_"] = bib::json::toJson(groupingsFile_);
	ret["samples_"] = bib::json::toJson(samples_);
	ret["missingSamples_"] = bib::json::toJson(missingSamples_);
	ret["missingMetaForSamples_"] = bib::json::toJson(missingMetaForSamples_);
	ret["groupData_"] = bib::json::toJson(groupData_);
	return ret;
}

MultipleGroupMetaData MultipleGroupMetaData::fromJson(
		const Json::Value & jsonValue) {
	bib::json::MemberChecker checker(jsonValue);
	checker.failMemberCheckThrow(VecStr { "groupingsFile_", "samples_",
			"missingSamples_", "missingMetaForSamples_", "groupData_" },
			__PRETTY_FUNCTION__);

	MultipleGroupMetaData ret(jsonValue["groupingsFile_"].asString());
	auto jsonToStr = [](const Json::Value & v) {return v.asString();};
	ret.samples_ = bib::json::jsonArrayToSet<std::string>(jsonValue["samples_"],
			jsonToStr);
	ret.missingSamples_ = bib::json::jsonArrayToSet<std::string>(
			jsonValue["missingSamples_"], jsonToStr);
	ret.missingMetaForSamples_ = bib::json::jsonArrayToSet<std::string>(
			jsonValue["missingMetaForSamples_"], jsonToStr);
	auto groupNames = jsonValue["groupData_"].getMemberNames();
	uint32_t groupCount = 0;
	for (const auto & group : jsonValue["groupData_"]) {
		ret.groupData_.emplace(groupNames[groupCount],
				std::make_unique<GroupMetaData>(GroupMetaData::fromJson(group)));
		++groupCount;
	}
	return ret;
}

MultipleGroupMetaData::GroupPopInfo::GroupPopInfo(const std::string & groupName,
		uint32_t numOfSamples) :
		groupName_(groupName), numOfSamples_(numOfSamples) {

}

void MultipleGroupMetaData::GroupPopInfo::increaseSubGroupCount(
		const std::string & subGroupName) {
	++subGroupCounts_[subGroupName];
}

std::string MultipleGroupMetaData::GroupPopInfo::groupCountsStr() const{
	std::string ret { "" };
	for (const auto & subGroupCount : subGroupCounts_) {
		ret += bib::pasteAsStr(subGroupCount.first, ":", subGroupCount.second, ";");
	}
	return ret;
}

std::string MultipleGroupMetaData::GroupPopInfo::groupFracsStr() const{
	std::string ret { "" };
	for (const auto & subGroupCount : subGroupCounts_) {
		ret += bib::pasteAsStr(subGroupCount.first, ":",
				subGroupCount.second / static_cast<double>(numOfSamples_), ";");
	}
	return ret;
}

std::vector<MultipleGroupMetaData::GroupPopInfo> MultipleGroupMetaData::getGroupPopInfos(const VecStr & samples){
	std::vector<GroupPopInfo> ret;
	for(const auto & group : groupData_){
		GroupPopInfo info(group.first, samples.size());
		for (const auto & samp : samples) {
			info.increaseSubGroupCount(group.second->getGroupForSample(samp));
		}
		ret.emplace_back(info);
	}
	return ret;
}


}  // namespace bibseq


