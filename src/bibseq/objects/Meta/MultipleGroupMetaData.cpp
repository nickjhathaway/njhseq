/*
 * MultipleGroupMetaData.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: nick
 */


#include "MultipleGroupMetaData.hpp"

namespace bibseq {

void MultipleGroupMetaData::setInfoWithTable(const table & groupsTab,
		const std::set<std::string> & availableSamples){
	resetInfo();
	if (!groupsTab.containsColumn("Sample")
			&& !groupsTab.containsColumn("sample")
			&& !groupsTab.containsColumn("Samples")
			&& !groupsTab.containsColumn("samples") ) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, input table"
		   << " must contain the column \"Sample\""
			 << " or \"sample\"  or \"Samples\"  or \"samples\""  << "\n";
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
	}else if(groupsTab.containsColumn("sample")){
		samples = groupsTab.getColumn("sample");
	}else if(groupsTab.containsColumn("Samples")){
		samples = groupsTab.getColumn("Samples");
	}else if(groupsTab.containsColumn("samples")){
		samples = groupsTab.getColumn("samples");
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


void MultipleGroupMetaData::setInfoWithTable(const table & groupsTab){
	resetInfo();
	if (!groupsTab.containsColumn("Sample")
			&& !groupsTab.containsColumn("sample")
			&& !groupsTab.containsColumn("Samples")
			&& !groupsTab.containsColumn("samples") ) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, input table"
		   << " must contain the column \"Sample\""
			 << " or \"sample\"  or \"Samples\"  or \"samples\""  << "\n";
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
	}else if(groupsTab.containsColumn("sample")){
		samples = groupsTab.getColumn("sample");
	}else if(groupsTab.containsColumn("Samples")){
		samples = groupsTab.getColumn("Samples");
	}else if(groupsTab.containsColumn("samples")){
		samples = groupsTab.getColumn("samples");
	}

	for (const auto & samp : samples) {
		if (bib::in(samp, samples_) || bib::in(samp, missingSamples_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, have sample: " << samp
					<< " entered more than once in Sample column in " << groupingsFile_
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		samples_.insert(samp);
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
	setInfoWithTable(table(groupingsFile, "\t", true));
}


MultipleGroupMetaData::MultipleGroupMetaData(const table & info, const std::set<std::string> & availableSamples) : groupingsFile_("") {
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



void MultipleGroupMetaData::checkForFieldsThrow(const VecStr & fields) const {
	bool fail = false;
	std::stringstream errorStream;
	for(const auto & field : fields){
		if (!hasMetaField(field) ) {
			errorStream << "Error, need column named " << field << "\n";
			fail = true;
		}
	}
	if (fail) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", errors" << "\n";
		ss << errorStream.str();
		throw std::runtime_error { ss.str() };
	}
}

void MultipleGroupMetaData::transformSubFields(const std::string & field,
		std::function<std::string(const std::string & subField)> transformer) {
	if(!hasMetaField(field)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, don't have field " << field << "\n";
		ss << "Options are : " << bib::conToStr(bib::getVecOfMapKeys(groupData_)) << "\n";
		throw std::runtime_error { ss.str() };
	}
	auto subGroupToSamplesOld = groupData_.at(field)->subGroupToSamples_;
	groupData_.at(field)->subGroupToSamples_.clear();
	groupData_.at(field)->sampleToSubGroup_.clear();

	for(const auto & subFields : subGroupToSamplesOld){
		for(const auto & samp : subFields.second){
			groupData_.at(field)->addSampGroup(samp, transformer(subFields.first));
		}
	}

	groupData_.at(field)->setSubGroupsLevels();
}

bool MultipleGroupMetaData::hasMetaField(const std::string & metaField) const{
	return bib::in(metaField, groupData_);
}

bool MultipleGroupMetaData::hasSample(const std::string & sample) const{
	return bib::in(sample, samples_);
}

MetaDataInName MultipleGroupMetaData::getMetaForSample(const std::string & name, const VecStr & fields) const{
	if(!bib::in(name, samples_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no smaple: " << name << "\n";
		ss << "Options are: " << bib::conToStr(samples_);
		throw std::runtime_error{ss.str()};
	}
	checkForFieldsThrow(fields);

	MetaDataInName ret;
	for(const auto & gMeta : groupData_){
		if(bib::in(gMeta.first, fields)){
			ret.addMeta(gMeta.first, gMeta.second->getGroupForSample(name));
		}
	}
	return ret;
}

MetaDataInName MultipleGroupMetaData::genNaMeta( const VecStr & fields) const{
	checkForFieldsThrow(fields);

	MetaDataInName ret;
	for(const auto & gMeta : groupData_){
		if(bib::in(gMeta.first, fields)){
			ret.addMeta(gMeta.first, "NA");
		}
	}
	return ret;
}




}  // namespace bibseq
