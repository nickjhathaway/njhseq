/*
 * MultipleGroupMetaData.cpp
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

#include "MultipleGroupMetaData.hpp"

#include <njhseq/IO/OutputStream.hpp>

namespace njhseq {

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
		if (njh::in(samp, samples_) || njh::in(samp, missingSamples_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, have sample: " << samp
					<< " entered more than once in Sample column in " << groupingsFile_
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		if (njh::in(samp, availableSamples)) {
			samples_.insert(samp);
		} else {
			missingSamples_.insert(samp);
		}
	}
	for (const auto & samp : availableSamples) {
		if (!njh::in(samp, samples)) {
			missingMetaForSamples_.insert(samp);
		}
	}
	for (const auto & col : groupsTab.columnNames_) {
		if ("Sample" != col && "sample" != col) {
			if (njh::in(col, groupData_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << " Error, have grouping: " << col
						<< " entered more than once in the column names in "
						<< groupingsFile_ << "\n";
				throw std::runtime_error { ss.str() };
			}
			groupData_.emplace(col, std::make_unique<GroupMetaData>(col));
			VecStr grouping = groupsTab.getColumn(col);
			for (const auto pos : iter::range(grouping.size())) {
				if (njh::in(samples[pos], samples_)) {
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
	std::string sampleNameColName = "Sample";
	if(groupsTab.containsColumn("Sample")){
		samples = groupsTab.getColumn("Sample");
		sampleNameColName = "Sample";
	}else if(groupsTab.containsColumn("sample")){
		samples = groupsTab.getColumn("sample");
		sampleNameColName = "sample";
	}else if(groupsTab.containsColumn("Samples")){
		samples = groupsTab.getColumn("Samples");
		sampleNameColName = "Samples";
	}else if(groupsTab.containsColumn("samples")){
		samples = groupsTab.getColumn("samples");
		sampleNameColName = "samples";
	}

	for (const auto & samp : samples) {
		if (njh::in(samp, samples_) || njh::in(samp, missingSamples_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " Error, have sample: " << samp
					<< " entered more than once in Sample column in " << groupingsFile_
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
		samples_.insert(samp);
	}
	for (const auto & col : groupsTab.columnNames_) {
		if (sampleNameColName != col) {
			if (njh::in(col, groupData_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << " Error, have grouping: " << col
						<< " entered more than once in the column names in "
						<< groupingsFile_ << "\n";
				throw std::runtime_error { ss.str() };
			}
			groupData_.emplace(col, std::make_unique<GroupMetaData>(col));
			VecStr grouping = groupsTab.getColumn(col);
			for (const auto pos : iter::range(grouping.size())) {
				if (njh::in(samples[pos], samples_)) {
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
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["groupingsFile_"] = njh::json::toJson(groupingsFile_);
	ret["samples_"] = njh::json::toJson(samples_);
	ret["missingSamples_"] = njh::json::toJson(missingSamples_);
	ret["missingMetaForSamples_"] = njh::json::toJson(missingMetaForSamples_);
	ret["groupData_"] = njh::json::toJson(groupData_);
	return ret;
}


void MultipleGroupMetaData::writeOutMetaFile(const OutOptions & outOpts, const std::set<std::string> & samples) const {
	OutputStream out(outOpts);
	auto metaFields = getVectorOfMapKeys(groupData_);
	njh::sort(metaFields);

	VecStr missingSamples;
	for(const auto & samp : samples) {
		if(njh::notIn(samp, samples_)) {
			missingSamples.emplace_back(samp);
		}
	}

	if(!missingSamples.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "missing the following samples when writing: " << njh::conToStr(missingSamples, ",") << "\n";
		ss << "options are: " << njh::conToStr(samples_, ",") << "\n";
		throw std::runtime_error { ss.str() };
	}

	out << "samples" << "\t" << njh::conToStr(metaFields, "\t") << std::endl;
	for(const auto & samp : samples) {
		out << samp;
		for(const auto & meta : metaFields) {
			out << "\t" << groupData_.at(meta)->sampleToSubGroup_.at(samp);
		}
		out << std::endl;
	}
}

void MultipleGroupMetaData::writeOutMetaFile(const OutOptions & outOpts) const {
	writeOutMetaFile(outOpts, samples_);
}



MultipleGroupMetaData MultipleGroupMetaData::fromJson(
		const Json::Value & jsonValue) {
	njh::json::MemberChecker checker(jsonValue);
	checker.failMemberCheckThrow(VecStr { "groupingsFile_", "samples_",
			"missingSamples_", "missingMetaForSamples_", "groupData_" },
			__PRETTY_FUNCTION__);

	MultipleGroupMetaData ret(jsonValue["groupingsFile_"].asString());
	auto jsonToStr = [](const Json::Value & v) {return v.asString();};
	ret.samples_ = njh::json::jsonArrayToSet<std::string>(jsonValue["samples_"],
			jsonToStr);
	ret.missingSamples_ = njh::json::jsonArrayToSet<std::string>(
			jsonValue["missingSamples_"], jsonToStr);
	ret.missingMetaForSamples_ = njh::json::jsonArrayToSet<std::string>(
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
		ret += njh::pasteAsStr(subGroupCount.first, ":", subGroupCount.second, ";");
	}
	return ret;
}

std::string MultipleGroupMetaData::GroupPopInfo::groupFracsStr() const{
	std::string ret { "" };
	for (const auto & subGroupCount : subGroupCounts_) {
		ret += njh::pasteAsStr(subGroupCount.first, ":",
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
		ss << "Options are : " << njh::conToStr(njh::getVecOfMapKeys(groupData_)) << "\n";
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
	return njh::in(metaField, groupData_);
}

bool MultipleGroupMetaData::hasSample(const std::string & sample) const{
	return njh::in(sample, samples_);
}


MetaDataInName MultipleGroupMetaData::getMetaForSample(const std::string & name) const{
	return getMetaForSample(name, getVectorOfMapKeys(groupData_));
}

MetaDataInName MultipleGroupMetaData::getMetaForSample(const std::string & name, const VecStr & fields) const{
	if(!njh::in(name, samples_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no sample: " << name << "\n";
		ss << "Options are: " << njh::conToStr(samples_);
		throw std::runtime_error{ss.str()};
	}
	checkForFieldsThrow(fields);

	MetaDataInName ret;
	for(const auto & gMeta : groupData_){
		if(njh::in(gMeta.first, fields)){
			ret.addMeta(gMeta.first, gMeta.second->getGroupForSample(name));
		}
	}
	return ret;
}

MetaDataInName MultipleGroupMetaData::genNaMeta( const VecStr & fields) const{
	checkForFieldsThrow(fields);

	MetaDataInName ret;
	for(const auto & gMeta : groupData_){
		if(njh::in(gMeta.first, fields)){
			ret.addMeta(gMeta.first, "NA");
		}
	}
	return ret;
}

table MultipleGroupMetaData::leftJoinWithMeta(const table & sampleTable,
		const std::string & sampleColumn) const{
	sampleTable.checkForColumnsThrow({sampleColumn}, __PRETTY_FUNCTION__);
	auto sampleColPos = sampleTable.getColPos(sampleColumn);
	auto metaLevels = getVectorOfMapKeys(groupData_);
	njh::sort(metaLevels);
	table outTab(concatVecs(sampleTable.columnNames_, metaLevels));
	for(auto row : sampleTable.content_){
		for(const auto & field : metaLevels){
			if(hasSample(row[sampleColPos])){
				row.emplace_back(groupData_.at(field)->getGroupForSample(row[sampleColPos]));
			}else{
				row.emplace_back("NA");
			}
		}
		outTab.addRow(row);
	}
	return outTab;
}


table MultipleGroupMetaData::leftJoinWithMeta(const table & sampleTable) const {
	std::regex sampPat(
			R"(sample(s\b|\b))", std::regex_constants::ECMAScript | std::regex_constants::icase );
	std::vector<uint32_t> possibleMatches = getPositionsMatchingPattern(sampleTable.columnNames_, sampPat);
	if(possibleMatches.size() == 0){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error found no column named that could match possible sample names, found only: "
				<< njh::conToStr(sampleTable.columnNames_, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}else if(possibleMatches.size()  > 1){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error found too many columns that match possible sample names, found : "
				<< njh::conToStr(getTargetsAtPositions(sampleTable.columnNames_, possibleMatches), ",") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return leftJoinWithMeta(sampleTable, sampleTable.columnNames_[possibleMatches.front()]);
}


}  // namespace njhseq
