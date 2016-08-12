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

MultipleGroupMetaData::MultipleGroupMetaData(const bfs::path & groupingsFile,
		const std::set<std::string> & availableSamples) :
		groupingsFile_(groupingsFile) {
	table groupsTab(groupingsFile_.string(), "\t", true);
	if (!groupsTab.containsColumn("Sample")) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, " << groupingsFile
				<< " must contain the column \"Sample\"" << "\n";
		ss << "Current columns are: " << vectorToString(groupsTab.columnNames_, ",")
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	if (groupsTab.nCol() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, " << groupingsFile
				<< " should contain at least two columns, current size is : "
				<< groupsTab.nCol() << "\n";
		throw std::runtime_error { ss.str() };
	}
	VecStr samples = groupsTab.getColumn("Sample");
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
	for (const auto & col : groupsTab.columnNames_) {
		if ("Sample" != col) {
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
			groupData_[col]->setSubGroupsLevels();
		}
	}
	for (const auto & samp : availableSamples) {
		if (!bib::in(samp, samples)) {
			missingMetaForSamples_.insert(samp);
		}
	}
}




}  // namespace bibseq


