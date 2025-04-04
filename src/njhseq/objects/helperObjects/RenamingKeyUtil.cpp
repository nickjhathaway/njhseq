//
// Created by Nicholas Hathaway on 4/4/25.
//

#include "RenamingKeyUtil.hpp"

namespace njhseq {

void RenamingKeyUtil::RenamingKeyUtilPars::setOptions(seqSetUp & setUp) {
	setUp.setOption(nameKeyFnp_, "--nameKeyFnp", "name Key Fnp, tab-delimited file either no column file with col1 being old name and col2 being new name or can supply which column names are old and new names with flags --oldNameColumnName and --newNameColumnName", true);
	setUp.setOption(old_name_column_name_, "--oldNameColumnName", "the name of a column to be the old name");
	setUp.setOption(new_name_column_name_, "--newNameColumnName", "the name of the column to be the new/replacement name", "" != old_name_column_name_);
	if (!new_name_column_name_.empty()) {
		setUp.addWarning(njh::pasteAsStr("if setting --newNameColumnName, also need to set --oldNameColumnName"));
		setUp.failed_ = true;
	}
}

RenamingKeyUtil::RenamingKeyUtil(const RenamingKeyUtilPars& pars): pars_(pars) {
	if (pars_.new_name_column_name_.empty()) {
		table keyTab(pars_.nameKeyFnp_, "\t", false);
		if (keyTab.nCol() != 2) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << pars_.nameKeyFnp_ <<
					" should have two columns, not: " << keyTab.nCol() << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto old_name_col_pos = 0;
		auto new_name_col_pos = 1;
		for (const auto& row: keyTab) {
			if (njh::in(row[old_name_col_pos], nameKeyMap_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have " << row[old_name_col_pos] << " in replacement map"
						<< "\n";
				throw std::runtime_error{ss.str()};
			}
			if (njh::in(row[new_name_col_pos], replacementNameToOriginalName_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have replacement name: " << row[new_name_col_pos] <<
						" for " << replacementNameToOriginalName_[row[new_name_col_pos]] << "\n";
				throw std::runtime_error{ss.str()};
			}
			nameKeyMap_[row[old_name_col_pos]] = row[new_name_col_pos];
			replacementNameToOriginalName_[row[new_name_col_pos]] = row[old_name_col_pos];
		}
	} else {
		table keyTab(pars_.nameKeyFnp_, "\t", true);
		keyTab.checkForColumnsThrow(VecStr{pars_.old_name_column_name_, pars_.new_name_column_name_}, __PRETTY_FUNCTION__);
		auto old_name_col_pos = keyTab.getColPos(pars_.old_name_column_name_);
		auto new_name_col_pos = keyTab.getColPos(pars_.new_name_column_name_);
		for (const auto& row: keyTab) {
			if (njh::in(row[old_name_col_pos], nameKeyMap_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have " << row[old_name_col_pos] << " in replacement map"
						<< "\n";
				throw std::runtime_error{ss.str()};
			}
			if (njh::in(row[new_name_col_pos], replacementNameToOriginalName_)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have replacement name: " << row[new_name_col_pos] <<
						" for " << replacementNameToOriginalName_[row[new_name_col_pos]] << "\n";
				throw std::runtime_error{ss.str()};
			}
			nameKeyMap_[row[old_name_col_pos]] = row[new_name_col_pos];
			replacementNameToOriginalName_[row[new_name_col_pos]] = row[old_name_col_pos];
		}
	}
}



} // njhseq