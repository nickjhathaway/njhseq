#pragma once

//
// Created by Nicholas Hathaway on 4/4/25.
//

#include "njhseq/common/allSystemIncludes.h"
#include "njhseq/programUtils/seqSetUp.hpp"



namespace njhseq {
/**
 * @brief A utility to help rename fields or identifiers
 */
class RenamingKeyUtil {

public:

	struct RenamingKeyUtilPars {
		bfs::path nameKeyFnp_;
		std::string old_name_column_name_;
		std::string new_name_column_name_;

		/**
		* @brief Set the default options, uses flags: --oldNameColumnName, --nameKeyFnp, --newNameColumnName
		 * @param setUp a program seq set up object to set the options
		 */
		void setOptions(seqSetUp & setUp);
	};

	/**
	 * @brief construct with default parameters the renaming key, will read in file and set up key, will throw if errors
	 * @param pars the parameters to setup the renaming key, tab delimited files, if pars.old_name_column_name_ and pars.new_name_column_name_ blank, will assume no column name 2 column file 1) old name, 2) new name
	 */
	RenamingKeyUtil(const RenamingKeyUtilPars& pars);

	RenamingKeyUtilPars pars_;
	std::unordered_map<std::string, std::string> nameKeyMap_;
	std::unordered_map<std::string, std::string> replacementNameToOriginalName_;
};

} // njhseq
