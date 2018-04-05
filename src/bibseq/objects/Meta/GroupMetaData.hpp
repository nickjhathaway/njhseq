#pragma once
/*
 * GroupMetaData.hpp
 *
 *  Created on: Jul 30, 2016
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
#include "bibseq/utils.h"
#include "bibseq/objects/dataContainers/tables/table.hpp"

namespace bibseq {
/**@brief class to hold samples for meta data group
 *
 */
class GroupMetaData {
public:

	/**@brief construct with the sample master group
	 *
	 * @param name
	 */
	GroupMetaData(const std::string & name);

	std::string name_; /**< name of the master group*/
	std::set<std::string> subGroupsLevels_; /**< all the names for the sub groups */
	std::unordered_map<std::string, std::set<std::string>> subGroupToSamples_; /**< a map with group for key and the samples that are in that group as value*/
	std::unordered_map<std::string, std::string> sampleToSubGroup_; /**< a map with sample as key and sub group as value*/

	/**@brief add a sample and it's corresponding group
	 *
	 * @param sample a sample name
	 * @param subGroup the sub group the sample belongs to
	 */
	void addSampGroup(const std::string & sample, const std::string & subGroup);
	/**@brief set the sub groups levels
	 *
	 */
	void setSubGroupsLevels();

	/**@brief get the sub group for the given sample
	 *
	 * @param sample the sample to get the group for
	 * @return the group the sample is in
	 */
	std::string getGroupForSample(const std::string & samp) const;

	/**@brief convert to json
	 *
	 * @return
	 */
	Json::Value toJson() const;

	/**@brief get the samples stored in meta data
	 *
	 * @return a vector of the sample names
	 */
	VecStr getSampleNames() const;

	/**@brief create a GroupMetaData object from json (not declared as constructor due to overload rules)
	 *
	 * @param jsonValue the json value holding the info for the object
	 * @return a GroupMetaData object
	 */
	static GroupMetaData fromJson(const Json::Value & jsonValue);

};



}  // namespace bibseq



