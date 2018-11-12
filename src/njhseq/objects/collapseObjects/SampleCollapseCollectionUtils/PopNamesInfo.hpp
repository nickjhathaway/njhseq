#pragma once
/*
 * PopNamesInfo.hpp
 *
 *  Created on: Jul 31, 2016
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
#include "njhseq/utils.h"

namespace njhseq {

/**@brief simple class to hold the population name and the samples that go in it
 *
 */
class PopNamesInfo {

public:
	/**@brief construct with pop name and samples
	 *
	 * @param populationName population name
	 * @param samples the samples
	 */
	PopNamesInfo(std::string populationName, std::set<std::string> samples);
	/**@brief construct with population name and samples
	 *
	 * @param populationName the population name
	 * @param samples samples in vector that will be converted into set to get rid of duplicate names
	 */
	PopNamesInfo(std::string populationName, VecStr samples);

	std::string populationName_; /**< population name*/
	std::set<std::string> samples_; /**< samples in set so no duplicates*/

	/**@brief check to see if sample is samples_
	 *
	 * @param sample the sample to check for
	 * @return true if sample is in samples_
	 */
	bool hasSample(const std::string & sample) const;

	/**@brief output json object representation of class
	 *
	 * @return Json::Value with info of class
	 */
	Json::Value toJson() const;
private:
	/**@brief check formating of population name, it cannot contain ".", will throw if fails format requirments
	 *
	 */
	void checkPopNameThrow() const;



};


}  // namespace njhseq






