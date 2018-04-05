#pragma once
//
//  baseReadObject.hpp
//
//  Created by Nicholas Hathaway on 8/31/13.
//
//
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
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"

namespace bibseq {

/**@brief Simple object that just has a seqInfo as a member,
 *
 * Most seq objects should inherit at least from this so that they all have seqBase_ as a member
 *
 */
class baseReadObject {

public:
	/**@brief Constructs an empty object, seqBase_ is also just an empty seqInfo object
	 *
	 */
	baseReadObject();
	/**@brief Construct with seqInfo object that will used to set seqBase_
	 *
	 * @param seqBase The seqInfo object to set as the base of the object
	 */
	baseReadObject(const seqInfo& seqBase);

	seqInfo seqBase_; /**< The actual workhorse of the class, has name, seq and count */

	/**@brief Convert to jsoncpp object
	 *
	 * @return a jsoncpp value object representing the class
	 */
	virtual Json::Value toJson() const;

	virtual ~baseReadObject();

	using size_type = seqInfo::size_type;
};

/**@brief template specfication of the len template function to return the length of the sequence stored in seqBase_
 *
 * @param read object to get the length of
 * @return The length of the sequence
 */
template<>
inline baseReadObject::size_type len(const baseReadObject & read){
	return read.seqBase_.seq_.size();
}

}  // namespace bibseq


