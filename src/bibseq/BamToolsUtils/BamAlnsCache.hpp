#pragma once
/*
 * BamAlnsCache.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: nick
 */
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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


#include <api/BamAlignment.h>
#include "bibseq/utils.h"

namespace bibseq {


/**@brief Simple class to hold BamTools::BamAlignment in an unordered_map
 *
 */
class BamAlnsCache{
private:
	std::unordered_map<std::string, std::shared_ptr<BamTools::BamAlignment>> cache_; /**< The cache of aligmnets */
public:

	/**@brief Add alignment to the cache
	 *
	 * @param aln The alignment to add
	 */
	void add(const BamTools::BamAlignment & aln);
	/**@brief Remove an alignment from the cache
	 *
	 * @param name The name of the alignment to remove
	 */
	void remove(const std::string & name);

	/**@brief Check to see if cache contains alignment by name
	 *
	 * @param name the name of the alignment to check for
	 * @return true if cache has alignmnet
	 */
	bool has(const std::string & name)const;

	/**@brief get a shared ptr to the alignment, nullptr is returned if cache doesn't contain alignment
	 *
	 * @param name The name of the alignment to get
	 * @return
	 */
	std::shared_ptr<BamTools::BamAlignment> get(const std::string & name);

	/**@brief Get a vector of the names of the alignments in the cache
	 *
	 * @return A vector of names
	 */
	VecStr getNames() const;

	/**@brief Get the min alignment position for each refId
	 *
	 * @return a map with the min pos by refId
	 */
	std::unordered_map<int32_t, int32_t> getMinPos() const;
	/**@brief Get the max alignment position for each refId
	 *
	 * @return a map with the max pos by refId
	 */
	std::unordered_map<int32_t, int32_t> getMaxEndPos() const;

	using size_type = std::unordered_map<std::string, BamTools::BamAlignment>::size_type;

	template<typename T>
	friend typename T::size_type len(const T & cache);

};


template<>
inline BamAlnsCache::size_type len(const BamAlnsCache & cache) {
	return cache.cache_.size();
}


}  // namespace bibseq


