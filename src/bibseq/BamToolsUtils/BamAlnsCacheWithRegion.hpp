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


#include "bibseq/BamToolsUtils/BamAlnsCache.hpp"
#include "bibseq/objects/BioDataObject/GenomicRegion.hpp"


namespace bibseq {


/**@brief Simple class to hold BamTools::BamAlignment in an unordered_map and the region associated with
 *
 */
class BamAlnsCacheWithRegion : public BamAlnsCache{
private:
	std::unordered_map<std::string, std::shared_ptr<GenomicRegion>> regionCache_; /**< The cache of region assoicated with bam alignment cache */
	template<typename T>
	struct fail: std::false_type {
	};
public:

	/**@brief Add alignment to the cache
	 *
	 * @param aln The alignment to add
	 * @param region The region associated with alignment
	 */
	virtual void addWithRegion(const BamTools::BamAlignment & aln, const GenomicRegion & region);

	/**@brief Add alignment to the cache
	 *
	 * @param aln The alignment to add
	 * @param region The region associated with alignment
	 */
	virtual void addWithRegion(const BamTools::BamAlignment & aln, const std::shared_ptr<GenomicRegion> & region);


	template<typename T=bool>
	void add(const BamTools::BamAlignment & aln){
		static_assert (fail<T>::value, "Do not use BamAlnsCacheWithRegion::add(const BamTools::BamAlignment & aln)!, use BamAlnsCacheWithRegion::addWithRegion(const BamTools::BamAlignment & aln, const GenomicRegion & region)");
	}

	/**@brief Remove an alignment from the cache
	 *
	 * @param name The name of the alignment to remove
	 */
	virtual void remove(const std::string & name);

	/**@brief get a shared ptr to the region associated with the alignment name, nullptr is returned if cache doesn't contain alignment
	 *
	 * @param name The name of the region associated with the alignment name
	 * @return shared_ptr to region
	 */
	std::shared_ptr<GenomicRegion> getRegion(const std::string & name);


	using size_type = std::unordered_map<std::string, BamTools::BamAlignment>::size_type;

	template<typename T>
	friend typename T::size_type len(const T & cache);

	virtual  ~BamAlnsCacheWithRegion();
};


template<>
inline BamAlnsCacheWithRegion::size_type len(const BamAlnsCacheWithRegion & cache) {
	return cache.cache_.size();
}


}  // namespace bibseq


