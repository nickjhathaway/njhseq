#pragma once
/*
 * MetaDataInName.hpp
 *
 *  Created on: Jan 27, 2017
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
#include "njhseq/common.h"

namespace njhseq {


class MetaDataInName {
public:

	MetaDataInName();
	MetaDataInName(const std::string & str);

	std::unordered_map<std::string, std::string> meta_;

	template<typename T>
	void addMeta(const std::string & key, const T & val, bool replace = false) {
		if (containsMeta(key) && !replace) {
			std::stringstream ss;
			ss << "Error in " << njh::bashCT::bold << __PRETTY_FUNCTION__
					<< njh::bashCT::reset << " attempting to add meta, "
					<< njh::bashCT::bold << key << njh::bashCT::reset
					<< ", that's already in meta_, use replace = true to replace"
					<< std::endl;
			throw std::runtime_error { ss.str() };
		} else {
			meta_[key] = estd::to_string(val);
		}
	}

	void addMeta(const MetaDataInName & otherMeta, bool replace);

	template<typename T>
	T getMeta(const std::string & key) const {
		return njh::lexical_cast<T>(getMeta(key));
	}



	void removeMeta(const std::string & metaField);

	void processNameForMeta(const std::string & name, bool replace);

	bool containsMeta(const std::string & key) const;
	void containsMetaThrow(const std::string & key, const std::string & funcName) const;

	std::string getMeta(const std::string & key) const;

	std::string createMetaName() const;
	std::string createMetaName(const std::function<bool(const std::string &, const std::string &)> & metaKeyPredSorter) const;


	std::string pasteLevels(const std::string & sep = "") const;
	std::string pasteLevels(const VecStr & metalevels, const std::string & sep = "") const;


	void resetMetaInName(std::string & name,
			size_t pos = std::numeric_limits<size_t>::max()) const;

	static void removeMetaDataInName(std::string & name);
	static std::string removeMetaDataInNameRet(std::string name);

	static bool nameHasMetaData(const std::string & name);


	Json::Value toJson() const;

	template<typename MAP>
	static MetaDataInName mapToMeta(const MAP & m){
		MetaDataInName ret;
		for(const auto & p : m){
			ret.addMeta(estd::to_string(p.first), p.second);
		}
		return ret;
	}

	static MetaDataInName genMetaFromJson(const Json::Value & val);
};

template<>
inline bool MetaDataInName::getMeta(const std::string & key) const {
	return "true" == njh::strToLowerRet(getMeta(key));
}


}  // namespace njhseq




