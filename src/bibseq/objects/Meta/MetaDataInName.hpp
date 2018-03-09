#pragma once
/*
 * MetaDataInName.hpp
 *
 *  Created on: Jan 27, 2017
 *      Author: nick
 */

#include "bibseq/common.h"

namespace bibseq {


class MetaDataInName {
public:

	MetaDataInName();
	MetaDataInName(const std::string & str);

	std::unordered_map<std::string, std::string> meta_;

	template<typename T>
	void addMeta(const std::string & key, const T & val, bool replace = false) {
		if (containsMeta(key) && !replace) {
			std::stringstream ss;
			ss << "Error in " << bib::bashCT::boldBlack(__PRETTY_FUNCTION__)
					<< " attempting to add meta, " << key << ", that's already in meta_, use replace = true to replace"
					<< std::endl;
			throw std::runtime_error { ss.str() };
		} else {
			meta_[key] = estd::to_string(val);
		}
	}

	void addMeta(const MetaDataInName & otherMeta, bool replace);

	template<typename T>
	T getMeta(const std::string & key) const {
		return bib::lexical_cast<T>(getMeta(key));
	}

	void removeMeta(const std::string & metaField);

	void processNameForMeta(const std::string & name, bool replace);

	bool containsMeta(const std::string & key) const;
	void containsMetaThrow(const std::string & key, const std::string & funcName) const;

	std::string getMeta(const std::string & key) const;

	std::string createMetaName() const;
	std::string createMetaName(const std::function<bool(const std::string &, const std::string &)> & metaKeyPredSorter) const;

	void resetMetaInName(std::string & name,
			size_t pos = std::numeric_limits<size_t>::max()) const;

	static void removeMetaDataInName(std::string & name);

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

};


}  // namespace bibseq




