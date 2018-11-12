#pragma once
/*
 * kmerMap.hpp
 *
 *  Created on: Dec 2, 2015
 *      Author: nick
 */
//
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

#include "njhseq/objects/kmer/kmer.hpp"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"



namespace njhseq {

/**@brief simple struct to hold kmer with it's position
 *
 */
struct kmerWithPos{
	std::string kmer_;/**< Kmer string */
	uint32_t kPos_;/**< Position of kmer */
};

/**@brief Return a kmer centered around a position and the position of the kmer (adjust for when pos can't be centered near the ends)
 *
 *@todo potential pitfall is you could get a different size kmer if kLength is larger than str and a potential seg fault if str is empty
 *
 * @param pos The position to center the kmer on
 * @param kLength the length of the kmer
 * @param str the string to take the kmer from
 * @return a kmerWithPos object that has the kmer and it's position
 */
inline kmerWithPos getKmerPos(uint32_t pos, uint32_t kLength, const std::string & str) {
	//Intentionally rounded down
	uint32_t halfSize = kLength / 2;
	//determine position that centers pos in the middle of the kmer
	if (halfSize > pos) {
		pos = 0;
	} else if (pos + halfSize >= str.size()) {
		pos = str.size() - kLength;
	} else {
		pos = pos - halfSize;
	}
	return kmerWithPos { str.substr(pos, kLength), pos };
}


/**@brief base class to build KmerMaps for determining frequency of kmers
 *
 */
class KmerMapBase {

public:
	KmerMapBase(size_t kLength, uint32_t runCutOff);
	KmerMapBase(size_t kLength);
	KmerMapBase();
	size_t kLength_;
	uint32_t runCutOff_;

	virtual void increaseCounts(const std::string & str) = 0;
	virtual void increaseCounts(const std::string & str, uint32_t readCount) = 0;

	virtual void increaseCountsExtra(const std::string & str) = 0;
	virtual void increaseCountsExtra(const std::string & str, const std::string & uid, uint32_t readCount) = 0;

	virtual bool isKmerLowFreq(const kmerWithPos & k) const = 0;
	virtual bool isKmerLowFreq(const kmerWithPos & k, uint32_t runCutOff_) const = 0;

	virtual uint32_t getKmerFreq(const kmerWithPos & k) const = 0;

	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;

	virtual ~KmerMapBase() {}
};


class KmerMap: public KmerMapBase {
public:

	KmerMap(size_t kLength, uint32_t runCutOff,
			const std::unordered_map<std::string, kmer> & kmers);
	KmerMap(size_t kLength, uint32_t runCutOff);
	KmerMap(size_t kLength);
	KmerMap();
	std::unordered_map<std::string, kmer> kmers_;
	virtual void increaseCounts(const std::string & str);
	virtual void increaseCounts(const std::string & str, uint32_t readCount);

	virtual void increaseCountsExtra(const std::string & str);
	virtual void increaseCountsExtra(const std::string & str, const std::string & uid,
			uint32_t readCount);

	virtual bool isKmerLowFreq(const kmerWithPos & k) const;
	virtual bool isKmerLowFreq(const kmerWithPos & k, uint32_t runCutOff) const;

	virtual uint32_t getKmerFreq(const kmerWithPos & k) const;

	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;


	virtual ~KmerMap();
};

class KmerMapByPos: public KmerMapBase {
public:

	KmerMapByPos(size_t kLength, uint32_t runCutOff,
			const std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> & kmers);
	KmerMapByPos(size_t kLength, uint32_t runCutOff);
	KmerMapByPos(size_t kLength);
	KmerMapByPos();
	std::unordered_map<uint32_t, std::unordered_map<std::string, kmer>> kmers_;
	static std::pair<uint32_t, uint32_t> determineExpandPos(uint32_t expandSize,
			uint32_t pos, const std::string & seq);
	virtual void increaseCounts(const std::string & str);
	virtual void increaseCounts(const std::string & str, uint32_t readCount);

	virtual void increaseCountsExtra(const std::string & str);
	virtual void increaseCountsExtra(const std::string & str, const std::string & uid,
			uint32_t readCount);

	virtual void increaseCountsExp(const std::string & str, uint32_t expandSize);
	virtual void increaseCountsExp(const std::string & str, uint32_t expandSize,
			uint32_t readCount);

	virtual void increaseCountsExtraExp(const std::string & str, uint32_t expandSize);
	virtual void increaseCountsExtraExp(const std::string & str, uint32_t expandSize,
			const std::string & uid, uint32_t readCount);

	virtual bool isKmerLowFreq(const kmerWithPos & k) const;
	virtual bool isKmerLowFreq(const kmerWithPos & k, uint32_t runCutOff) const;

	virtual uint32_t getKmerFreq(const kmerWithPos & k) const;

	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;

	virtual ~KmerMapByPos();
};

class KmerMaps {
public:

	KmerMaps(size_t kLength, uint32_t runCutOff);
	KmerMaps(size_t kLength);
	KmerMaps();
	uint32_t kLength_;
	uint32_t runCutOff_;
	bool kmerSet_ = false;
	bool kmersByPosition_ = true;
	std::shared_ptr<KmerMapBase> kmers_;
	std::shared_ptr<KmerMap> kmersNoPos_;
	std::shared_ptr<KmerMapByPos> kmersByPos_;

	void setKmers(bool kmersByPosition);
	bool isKmerLowFreq(const kmerWithPos & k) const;
	bool isKmerLowFreq(const kmerWithPos & k, uint32_t runCutOff) const;

	void setAllRunCutOffs(uint32_t runCutOff);

	void resetMaps();

	/**@brief convert to json representation
	 *
	 * @return Json::Value object
	 */
	Json::Value toJson() const;

	template<typename T>
	void setAllKmers(const std::vector<T> & reads, bool kmersByPosition,
			uint32_t kmerExpandSize) {
		resetMaps();
		for (const auto & read : reads) {
			kmersNoPos_->increaseCounts(getSeqBase(read).seq_, getSeqBase(read).cnt_);
			kmersByPos_->increaseCountsExp(getSeqBase(read).seq_, kmerExpandSize,
					getSeqBase(read).cnt_);
		}
		setKmers(kmersByPosition);
	}

	template<typename T>
	void setAllKmers(const std::vector<T> & reads, bool kmersByPosition) {
		resetMaps();
		for (const auto & read : reads) {
			kmersNoPos_->increaseCounts(getSeqBase(read).seq_, getSeqBase(read).cnt_);
			kmersByPos_->increaseCounts(getSeqBase(read).seq_, getSeqBase(read).cnt_);
		}
		setKmers(kmersByPosition);
	}

	template<typename T>
	void setKmersByPos(const std::vector<T> & reads, bool kmersByPosition) {
		resetMaps();
		for (const auto & read : reads) {
			kmersByPos_->increaseCounts(getSeqBase(read).seq_, getSeqBase(read).cnt_);
		}
		setKmers(true);
	}

	template<typename T>
	void setKmersAnywhere(const std::vector<T> & reads) {
		resetMaps();
		for (const auto & read : reads) {
			kmersNoPos_->increaseCounts(getSeqBase(read).seq_, getSeqBase(read).cnt_);
		}
		setKmers(false);
	}

};

template<typename T>
KmerMaps indexKmers(const std::vector<T> & reads, uint32_t kLenth,
		uint32_t runCutOff, bool kmersByPosition, bool expandKmerPos,
		uint32_t expandKmerSize) {
	KmerMaps otherMap = KmerMaps(kLenth, runCutOff);
	if (kmersByPosition) {
		otherMap.setAllKmers(reads, kmersByPosition, expandKmerSize);
	} else {
		otherMap.setAllKmers(reads, kmersByPosition);
	}
	return otherMap;
}


}  // namespace njhseq


