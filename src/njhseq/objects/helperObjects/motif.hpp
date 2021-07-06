#pragma once
/*
 * motif.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: nickhathaway
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

#include "njhseq/utils.h"

namespace njhseq {




class motif {
	struct motifSubUnit{
		//constructors
		motifSubUnit();
		motifSubUnit(const std::vector<char> & aas, bool inclusive);
		motifSubUnit(const std::string & motifSub, bool inclusive);
	private:
		//members
		std::vector<char> aas_;
		bool inclusive_;
		std::array<uint32_t, 26> score_;
		//functions
		void setScoreArray();
		uint32_t otherScore_{std::numeric_limits<uint32_t>::max()};

	public:
		uint32_t scoreCharNoCheck(char c) const;
		uint32_t scoreCharCheck(char c) const;
		friend class motif;
	};
public:
	/**@brief Constructor, in string should be a format similar to  N{P}[ST]{P},
	 *  where [ST] means either S or T and {P} means anything but P
	 *
	 * @param inMotif In protein string
	 */
	motif(const std::string & inMotif);

	//members
	std::string motifOriginal_;
private:
	std::map<uint32_t, motifSubUnit> motifUnits_;
	//functions
	void processMotif();
	motifSubUnit processInclusion(uint32_t start, uint32_t stop);
	motifSubUnit processExclusion(uint32_t start, uint32_t stop);
public:
	uint32_t scoreMotif(const std::string & possibleMotif)const;
	uint32_t scoreMotif(const std::string::const_iterator & targetBegin,
			const std::string::const_iterator & targetEnd )const;

	bool passMotifParameter(const std::string & possibleMotif,
			uint32_t scoreCutOff) const;

	std::vector<size_t> findPositions(const std::string & wholeProtein,
			uint32_t scoreCutOff) const;

	std::vector<size_t> findPositionsFull(const std::string & wholeProtein,
			uint32_t allowableErrors) const;
	/**@brief Look for motif between these positions
	 *
	 * @param wholeProtein The string to search for the motif in
	 * @param allowableErrors How many errors to allow in the motif
	 * @param start The position from which to start to search
	 * @param stop The end position to stop searching in, the whole motif should come before this postion
	 * @return All the positions the motif was found within
	 */
	std::vector<size_t> findPositionsFull(const std::string & wholeProtein,
			uint32_t allowableErrors, size_t start, size_t stop) const;

	bool frontPassNoCheck(const std::string & wholeProtein,
			uint32_t allowableErrors) const;

	std::vector<size_t> findPositionsSubSets(
			const std::string & wholeProtein,
			uint32_t allowableErrors,
			size_t start, size_t stop,
			uint32_t motifStart, uint32_t motifEnd) const;

	std::vector<size_t> findPositionsSubSetsBest(
			const std::string & wholeProtein,
			uint32_t allowableErrors,
			size_t start, size_t stop,
			uint32_t motifStart, uint32_t motifEnd) const;

	size_t size()const;
};

} /* namespace njhseq */

