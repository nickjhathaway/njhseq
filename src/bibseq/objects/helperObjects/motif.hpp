#pragma once
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
/*
 * motif.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: nickhathaway
 */


#include "bibseq/utils.h"

namespace bibseq {


struct motifSubUnit{
	//constructors
	motifSubUnit():aas_(getUpperCaseLetters()), inclusive_(true){
		setScoreArray();
	}
	motifSubUnit(const std::vector<char> & aas, bool inclusive): aas_(aas), inclusive_(inclusive){
		setScoreArray();
	}
	motifSubUnit(const std::string & motifSub, bool inclusive): inclusive_(inclusive){
		for(const auto & c : motifSub){
			aas_.emplace_back(c);
		}
		setScoreArray();
	}
	//members
	std::vector<char> aas_;
	bool inclusive_;
	std::array<uint32_t, 26> score_;
	//functions
	void setScoreArray();
	uint32_t scoreChar(char c) const;
};

class motif {
public:
	/**@brief Constructor, in string should be a format similar to  N{P}[ST]{P},
	 *  where [ST] means either S or T and {P} means anything but P
	 *
	 * @param inMotif In protein string
	 */
	motif(const std::string & inMotif): motifOriginal_(inMotif){
		processMotif();
	}

	//members
	std::string motifOriginal_;

	std::map<uint32_t, motifSubUnit> motifUnits_;

	//functions
	void processMotif();
	motifSubUnit processInclusion(uint32_t start, uint32_t stop);
	motifSubUnit processExclusion(uint32_t start, uint32_t stop);
	uint32_t scoreMotif(const std::string & possibleMotif);
	uint32_t scoreMotif(const std::string::const_iterator & targetBegin,
			const std::string::const_iterator & targetEnd );

	bool passMotifParameter(const std::string & possibleMotif, uint32_t scoreCutOff);
	std::vector<uint32_t> findPositions(const std::string & wholeProtein, uint32_t scoreCutOff);
	std::vector<uint32_t> findPositionsFull(const std::string & wholeProtein,
			uint32_t allowableErrors = 0);
	std::vector<uint32_t> findPositionsFull(const std::string & wholeProtein,
			uint32_t allowableErrors, uint32_t start, uint32_t stop) const;

};

} /* namespace bib */
#ifndef NOT_HEADER_ONLY
#include "motif.cpp"
#endif
