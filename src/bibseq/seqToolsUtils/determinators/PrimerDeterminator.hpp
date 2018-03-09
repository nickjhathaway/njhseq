#pragma once
/*
 * PrimerDeterminator.hpp
 *
 *  Created on: Jun 10, 2015
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
#include "bibseq/objects/dataContainers/tables/table.hpp"
#include "bibseq/alignment/aligner.h"
#include "bibseq/objects/seqObjects/Paired/PairedRead.hpp"

namespace bibseq {

class PrimerDeterminator{
public:
	struct primerInfo {
		primerInfo() {}
		primerInfo(const std::string & primerPairName,
				const std::string & forwardPrimer, const std::string &reversePrimer );
		std::string primerPairName_;
		std::string forwardPrimer_;      /**< 5`-3` direction */
		seqInfo forwardPrimerInfo_;      /**< 5`-3` direction */
		seqInfo forwardPrimerInfoRevDir_;/**< 3`-5` direction */
		std::string reversePrimer_;      /**< 5`-3` direction */
		seqInfo reversePrimerInfo_;			 /**< 3`-5` direction */
		seqInfo reversePrimerInfoForDir_;/**< 5`-3` direction */
	};

	explicit PrimerDeterminator(const table & primers);

	explicit PrimerDeterminator(const std::unordered_map<std::string, primerInfo> & primers);

//	void addPrimerInfo(const std::string & primerName,
//			const std::string & forwardPrimer, const std::string & reversePrimer);

	bool containsTarget(const std::string & targetName) const;


	std::map<std::string, primerInfo> primers_;

	size_t getMaxPrimerSize() const;

	template<typename T>
	std::string determineForwardPrimer(T & read, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool primerToLowerCase){
		return determineForwardPrimer(getSeqBase(read), withinPos, alignerObj, allowable, primerToLowerCase);
	}
	std::string determineForwardPrimer(seqInfo & info, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool primerToLowerCase);

	template<typename T>
	std::string determineWithReversePrimer(T & read, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool primerToLowerCase){
		return determineWithReversePrimer(getSeqBase(read), withinPos, alignerObj, allowable, primerToLowerCase);
	}

	std::string determineWithReversePrimer(seqInfo & info, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool primerToLowerCase);


	template<typename T>
	bool checkForReversePrimer(T & read, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
			uint32_t within, bool trimExtra) {
		return checkForReversePrimer(getSeqBase(read), primerName, alignObj, allowable, reversePrimerToLowerCase, within, trimExtra);
	}

	bool checkForReversePrimer(seqInfo & info, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
      uint32_t within, bool trimExtra);

	template<typename T>
	bool checkForForwardPrimerInRev(T & read, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
			uint32_t within, bool trimExtra) {
		return checkForForwardPrimerInRev(getSeqBase(read), primerName, alignObj, allowable, reversePrimerToLowerCase, within,trimExtra);
	}

	bool checkForForwardPrimerInRev(seqInfo & info, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
      uint32_t within, bool trimExtra);
};

} /* namespace bibseq */


