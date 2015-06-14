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
 * PrimerDeterminator.hpp
 *
 *  Created on: Jun 10, 2015
 *      Author: nick
 */

#include "bibseq/objects/dataContainers/table.hpp"
#include "bibseq/alignment/aligner.hpp"


namespace bibseq {

class PrimerDeterminator{
public:
	PrimerDeterminator(const table & primers) ;

	struct primerInfo {
		primerInfo() {}
		primerInfo(const std::string & primerPairName,
				const std::string & forwardPrimer, const std::string &reversePrimer );
		std::string primerPairName_;
		std::string forwardPrimer_;
		seqInfo forwardPrimerInfo_;
		seqInfo forwardPrimerInfoRevDir_;
		std::string reversePrimer_;
		seqInfo reversePrimerInfo_;
		seqInfo reversePrimerInfoForDir_;
	};

	std::map<std::string, primerInfo> primers_;

	template<typename T>
	std::string determineForwardPrimer(T & read, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool forwardPrimerToLowerCase,
			bool weighHomopolyers){
		return determineForwardPrimer(read.seqBase_, withinPos, alignerObj, allowable, forwardPrimerToLowerCase, weighHomopolyers);
	}
	std::string determineForwardPrimer(seqInfo & info, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool forwardPrimerToLowerCase,
			bool weighHomopolyers);

	template<typename T>
	std::string determineWithReversePrimer(T & read, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool forwardPrimerToLowerCase,
			bool weighHomopolyers){
		return determineWithReversePrimer(read.seqBase_, withinPos, alignerObj, allowable, forwardPrimerToLowerCase, weighHomopolyers);
	}

	std::string determineWithReversePrimer(seqInfo & info, uint32_t withinPos,
			aligner & alignerObj, const comparison & allowable, bool forwardPrimerToLowerCase,
			bool weighHomopolyers);

	template<typename T>
	bool checkForReversePrimer(T & read, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
			bool weighHomopolyers) {
		return checkForReversePrimer(read.seqBase_, primerName, alignObj, allowable, reversePrimerToLowerCase, weighHomopolyers);
	}

	bool checkForReversePrimer(seqInfo & info, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
      bool weighHomopolyers);

	template<typename T>
	bool checkForForwardPrimerInRev(T & read, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
			bool weighHomopolyers) {
		return checkForForwardPrimerInRev(read.seqBase_, primerName, alignObj, allowable, reversePrimerToLowerCase, weighHomopolyers);
	}

	bool checkForForwardPrimerInRev(seqInfo & info, const std::string & primerName,
			aligner & alignObj, const comparison & allowable, bool reversePrimerToLowerCase,
      bool weighHomopolyers);
};

} /* namespace bibseq */


