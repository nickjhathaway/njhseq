#pragma once
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


