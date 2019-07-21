#pragma once
/*
 * PrimerDeterminator.hpp
 *
 *  Created on: Jun 10, 2015
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
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/alignment/aligner.h"
#include "njhseq/objects/seqObjects/Paired/PairedRead.hpp"

namespace njhseq {

class PrimerDeterminator {
public:

	struct PrimerPositionScore {
		PrimerPositionScore() {
		}
		PrimerPositionScore(const uint32_t start, const uint32_t end,
				const std::string & primerName, const comparison & comp) :
				start_(start), end_(end), primerName_(primerName), comp_(comp) {
		}
		uint32_t start_ { std::numeric_limits<uint32_t>::max()-1 };
		uint32_t end_ { std::numeric_limits<uint32_t>::max() };
		std::string primerName_;
		comparison comp_;
		double getNormalizedScore() {
			return comp_.alnScore_ / (end_ - start_);
		}
	};

	struct PrimerDeterminatorPars {
		comparison allowable_;

		bool primerToLowerCase_ { true };
		uint32_t primerWithin_ { 0 };
		uint32_t primerStart_ { 0 };
		bool trimExtra_ { false };
		bool checkComplement_ { false };
	};

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

		charCounter forwardPrimerInfoLetCounter_;
		charCounter forwardPrimerInfoRevDirLetCounter_;
		charCounter reversePrimerInfoLetCounter_;
		charCounter reversePrimerInfoForDirLetCounter_;

	};

	explicit PrimerDeterminator(const table & primers);

	explicit PrimerDeterminator(const std::unordered_map<std::string, primerInfo> & primers);

//	void addPrimerInfo(const std::string & primerName,
//			const std::string & forwardPrimer, const std::string & reversePrimer);

	bool containsTarget(const std::string & targetName) const;


	std::map<std::string, primerInfo> primers_;

	size_t getMaxPrimerSize() const;



	template<typename T>
	std::string determineForwardPrimer(T & read, const PrimerDeterminatorPars & pars, aligner & alignerObj){
		return determineForwardPrimer(getSeqBase(read), pars, alignerObj);
	}
	std::string determineForwardPrimer(seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj);

	template<typename T>
	std::string determineWithReversePrimer(T & read, const PrimerDeterminatorPars & pars, aligner & alignerObj){
		return determineWithReversePrimer(getSeqBase(read), pars, alignerObj);
	}
	std::string determineWithReversePrimer(seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj);


	template<typename T>
	bool checkForReversePrimer(T & read, const std::string & primerName, const PrimerDeterminatorPars & pars, aligner & alignObj) {
		return checkForReversePrimer(getSeqBase(read), primerName,pars, alignObj);
	}
	bool checkForReversePrimer(seqInfo & info, const std::string & primerName, const PrimerDeterminatorPars & pars, aligner & alignObj);

	template<typename T>
	bool checkForForwardPrimerInRev(T & read, const std::string & primerName, const PrimerDeterminatorPars & pars, aligner & alignObj) {
		return checkForForwardPrimerInRev(getSeqBase(read), primerName, pars, alignObj);
	}

	bool checkForForwardPrimerInRev(seqInfo & info, const std::string & primerName, const PrimerDeterminatorPars & pars, aligner & alignObj);


	//just getting positions
	template<typename T>
	PrimerPositionScore determineBestForwardPrimerPosFront(const T & read, const PrimerDeterminatorPars & pars, aligner & alignerObj){
		return determineBestForwardPrimerPosFront(getSeqBase(read), pars, alignerObj);
	}
	PrimerPositionScore determineBestForwardPrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj);


	template<typename T>
	PrimerPositionScore determineBestReversePrimerPosFront(const T & read, const PrimerDeterminatorPars & pars, aligner & alignerObj){
		return determineBestReversePrimerPosFront(getSeqBase(read), pars, alignerObj);
	}
	PrimerPositionScore determineBestReversePrimerPosFront(const seqInfo & info, const PrimerDeterminatorPars & pars, aligner & alignerObj);




};

} /* namespace njhseq */


