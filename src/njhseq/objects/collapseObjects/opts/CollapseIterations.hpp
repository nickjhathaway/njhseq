#pragma once
/*
 * CollapseIterations.hpp
 *
 *  Created on: May 6, 2016
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

#include "njhseq/objects/collapseObjects/opts/IterPar.hpp"

namespace njhseq {

class CollapseIterations {
public:

	CollapseIterations(const std::string & parametersFilename, bool onPerId);

	CollapseIterations(std::map<uint32_t, IterPar> pars, bool onPerId);

	CollapseIterations(bool onPerId);

	CollapseIterations();

	std::map<uint32_t, IterPar> iters_;
	std::map<uint32_t, uint32_t> itersNums_;
	bool onPerId_ = false;
	bool onHqPerId_ = false;


	CollapseIterations& operator=(const CollapseIterations & other);
	CollapseIterations& operator=(CollapseIterations&& other);
	CollapseIterations(const CollapseIterations& other);

	void addIteration(uint32_t iterNum, std::vector<double> pars);
	void addIteration(uint32_t iterNum, const IterPar & par);
	void addIteration(const IterPar & par);

	void addIterations(std::vector<std::vector<double>> pars);

	void writePars(std::ostream & out)const;

	static CollapseIterations gen454ItDefaultPars(uint32_t stopCheck);
	static CollapseIterations genIlluminaDefaultPars(uint32_t stopCheck);
	static CollapseIterations genIlluminaDefaultParsCollapseHomopolymers(uint32_t stopCheck);

	static CollapseIterations gen454ItDefaultParsWithHqs(uint32_t stopCheck,
			uint32_t hqMismatches);
	static CollapseIterations genIlluminaDefaultParsWithHqs(uint32_t stopCheck,
			uint32_t hqMismatches);
	static CollapseIterations genIlluminaDefaultParsWithHqsCollapseHomopolymers(uint32_t stopCheck,
			uint32_t hqMismatches);

	static CollapseIterations genStrictNoErrorsDefaultPars(uint32_t stopCheck);
	static CollapseIterations genStrictNoErrorsDefaultParsWithHqs(uint32_t stopCheck, uint32_t hqMismatches);
	static CollapseIterations genStrictDefaultPars(uint32_t stopCheck);


	static CollapseIterations genStrictDefaultParsWithHqs(uint32_t stopCheck, uint32_t hqMismatches, bool illumina = false);

	static CollapseIterations genOtuPars(uint32_t stopCheck, double perId, bool onHqPerId = false);

	static CollapseIterations genOtu99To97(uint32_t stopCheck, bool onHqPerId = false);

	static CollapseIterations genOtu99(uint32_t stopCheck, bool onHqPerId = false);

	auto begin() const {
		return iters_.begin();
	}
	auto end() const {
		return iters_.end();
	}

	auto begin() {
		return iters_.begin();
	}
	auto end() {
		return iters_.end();
	}

};




}  // namespace njhseq





