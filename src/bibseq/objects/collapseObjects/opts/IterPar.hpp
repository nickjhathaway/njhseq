#pragma once
//
//  iterPars.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/13/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
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
#include "bibseq/utils.h"
#include "bibseq/alignment/alignerUtils/comparison.hpp"


namespace bibseq {

struct CompareProfile {
	virtual bool passErrors(const comparison & threshold,
			const comparison & generated) = 0;
	virtual ~CompareProfile() {
	}
};

struct CompareIDProfile: public CompareProfile {
	virtual bool passErrors(const comparison & threshold,
			const comparison & generated);
	virtual ~CompareIDProfile();
};

struct CompareIDHqProfile: public CompareProfile {
	virtual bool passErrors(const comparison & threshold,
			const comparison & generated);
	virtual ~CompareIDHqProfile();
};

struct CompareErrorProfile: public CompareProfile {
	virtual bool passErrors(const comparison & threshold,
			const comparison & generated);
	virtual ~CompareErrorProfile();
};

struct CompareIDAndErrorProfile: public CompareProfile {
	virtual bool passErrors(const comparison & threshold,
			const comparison & generated);
	virtual ~CompareIDAndErrorProfile();
};


class IterPar {
public:
	struct PerIdPars{
		bool onPerId_; //should be on if doing percent identity OTU clustering
		bool onHqPerId_;// if this is on onPerId_ should be on as well, this being on indicates to do only high quality events
	};
	IterPar();
	IterPar(const std::vector<double>& parameter, uint32_t iterNumber,
			const PerIdPars & perIdPars);

	IterPar& operator=(const IterPar & other);
	IterPar& operator=(IterPar&& other);
	IterPar(const IterPar& other);

	// procedure parameters
	uint32_t stopCheck_;
	uint32_t smallCheckStop_;
	uint32_t iterNumber_;
	bool onPerId_ = false;
	bool onHqPerId_ = false;
	// error parameters
	comparison errors_;
private:
	std::shared_ptr<CompareProfile> comp_;
public:
	bool passErrorCheck(const comparison & errors) const;

	void setCompFunc();

	void printIterInfo(std::ostream & out, bool colorFormat) const;

	VecStr outputIterInfo() const;
	static VecStr outputIterInfoHeader(bool onPerId);

};




}  // namespace bib


