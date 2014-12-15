//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#pragma once
/*

 * mip.h
 *
 *  Created on: Mar 31, 2014
 *      Author: nickhathaway
 */
#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/readObject.hpp"

namespace bibseq {

class mip {
public:
	mip(){}
	mip(const std::string & forwardPrimer,
			const std::string & reversePrimer,
			uint32_t mbLen,
			const std::string & ligationArm,
			const std::string & extentionArm,
			const std::string & name): forwardPrimer_(forwardPrimer), reversePrimer_(reversePrimer),
			mbLen_(mbLen), ligationArm_(ligationArm), extentionArm_(extentionArm), name_(name){
		setArms();
	}
	//members
	std::string forwardPrimer_;
	std::string reversePrimer_;
	readObject forwardPrimerObj_;
	readObject reversePrimerObj_;
	readObject forwardPrimerObjRc_;
	readObject reversePrimerObjRc_;
	uint32_t mbLen_;
	std::string ligationArm_;
	std::string extentionArm_;
	std::string ligationArmRc_;
	std::string extentionArmRc_;

	std::string name_;
	readObject ligationArmObj_;

	readObject extentionArmObj_;
	readObject ligationArmObjRc_;

	readObject extentionArmObjRc_;
	std::string forward_;
	readObject forwardInfo_;
	readObject forwardInfoRC_;
	std::string reverse_;
	readObject reverseInfo_;
	readObject reverseInfoRC_;
	void setArms();
	virtual void printDescription(std::ostream & out, bool deep);
};

} /* namespace bib */

#ifndef NOT_HEADER_ONLY
#include "mip.cpp"
#endif
