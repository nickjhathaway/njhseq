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

/*
 * mip.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: nickhathaway
 */

#include "mip.hpp"

namespace bibseq {
void mip::setArms(){
	//forward_ = forwardPrimer_ + std::string(mbLen_, 'N') + seqUtil::reverseComplement(extentionArm_, "DNA");
	forward_ = forwardPrimer_ + std::string(mbLen_, 'N') + extentionArm_;

	forwardInfo_ = readObject(seqInfo("forward", forward_));
	forwardInfoRC_ = readObject(seqInfo("forwardRc", seqUtil::reverseComplement(forward_, "DNA")));

	reverse_ = reversePrimer_ + std::string(mbLen_, 'N') + seqUtil::reverseComplement(ligationArm_, "DNA");

	ligationArmRc_ =  seqUtil::reverseComplement(ligationArm_, "DNA");
	extentionArmRc_ =  seqUtil::reverseComplement(extentionArm_, "DNA");
	reverseInfo_ = readObject(seqInfo("reverse", reverse_));
	reverseInfoRC_ = readObject(seqInfo("reverseRc", seqUtil::reverseComplement(reverse_, "DNA")));

	forwardPrimerObj_ = readObject(seqInfo("forwardPrimer", forwardPrimer_));
	reversePrimerObj_ = readObject(seqInfo("reversePrimer", reversePrimer_));
	forwardPrimerObjRc_ = readObject(seqInfo("forwardPrimerRc", seqUtil::reverseComplement(forwardPrimer_, "DNA")));
	reversePrimerObjRc_ = readObject(seqInfo("reversePrimerRc", seqUtil::reverseComplement(reversePrimer_, "DNA")));

	extentionArmObj_ = readObject(seqInfo("extentionArm", extentionArm_));
	ligationArmObj_ = readObject(seqInfo("ligationArm", ligationArm_));
	extentionArmObjRc_ = readObject(seqInfo("extentionArmRc", seqUtil::reverseComplement(extentionArm_, "DNA")));
	ligationArmObjRc_ = readObject(seqInfo("ligationArmRc", seqUtil::reverseComplement(ligationArm_, "DNA")));
}
void mip::printDescription(std::ostream & out, bool deep){
	std::cout << "mip{" << std::endl;

	std::cout << "}" << std::endl;

}
} /* namespace bib */
