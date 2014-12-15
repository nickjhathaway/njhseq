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
 * swisProt.cpp
 *
 *  Created on: Apr 8, 2014
 *      Author: nickhathaway
 */

#include "swisProt.hpp"

namespace bibseq {
void swisProt::readIdLine(const std::string & line){
	/*
	 * Prior to release 51, included with MoleculeType:
	 * ID   EntryName DataClass; MoleculeType; SequenceLength AA.
	 *
	 * Newer files lack the MoleculeType:
	 * ID   EntryName DataClass; SequenceLength AA.
	 */
	auto toks = tokenizeString(line, "whitespace");
	if(len(toks) == 4){
		entryName_ = toks[0];
		dataClass_ = rstripReturn(toks[1], ';');
		moleculeType_ = "";
		sequenceLength_ = std::stoi(toks[3]);
	}else if ( len(toks) == 5){
		entryName_ = toks[0];
		dataClass_ = rstripReturn(toks[1], ';');
		moleculeType_ = rstripReturn(toks[2], ';');
		sequenceLength_ = std::stoi(toks[3]);
	}else{
		std::cout << "ID line has unrecognised format: " << std::endl;
		std::cout << line << std::endl;
	}
	VecStr dataClassPossible = VecStr {"STANDARD", "PRELIMINARY", "IPI", "Reviewed", "Unreviewed"};

	if(in(dataClass_, dataClassPossible)){
		std::cout << "Unrecognized data class" << std::endl;
		std::cout << dataClass_ << std::endl;
		exit(1);
	}
	if("" != moleculeType_ && "PRT" != moleculeType_){
		std::cout << "Unrecognized molecular type, should be PRT or blank" << std::endl;
		std::cout << moleculeType_ << std::endl;
	}
}
void swisProt::printDescription(std::ostream & out, bool deep){
	out << "swisProt{" << std::endl;
	out << "entryName_: " << entryName_ << std::endl;
	out << "dataClass_: " << dataClass_ << std::endl;
	out << "moleculeType_: " << moleculeType_ << std::endl;
	out << "sequenceLength_: " << sequenceLength_ << std::endl;
	out << "}" << std::endl;

}
} /* namespace bib */
