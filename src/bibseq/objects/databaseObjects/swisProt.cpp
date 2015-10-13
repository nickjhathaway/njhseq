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
	if(toks.size() == 4){
		entryName_ = toks[0];
		dataClass_ = rstripReturn(toks[1], ';');
		moleculeType_ = "";
		sequenceLength_ = std::stoi(toks[3]);
	}else if ( toks.size() == 5){
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
		std::stringstream ss;
		ss << "Unrecognized data class" << std::endl;
		ss << dataClass_ << std::endl;
		throw std::runtime_error{ss.str()};
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
