//
// Created by Nicholas Hathaway on 1/28/25.
//

#include "PeptideLibraryReducer.hpp"

namespace njhseq {


std::string PeptideLibraryReducer::reducePeptide(const std::string & peptide) {
	std::string reducedPeptide(peptide.size(), ' ');
	for (const auto pos : iter::range(peptide.size())) {
		reducedPeptide[pos]= static_cast<char>(residueToReduction[peptide[pos]]);
	}

	return reducedPeptide;
}

std::string PeptideLibraryReducer::reverseSimpleFirst(const std::string & reducedPeptide) {
	std::string peptide(reducedPeptide.size(), ' ');
	for (const auto pos : iter::range(reducedPeptide.size())) {
		peptide[pos] = reductionToResidues[reducedPeptide[pos]].front();
	}
	return peptide;
}

void PeptideLibraryReducer::clear() {
	residueToReduction.clear();
	reductionToResidues.clear();
}

void PeptideLibraryReducer::setReductionKeys(const VecStr & clusterings) {
	clear();
	for (const auto & group : iter::enumerate(clusterings)) {
		for (const auto residue : group.element) {
			if (njh::in(residue, residueToReduction)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "already have " << residue << " when processing group: " << group.element << "\n";
				throw std::runtime_error{ss.str()};
			}
			residueToReduction[residue] = static_cast<uint8_t>(group.index + 33);
		}
		reductionToResidues[static_cast<uint8_t>(group.index + 33)] = std::vector<char>(group.element.begin(), group.element.end());
	}
}


void PeptideLibraryReducer::setReduction(const std::string & reduction) {
	if (njh::notIn(reduction, availableReductions_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " no reduction: " << reduction << "available, options are:" << "\n";
		ss << njh::conToStr(availableReductions_, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}
	if (reduction == "UNIPROT18") {
		set_reduction_UNIPROT18();
	} else if (reduction == "HSDM17") {
		set_reduction_HSDM17();
	} else if (reduction == "WASS14") {
		set_reduction_WASS14();
	} else if (reduction == "MMSEQS12") {
		set_reduction_MMSEQS12();
	} else if (reduction == "SDM12") {
		set_reduction_SDM12();
	} else if (reduction == "GBMR7") {
		set_reduction_GBMR7();
	} else if (reduction == "WWMJ5") {
		set_reduction_WWMJ5();
	} else if (reduction == "GBMR4") {
		set_reduction_GBMR4();
	}
}

void PeptideLibraryReducer::set_reduction_UNIPROT18() {
	setReductionKeys(UNIPROT18_);
}
void PeptideLibraryReducer::set_reduction_HSDM17() {
	setReductionKeys(HSDM17_);
}
void PeptideLibraryReducer::set_reduction_MMSEQS12() {
	setReductionKeys(MMSEQS12_);
}
void PeptideLibraryReducer::set_reduction_WASS14() {
	setReductionKeys(WASS14_);
}
void PeptideLibraryReducer::set_reduction_SDM12() {
	setReductionKeys(SDM12_);
}
void PeptideLibraryReducer::set_reduction_GBMR7() {
	setReductionKeys(GBMR7_);
}
void PeptideLibraryReducer::set_reduction_WWMJ5() {
	setReductionKeys(WWMJ5_);
}
void PeptideLibraryReducer::set_reduction_GBMR4() {
	setReductionKeys(GBMR4_);
}


const VecStr PeptideLibraryReducer::availableReductions_ = VecStr{"UNIPROT18", "HSDM17", "WASS14", "MMSEQS12", "SDM12", "GBMR7", "WWMJ5", "GBMR4"};
const VecStr PeptideLibraryReducer::UNIPROT18_ = VecStr{"A", "R", "N", "D", "C", "Q", "EP", "G", "HL", "I", "K", "M", "F", "S", "T", "W", "Y", "V"};
const VecStr PeptideLibraryReducer::HSDM17_ = VecStr{"A", "D", "KE", "R", "N", "T", "S", "Q", "Y", "F", "LIV", "M", "C", "W", "H", "G", "P"};
const VecStr PeptideLibraryReducer::WASS14_ = VecStr{"WM", "DI", "P", "C", "AV", "K", "T", "RE", "G", "L", "Y", "SH", "F", "NQ"};
const VecStr PeptideLibraryReducer::MMSEQS12_ = VecStr{"AST", "LM", "IV", "KR", "EQ", "ND", "FY", "C", "G", "H", "P", "W"};
const VecStr PeptideLibraryReducer::SDM12_ = VecStr{"A", "D", "KER", "N", "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"};
const VecStr PeptideLibraryReducer::GBMR7_ = VecStr{"DN", "AEFIKLMQRVWY", "CH", "T", "S", "G", "P"};
const VecStr PeptideLibraryReducer::WWMJ5_ = VecStr{"CMFILVWY", "ATH", "GP", "DE", "SNQRK"};
const VecStr PeptideLibraryReducer::GBMR4_ = VecStr{"ADKERNTSQ", "YFLIVMCWH", "G", "P"};



}  // namespace njh

