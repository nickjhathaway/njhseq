#pragma once

//
// Created by Nicholas Hathaway on 1/28/25.
//


#include "njhseq/utils.h"

namespace njhseq {

class PeptideLibraryReducer {

	// https://pmc.ncbi.nlm.nih.gov/articles/PMC10872054/ table 1
	/*
	Alphabet=Clustering
		UNIPROT20=A R N D C Q E G H I L K M F P S T W Y V
		UNIPROT18=A R N D C Q EP G HL I K M F S T W Y V
		HSDM17=A D KE R N T S Q Y F LIV M C W H G P
		MMSEQS12=AST LM IV KR EQ ND FY C G H P W
		WASS14=WM DI P C AV K T RE G L Y SH F NQ
		SDM12=A D KER N TSQ YF LIVM C W H G P
		GBMR7=DN AEFIKLMQRVWY CH T S G P
		WWMJ5=CMFILVWY ATH GP DE SNQRK
		GBMR4=ADKERNTSQ YFLIVMCWH G P
	 *
	 */
public:

	PeptideLibraryReducer() = default;


	/**
	 * @brief Reduce a peptide to a coded reduced library, each residue will be recoded to a char corresponding to its numeric value in the ASCII
	 * @param peptide the peptide to reduce
	 * @return a reduced string encoded with the reduction, convert into real amio acid again with reverseSimpleFirst but cannot recover what the original amion acids were
	 */
	std::string reducePeptide(const std::string & peptide);


	/**
	 * @brief Convert a reduced encode peptide into real amino acids, cannot get back original sequence
	 * @param reducedPeptide the reduced encoded peptide to convert back into amino acids, will default to the first amino acid listed for the reduction symbol
	 * @return amino acids decoded from its reduced library symbol
	 */
	std::string reverseSimpleFirst(const std::string & reducedPeptide);


	std::unordered_map<char, uint8_t> residueToReduction;
	std::unordered_map<uint8_t, std::vector<char>> reductionToResidues;

	void clear();


	/**
	 * @brief Create reduction library by giving grouped amino acids in strings, checks to see if giving the same amino acids into multiple groups
	 * @param clustering a vector of strings with the different groupings, each word will be grouped together, e.g. giving VecStr{"ADKERNTSQ","YFLIVMCWH","G,"P"} will create 4 groups
	 */
	void setReductionKeys(const VecStr & clustering);

	// pre-built reduction libraries
	const static VecStr UNIPROT18_;
	const static VecStr HSDM17_;
	const static VecStr WASS14_;
	const static VecStr MMSEQS12_;
	const static VecStr SDM12_;
	const static VecStr GBMR7_;
	const static VecStr WWMJ5_;
	const static VecStr GBMR4_;

	void set_reduction_UNIPROT18();
	void set_reduction_HSDM17();
	void set_reduction_MMSEQS12();
	void set_reduction_WASS14();
	void set_reduction_SDM12();
	void set_reduction_GBMR7();
	void set_reduction_WWMJ5();
	void set_reduction_GBMR4();


	/**
	 * @brief set the reduction library by the name of a prebuilt reduction
	 * @param reduction the prebuilt reduction library to use
	 */
	void setReduction(const std::string & reduction);

	const static VecStr availableReductions_;
};



} //namespace njhseq
