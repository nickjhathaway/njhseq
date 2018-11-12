#pragma once
/*
 * GenomicAminoAcidPositionTyper.hpp
 *
 *  Created on: Sep 29, 2018
 *      Author: nick
 */



#include "bibseq/objects/Gene/GeneSeqInfo.hpp"
#include "bibseq/objects/Gene/GeneFromGffs.hpp"



namespace bibseq {

class GenomicAminoAcidPositionTyper{
public:
	struct GeneAminoTyperInfo {
		GeneAminoTyperInfo(const std::string & geneId);
		/**@brief positions should be 0-based
		 *
		 * @param geneId
		 * @param aminos
		 */
		GeneAminoTyperInfo(const std::string & geneId,
				const std::map<uint32_t, char> & aminos);
		std::string geneId_;
		std::string altId_;
		std::map<uint32_t, char> aminos_; //positions are 0-based
	};

	GenomicAminoAcidPositionTyper(const bfs::path & proteinMutantTypingFnp,
			bool inputZeroBased = false);

	bfs::path proteinMutantTypingFnp_;
	bool inputZeroBased_;

	table proteinMutantTypingTab_;
	std::unordered_map<std::string, std::vector<uint32_t>> aminoPositionsForTyping_;
	std::unordered_map<std::string, std::set<std::string>> altNamesForIds_;
	std::unordered_map<std::string, GeneAminoTyperInfo> aminoPositionsForTypingWithInfo_;


	std::set<std::string> getGeneIds() const;

	/**@brief Type the amino acids for the alignment
	 *
	 * @param bAln the alignment to type
	 * @param currentGene the gene that this alignment intersects
	 * @param currentGeneInfo the base information for this gene
	 * @param tReader a twobit reader
	 * @param alignerObj aligner
	 * @param refData the chromosome index to chrom name from bam reader
	 * @return a std::map with key zero based position in protein and value is amino acid for that position
	 */
	std::map<uint32_t, char> typeAlignment(
			const BamTools::BamAlignment & bAln,
			const GeneFromGffs & currentGene,
			const GeneSeqInfo & currentGeneInfo,
			TwoBit::TwoBitFile & tReader,
			aligner & alignerObj,
			const BamTools::RefVector & refData);
};





}  // namespace bibseq


