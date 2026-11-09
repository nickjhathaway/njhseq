#pragma once
//
// Adding the variant calls of collapseAndCallVariants to a PMO
//

#include "njhseq/PopulationGeneticsUtils/PopGenUtils.hpp"
#include "njhseq/objects/BioDataObject/pmo/pmo.hpp"
#include "njhseq/programUtils/seqSetUp.hpp"

namespace njhseq {

/**@brief What goes in a PMO's targeted_genomes for the genome the variants were called against
 *
 */
struct PmoGenomePars {
	std::string name_; /**< defaults to the genome file name without its extension*/
	std::string version_; /**< defaults to the name*/
	std::vector<uint32_t> taxonIds_;
	std::string url_{"NA"};
	std::string gffUrl_;

	void setOptions(seqSetUp & setUp);

	/**@brief Generate the genome info
	 *
	 * @param genomeFnp the genome file, used for the name when none was given
	 * @return the genome info
	 */
	[[nodiscard]] pmo::GenomeInfo genGenomeInfo(const bfs::path & genomeFnp) const;
};

/**@brief Counts of what addVariantCallsToPmo added
 *
 */
struct AddVariantCallsToPmoResults {
	uint32_t mhapsWithCalls_{0};
	uint32_t mhapsWithoutCalls_{0}; /**< haplotypes that weren't in the variant calling or couldn't be mapped to the genome*/
	uint32_t proteinVariants_{0};
	uint32_t seqVariants_{0};
};

/**@brief Add the variant calls from collapseAndCallVariants to the representative microhaplotypes of a target and the genome to targeted_genomes
 *
 * Amino acid positions (protein_location) are 1-based with the end included, so start and end are the same for one amino acid,
 * genomic positions (codon_genomic_location and the sequence variants) are 0-based with the end not included.
 * Protein variants are the amino acid of the haplotype at every typed position, known mutations and positions variable in the population, including ones matching the reference, the same positions as seqsAATyped.tab.txt.gz.
 * Each has the gene id, the gene's Name from the gff as alternative_gene_name when it has one (characters the schema doesn't allow become _), and the codon's genomic location, with the codons on the gene's strand.
 * Sequence variants are the haplotype's SNPs and indels against the genome, the same as variantsPerSeqAln.tab.txt.gz, on the + strand.
 * Haplotypes are matched by sequence to the uniqueSeqs.fasta.gz collapseAndCallVariants wrote, and any variants they already have are replaced.
 *
 * @param pmoObj the PMO to update
 * @param targetName the target the haplotypes are for
 * @param varCallPars the parameters collapseAndCallVariants was run with
 * @param calls what collapseAndCallVariants returned
 * @param genomePars the genome info
 * @return counts of what was added
 */
AddVariantCallsToPmoResults addVariantCallsToPmo(pmo::PortableMicrohaplotypeObject & pmoObj,
		const std::string & targetName,
		const CollapseAndCallVariantsPars & varCallPars,
		const TranslatorByAlignment::TranslatorByAlignmentResult & calls,
		const PmoGenomePars & genomePars);

}  // namespace njhseq
