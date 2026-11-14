#pragma once
//
// Operations on whole PMOs, pmo.hpp has the classes for the format itself
//

#include "njhseq/objects/BioDataObject/pmo/pmo.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/IO/OutputStream.hpp"

namespace njhseq::pmo {

/**@brief Read and validate a PMO file
 *
 * @param fnp the file, can be gzipped
 * @return the PMO
 */
[[nodiscard]] PortableMicrohaplotypeObject readPmo(const bfs::path & fnp);

/**@brief Write a PMO
 *
 * @param pmoObj the PMO to write
 * @param outOpts where to write, gzipped if it ends in .gz
 */
void writePmo(const PortableMicrohaplotypeObject & pmoObj, const OutOptions & outOpts);

/**@brief Merge minimum PMOs one at a time so they don't all have to be in memory, targets, panels, specimens and library samples with the same name become one and everything else is added
 *
 * This is for combining PMOs per target over the same samples as well as PMOs over different samples.
 * Genomes with the same name and version become one, and a haplotype in more than one PMO keeps the variants from all of them.
 * add throws if the same name has conflicting information, if a library sample has results for the same target more than once,
 * or if a PMO has project, sequencing, bioinformatics or read count sections.
 */
class MinimumPmoMerger {
public:
	/**
	 * @param generationMethod the program doing the merging
	 */
	explicit MinimumPmoMerger(const PmoGenerationMethod & generationMethod);

	void add(const PortableMicrohaplotypeObject & pmoObj);

	/**@brief Get the merged PMO, the merger shouldn't be used after this
	 *
	 * @return the validated merged PMO
	 */
	[[nodiscard]] PortableMicrohaplotypeObject finish();

private:
	PortableMicrohaplotypeObject merged_;
	DetectedMicrohaplotypes detected_;
	std::unordered_map<std::string, uint32_t> genomeIds_;
	std::unordered_map<std::string, uint32_t> targetIds_;
	std::unordered_map<std::string, uint32_t> panelIds_;
	std::unordered_map<std::string, uint32_t> specimenIds_;
	std::unordered_map<std::string, uint32_t> libraryIds_;
	std::unordered_map<uint32_t, uint32_t> repPosByTargetId_;
	std::unordered_map<uint32_t, std::unordered_map<std::string, uint32_t>> mhapIdBySeq_;
	std::unordered_map<uint32_t, std::unordered_set<std::string>> mhapNames_;
	std::unordered_map<uint32_t, uint32_t> detectedPosByLibraryId_;
	std::unordered_map<uint32_t, std::set<uint32_t>> repPosDoneByLibraryId_;
};

/**@brief Merge minimum PMOs, see MinimumPmoMerger
 *
 * @param pmos the PMOs to merge
 * @param generationMethod the program doing the merging
 * @return the validated merged PMO
 */
[[nodiscard]] PortableMicrohaplotypeObject mergeMinimumPmos(const std::vector<PortableMicrohaplotypeObject> & pmos,
		const PmoGenerationMethod & generationMethod);

/**@brief Results of adding meta to a PMO
 *
 */
struct AddSpecimenMetaResults {
	uint32_t specimensUpdated_{0};
	VecStr metaSamplesNotInPmo_;
};

/**@brief Add columns of a table as fields of specimens, columns that aren't specimen fields of the schema become custom fields
 *
 * Rows are matched on library sample name and then specimen name, empty and NA values are skipped.
 *
 * @param pmoObj the PMO to update
 * @param meta the table to add
 * @param sampleCol the column with the sample names
 * @return how many were updated and any samples in the table not in the PMO
 */
AddSpecimenMetaResults addSpecimenMeta(PortableMicrohaplotypeObject & pmoObj, const table & meta, const std::string & sampleCol);

/**@brief Where to find the target locations, genes and reference sequences for addTargetGenomeInfo
 *
 */
struct AddTargetGenomeInfoPars {
	bfs::path locationsDir_; /**< the locationsByGenome directory of SeekDeep genTargetInfoFromGenomes, its <genome>_infos.tab.txt tables are read*/
	bfs::path gffDir_; /**< optional, <genome>.gff for the genes each insert overlaps, genomes without one get no gene_name*/
	bfs::path genomeDir_; /**< optional, <genome>.2bit for the reference sequence of each insert*/
	VecStr preferredGenomes_; /**< a target found in several genomes goes to the first of these it's in*/
};

/**@brief What addTargetGenomeInfo added
 *
 */
struct AddTargetGenomeInfoResults {
	std::map<std::string, uint32_t> targetsPerGenome_;
	uint32_t targetsWithGenes_{0};
	VecStr targetsNotFound_;
};

/**@brief Add each target's genome, insert and primer locations, and the genes the insert overlaps
 *
 * A target found in more than one genome goes to the first of preferredGenomes_ it's in, then the genome with the fewest primer errors,
 * then a genome whose name starts the target name, then by genome name. Within a genome the hit with the fewest primer errors is used.
 * The genome is added to targeted_genomes with its name as the name and version, reusing an entry with the same name.
 * Targets clustered as separate mates (<target>-R1, <target>-R2) get the location of <target>.
 * gene_name is GENE_ID:GENE_NAME of each gene in the gff the insert overlaps, comma separated, with NA for genes without a Name.
 * Locations are 0-based with the end not included and on the target's strand, the reverse primer on the opposite strand.
 *
 * @param pmoObj the PMO to update
 * @param pars where the target info is
 * @return the targets added per genome and the targets that weren't found
 */
AddTargetGenomeInfoResults addTargetGenomeInfo(PortableMicrohaplotypeObject & pmoObj, const AddTargetGenomeInfoPars & pars);

}  // namespace njhseq::pmo
