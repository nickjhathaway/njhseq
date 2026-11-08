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
 * add throws if the same name has conflicting information, if a library sample has results for the same target more than once,
 * or if a PMO has project, sequencing, bioinformatics, genome or read count sections, which is what makes them more than minimum PMOs.
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

}  // namespace njhseq::pmo
