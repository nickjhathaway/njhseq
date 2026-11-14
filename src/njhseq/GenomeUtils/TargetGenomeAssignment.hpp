#pragma once
//
// Choosing the one genome each target is placed in when SeekDeep genTargetInfoFromGenomes found it in several genomes
//

#include "njhseq/common.h"
#include "njhseq/objects/BioDataObject/BedRecordCore.hpp"

namespace njhseq {

/**@brief Where a target was found in a genome, a row of genTargetInfoFromGenomes' <genome>_infos.tab.txt
 *
 * Positions are 0-based with the end not included.
 */
struct TargetGenomeHit {
	std::string genome_;
	std::string chrom_;
	char strand_{'+'};
	uint32_t insertStart_{0};
	uint32_t insertEnd_{0};
	uint32_t fwdStart_{0};
	uint32_t fwdEnd_{0};
	uint32_t revStart_{0};
	uint32_t revEnd_{0};
	uint32_t primerErrors_{0};
};

/**@brief The target locations of a genTargetInfoFromGenomes locationsByGenome directory
 *
 */
class TargetGenomeLocations {
public:
	/**@brief Read the <genome>_infos.tab.txt and <genome>_inner.bed files of a locationsByGenome directory
	 *
	 * @param locationsDir the directory
	 * @return the locations
	 */
	static TargetGenomeLocations read(const bfs::path & locationsDir);

	std::unordered_map<std::string, std::vector<TargetGenomeHit>> hitsByTarget_;
	/**< key = genome, the regions of <genome>_inner.bed in file order, named by target, mates of targets analyzed as separate R1 and R2 targets are NAME-R1 and NAME-R2*/
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> innerRegionsByGenome_;

	/**@brief The target a name is for, the name itself or for NAME-R1 and NAME-R2 that aren't targets themselves NAME
	 *
	 * @param name the target or mate name
	 * @return the target's name
	 */
	[[nodiscard]] std::string baseTarget(const std::string & name) const;

	/**@brief The inner region of a target, or of one of its mates, at a hit
	 *
	 * @param name the target or mate name
	 * @param hit the target's hit
	 * @return the region named name in the hit's genome on the hit's chromosome overlapping the hit's insert, null if there isn't one
	 */
	[[nodiscard]] std::shared_ptr<Bed6RecordCore> innerRegion(const std::string & name, const TargetGenomeHit & hit) const;
};

/**@brief How a target's genome is chosen
 *
 */
struct TargetGenomeChoicePars {
	VecStr preferredGenomes_; /**< a target goes to the first of these it's in*/
	bool onlyPreferred_{false}; /**< targets not in any of preferredGenomes_ get no genome rather than the best of the others*/
};

/**@brief Choose the genome of each target
 *
 * A target goes to the first of the preferred genomes it's in, then to the genome with the fewest primer errors,
 * then to a genome whose name starts the target name, then by genome name. Within a genome the hit with the fewest primer errors is used.
 *
 * @param locations the target locations
 * @param pars the preferred genomes
 * @return the chosen hit of each target, targets that weren't found (or with onlyPreferred_ weren't in a preferred genome) are left out
 */
std::map<std::string, TargetGenomeHit> chooseTargetGenomes(const TargetGenomeLocations & locations, const TargetGenomeChoicePars & pars);

}  // namespace njhseq
