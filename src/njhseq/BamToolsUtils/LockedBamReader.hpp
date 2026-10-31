#pragma once

/*
 * LockedBamReader.hpp
 *
 *  Created on: Oct 31, 2026
 *      Author: nick
 */


#include "njhseq/BamToolsUtils/BamToolsUtils.hpp"

namespace njhseq {

/**
 * @brief A BamTools::BamReader that several threads can share, each taking batches of alignments in file order under a lock
 */
class LockedBamReader {
public:
	/**
	 * @brief Open a bam file, throwing if it can't be opened
	 * @param bamFnp the bam file
	 */
	explicit LockedBamReader(const bfs::path & bamFnp);

	/**
	 * @brief Replace alns with the next batch of alignments
	 * @param alns cleared then filled with up to batchSize alignments
	 * @param batchSize the most alignments to read
	 * @return false once the file has no more alignments
	 */
	bool readNextAlnBatch(std::vector<BamTools::BamAlignment> & alns, uint32_t batchSize);

private:
	BamTools::BamReader bReader_;
	std::mutex mut_;
};

}  // namespace njhseq
