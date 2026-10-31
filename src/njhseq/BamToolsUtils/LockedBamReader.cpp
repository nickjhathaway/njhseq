/*
 * LockedBamReader.cpp
 *
 *  Created on: Oct 31, 2026
 *      Author: nick
 */


#include "LockedBamReader.hpp"

namespace njhseq {

LockedBamReader::LockedBamReader(const bfs::path & bamFnp) {
	bReader_.Open(bamFnp.string());
	checkBamOpenThrow(bReader_, bamFnp);
}

bool LockedBamReader::readNextAlnBatch(std::vector<BamTools::BamAlignment> & alns, uint32_t batchSize) {
	std::lock_guard<std::mutex> lock(mut_);
	alns.clear();
	BamTools::BamAlignment bAln;
	while (alns.size() < batchSize && bReader_.GetNextAlignment(bAln)) {
		alns.emplace_back(bAln);
	}
	return !alns.empty();
}

}  // namespace njhseq
