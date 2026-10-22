#pragma once
/*
 * SampleStateStore.hpp
 *
 *  Created on: Oct 22, 2026
 *      Author: nick
 */
// njhseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of njhseq.
//
// njhseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// njhseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with njhseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "njhseq/utils.h"

#include <mutex>
#include <optional>
#include <unordered_map>

namespace njhseq {
namespace collapse {

/**@brief The serialized state of a sample, the same text SampleCollapseCollection::dumpSample writes to samplesOutput/<sample>/, a file that wouldn't be written is left unset
 *
 */
struct SampleState {
	std::optional<std::string> final_; /**< final/<sample>.<ext>, the final clusters*/
	std::optional<std::string> finalInputClusters_; /**< final/clusters/inputClusters.<ext>, the sub clusters of the final clusters*/
	std::optional<std::string> excluded_; /**< excluded/<sample>.<ext>, the excluded clusters*/
	std::optional<std::string> excludedInputClusters_; /**< excluded/clusters/inputClusters.<ext>, the sub clusters of the excluded clusters*/
	std::optional<std::string> refCompInfos_; /**< refCompInfos.tab.txt*/
	std::optional<std::string> refCompInfosExcluded_; /**< refCompInfosExcluded.tab.txt*/

	/**@brief the number of bytes of text held
	 *
	 * @return the size of the text
	 */
	size_t sizeInBytes() const;
};

/**@brief A thread safe in memory store of sample states so the per sample files don't have to be written and re-read between processing steps, each sample is held as one compressed blob
 *
 */
class SampleStateStore {
public:
	/**@brief Set the state for a sample, replacing any previous state
	 *
	 * @param sampleName the sample
	 * @param state the state
	 */
	void set(const std::string & sampleName, SampleState state);

	/**@brief Get the state for a sample
	 *
	 * @param sampleName the sample
	 * @return the state or nullopt if there isn't one stored
	 */
	std::optional<SampleState> get(const std::string & sampleName) const;

	bool has(const std::string & sampleName) const;
	void erase(const std::string & sampleName);

	/**@brief The samples with a stored state
	 *
	 * @return the sample names, sorted
	 */
	VecStr samples() const;

	/**@brief the number of compressed bytes held for all samples
	 *
	 * @return the size held
	 */
	size_t sizeInBytes() const;

	/**@brief the number of bytes of text for all samples before compression
	 *
	 * @return the uncompressed size
	 */
	size_t uncompressedSizeInBytes() const;

private:
	/**@brief a sample's state serialized and compressed
	 *
	 */
	struct CompressedState {
		std::string data_;
		size_t uncompressedSize_{0};
	};
	static CompressedState compress(const SampleState & state);
	static SampleState decompress(const CompressedState & compressed);

	mutable std::mutex mut_;
	std::unordered_map<std::string, CompressedState> states_;
};

}  // namespace collapse
}  // namespace njhseq
