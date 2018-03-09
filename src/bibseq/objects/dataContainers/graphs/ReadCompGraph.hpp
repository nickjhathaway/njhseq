#pragma once
/*
 * ReadCompGraph.hpp
 *
 *  Created on: Dec 1, 2015
 *      Author: nick
 */

//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of bibseq.
//
// bibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// bibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bibseq.  If not, see <http://www.gnu.org/licenses/>.
//
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "bibseq/objects/dataContainers/graphs/readDistGraph.hpp"
#include "bibseq/alignment/alignerUtils/comparison.hpp"
#include "bibseq/alignment/aligner.h"

#include "bibseq/seqToolsUtils/RefDetermination/BestRefDetector.hpp"
#include "bibseq/objects/Meta/MultipleGroupMetaData.hpp"
#include "bibseq/concurrency/pools/AlignerPool.hpp"



namespace bibseq {

class ReadCompGraph: public njhUndirWeightedGraph<comparison,
		std::shared_ptr<seqInfo>> {
public:
	/**@brief Construct with a distance matrix and a vector of reads that were used to create the distances
	 *
	 * @param distances The distance matrix, the matrix should at least have the diagonal values (each row has as many elements as it's row position
	 * @param reads the reads the distance graph is describing
	 *
	 * @todo Add safety checks for size of distances input matrix
	 */
	template<typename T>
	ReadCompGraph(const std::vector<std::vector<comparison>> & distances,
			const std::vector<T> & reads) {
		for (const auto & pos : iter::range(reads.size())) {
			this->addNode(getSeqBase(reads[pos]).name_,
					std::make_shared<seqInfo>(getSeqBase(reads[pos])));
		}
		for (const auto & pos : iter::range(distances.size())) {
			for (const auto & subPos : iter::range<uint64_t>(distances[pos].size())) {
				this->addEdge(getSeqBase(reads[pos]).name_,
						getSeqBase(reads[subPos]).name_, distances[pos][subPos]);
			}
		}
	}

	/**@brief Construct with just reads for the edges to be added latter
	 *
	 * @param reads the sequences to add, a copy of the seqBases will be copied and used
	 */
	template<typename T>
	ReadCompGraph(const std::vector<T> & reads) {
		for (const auto & pos : iter::range(reads.size())) {
			this->addNode(getSeqBase(reads[pos]).name_,
					std::make_shared<seqInfo>(getSeqBase(reads[pos])));
		}
	}

	std::map<uint32_t, std::vector<char>> getVariantSnpLociMap(
			const std::string & name, VecStr names, uint32_t expand = 0) const;
	std::map<uint32_t, std::vector<gap>> getVariantIndelLociMap(
			const std::string & name, VecStr names, uint32_t expand = 0) const;

	comparison setMinimumConnections(std::function<void(comparison &)> modFunc,
			std::function<bool(const comparison &, const comparison &)> compFunc);
	comparison setMinimumEventConnections();
	comparison setMinimumHqMismatchConnections();

	void setJustBestConnection(bool doTies);

	Json::Value toD3Json(bib::color backgroundColor,
			const std::unordered_map<std::string, bib::color> & nameColors);




	struct ConnectedHaplotypeNetworkPars {
		BestRefDetector::FindBestRefPars matchPars;
		uint32_t numThreads = 1;
		uint32_t minNumberOfEvents = std::numeric_limits<uint32_t>::max();
		double minId = std::numeric_limits<double>::lowest();
		bool setJustBest = false;
		bool doTies = false;
		std::string colorField = "";
		std::string labelField = "";
		bool noLabel = false;
		std::unordered_map<std::string, bib::color> colorLookup;
		std::unique_ptr<MultipleGroupMetaData> seqMeta;
		uint32_t minlength = 50;
		uint32_t kLen = 47;
		bool verbose = false;
		void loadMeta(const bfs::path & metaFnp){
			seqMeta = std::make_unique<MultipleGroupMetaData>(metaFnp);
		}

		void checkForColorField(const std::string & funcName){
			if("" != colorField){
				if(nullptr == seqMeta){
					std::stringstream ss;
					ss << funcName << ", supplied a color field but no meta data" << "\n";
					throw std::runtime_error{ss.str()};
				}
				seqMeta->checkForFieldsThrow({colorField});
			}
		}

		void checkForLabelField(const std::string & funcName){
			if("" != labelField){
				if(nullptr == seqMeta){
					std::stringstream ss;
					ss << funcName << ", supplied a color field but no meta data" << "\n";
					throw std::runtime_error{ss.str()};
				}
				seqMeta->checkForFieldsThrow({labelField});
			}
		}



		std::unordered_map<std::string, bib::color> generateColorLookup(const bfs::path & fieldColorFnp){
			table colorTab;
			std::unordered_map<std::string, bib::color> colorLookup;
			if("" != fieldColorFnp){
				colorTab  = table(fieldColorFnp, "\t", true);
				if("" != colorField){
					colorTab.checkForColumnsThrow(VecStr{colorField, "color"}, __PRETTY_FUNCTION__);
					for(const auto & row : colorTab.content_){
						colorLookup[row[colorTab.getColPos(colorField)]] = bib::color(row[colorTab.getColPos("color")]);
					}
				}
			} else if("" != colorField){
				colorLookup = getColorsForNames(getVectorOfMapKeys(seqMeta->groupData_.at(colorField)->subGroupToSamples_));
			}
			return colorLookup;
		}

		static std::string htmlPageForConnectHpaNet;
	};

	Json::Value getSingleLineJsonOut(const ConnectedHaplotypeNetworkPars & netPars,
			const std::unordered_map<std::string, bib::color> & colorLookup);

	void writeAdjListPerId(const OutOptions & linksOutOpts,
			const ConnectedHaplotypeNetworkPars & netPars);

	void addEdgesBasedOnIdOrMinDif(const ConnectedHaplotypeNetworkPars & netPars,
			concurrent::AlignerPool & alnPool);

};

template<typename T>
ReadCompGraph genReadComparisonGraph(const std::vector<T> & reads,
		aligner & alignerObj,
		std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
		std::mutex & alignerLock, uint32_t numThreads) {
	std::function<
			comparison(const T &, const T &,
					std::unordered_map<std::string, std::unique_ptr<aligner>>&, aligner&)> getMismatchesFunc =
			[&alignerLock](const T & read1, const T & read2,
					std::unordered_map<std::string, std::unique_ptr<aligner>>& aligners,
					aligner &alignerObj) {
				alignerLock.lock();
				auto threadId = estd::to_string(std::this_thread::get_id());
				//std::cout << threadId<< std::endl;
				if(aligners.find(threadId) == aligners.end()) {
					aligners.emplace(threadId, std::make_unique<aligner>(alignerObj));
				}
				alignerLock.unlock();
				aligners.at(threadId)->alignCache(getSeqBase(read1),getSeqBase(read2), false);
				aligners.at(threadId)->profilePrimerAlignment(getSeqBase(read1), getSeqBase(read2));
				return aligners.at(threadId)->comp_;
			};
	auto distances = getDistanceNonConst(reads, numThreads, getMismatchesFunc,
			aligners, alignerObj);
	return ReadCompGraph(distances, reads);
}

template<typename T>
ReadCompGraph genReadComparisonGraph(const std::vector<T> & reads,
		aligner & alignerObj, uint32_t numThreads) {
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	std::mutex alignerLock;
	return genReadComparisonGraph(reads, alignerObj, aligners, alignerLock,
			numThreads);
}

}  // namespace bibseq

