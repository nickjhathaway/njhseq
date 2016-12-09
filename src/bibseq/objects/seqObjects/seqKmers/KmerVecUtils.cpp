/*
 * KmerVecUtils.cpp
 *
 *  Created on: May 24, 2016
 *      Author: nick
 */

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


#include "KmerVecUtils.hpp"
#include "bibseq/seqToolsUtils/distCalc.hpp"
#include "bibseq/IO/SeqIO/SeqInput.hpp"

namespace bibseq {


void allSetKmers(std::vector<std::unique_ptr<seqWithKmerInfo>> & reads, uint32_t kLength, bool setReverse){
	bib::for_each(reads, [&](std::unique_ptr<seqWithKmerInfo> & read){ read->setKmers(kLength, setReverse);});
}

void allSetKmers(std::vector<seqWithKmerInfo> & reads, uint32_t kLength, bool setReverse){
	bib::for_each(reads, [&](seqWithKmerInfo & read){ read.setKmers(kLength, setReverse);});
}


std::vector<std::vector<double>> getKmerAccerDistance(
		std::vector<std::unique_ptr<seqWithKmerInfo>>& reads, uint32_t kmerStart,
		uint32_t kmerStop, uint32_t numThreads, bool useKNumber, bool verbose) {
	//kmerStop is inclusive
	if(kmerStart > kmerStop){
		std::stringstream ss;
		ss << "getKmerAccerDistance: kmerStop is less than kmerStart, kmerStart: " << kmerStart <<
				", KmerStop: " << kmerStop;
		throw std::runtime_error{ss.str()};
	}
	std::vector<std::vector<double>> distances;
	if (useKNumber) {
		std::function<
				uint32_t(const std::unique_ptr<seqWithKmerInfo> &,
						const std::unique_ptr<seqWithKmerInfo> &)> disFun =
				[](const std::unique_ptr<seqWithKmerInfo> & read1,
						const std::unique_ptr<seqWithKmerInfo> & read2) {
					auto dist = read1->compareKmers(*read2);
					return dist.first;
				};
		std::unordered_map<uint32_t, std::vector<std::vector<uint32_t>>>distanceMaps;
		for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
			if(verbose){
				std::cout << "K: " << k << std::endl;
				std::cout << "\tIndexing Kmers" << std::endl;
			}
			bib::stopWatch watch;
			allSetKmers(reads, k, false);

			if(verbose){
				std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
				std::cout << "\tCalculating Distances" << std::endl;
			}
			watch.reset();
			distanceMaps[k] = getDistance(reads, numThreads, disFun);
			if(verbose){
				std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
			}
		}
		for (const auto & rowPos : iter::range(distanceMaps[kmerStart].size())) {
			std::vector<double> temp;
			for (uint32_t i = 0; i < distanceMaps[kmerStart][rowPos].size(); ++i) {
				temp.emplace_back(0.00);
			}
			distances.emplace_back(temp);
		}
		for (const auto & rowPos : iter::range(distances.size())) {
			for (const auto & colPos : iter::range(distances[rowPos].size())) {
				std::vector<double> differences;
				for (uint32_t k = kmerStart; k < kmerStop; ++k) {
					differences.emplace_back(
							uAbsdiff(distanceMaps[k][rowPos][colPos],
									distanceMaps[k + 1][rowPos][colPos]));
				}
				distances[rowPos][colPos] = vectorMean(differences);
			}
		}
	} else {
		std::function<
				double(const std::unique_ptr<seqWithKmerInfo> &,
						const std::unique_ptr<seqWithKmerInfo> &)> disFun =
				[](const std::unique_ptr<seqWithKmerInfo> & read1,
						const std::unique_ptr<seqWithKmerInfo> & read2) {
					auto dist = read1->compareKmers(*read2);
					return dist.second;
				};
		std::unordered_map<uint32_t, std::vector<std::vector<double>>>distanceMaps;
		for (uint32_t k = kmerStart; k < kmerStop + 1; ++k) {
			if(verbose){
				std::cout << "K: " << k << std::endl;
				std::cout << "\tIndexing Kmers" << std::endl;
			}
			bib::stopWatch watch;
			allSetKmers(reads, k, false);

			if(verbose){
				std::cout << "\tIndexing Time: " << watch.totalTimeFormatted(0) << std::endl;
				std::cout << "\tCalculating Distances" << std::endl;
			}
			watch.reset();
			distanceMaps[k] = getDistance(reads, numThreads, disFun);
			if(verbose){
				std::cout << "\tCalculating Distances Time: " << watch.totalTimeFormatted(0) << std::endl;
			}
		}
		distances = distanceMaps[kmerStart];
		for (const auto & rowPos : iter::range(distances.size())) {
			for (const auto & colPos : iter::range(distances[rowPos].size())) {
				std::vector<double> differences;
				for (uint32_t k = kmerStart; k < kmerStop; ++k) {
					differences.emplace_back(
							std::abs(
									distanceMaps[k][rowPos][colPos]
											- distanceMaps[k + 1][rowPos][colPos]));
				}
				distances[rowPos][colPos] = vectorMean(differences);
			}
		}
	}
	return distances;
}

std::vector<kmerCluster> greedyKmerSimCluster(const SeqIOOptions & inReadsOpts,
		uint32_t kLength, double kmerSimCutOff, bool checkComplement,
		bool verbose) {

	SeqInput reader(inReadsOpts);
	reader.openIn();
	return greedyKmerSimCluster(reader.readAllReads<readObject>(), kLength,
			kmerSimCutOff, checkComplement, verbose);
}

}  // namespace bibseq
