#pragma once
//
//  collapser.h
//
//  Created by Nicholas Hathaway on 1/1/14.
//

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
#include "njhseq/alignment.h"
#include "njhseq/seqToolsUtils.h"
#include "njhseq/readVectorManipulation.h"
#include "njhseq/objects/kmer/kmerCalculator.hpp"
#include "njhseq/objects/seqObjects/Clusters/cluster.hpp"
#include "njhseq/objects/collapseObjects/opts.h"
#include "njhseq/objects/dataContainers/tables/table.hpp"


namespace njhseq {



struct SnapShotsOpts {
	SnapShotsOpts(bool snapShots, const std::string & snapShotsDirName);
	SnapShotsOpts();
	bool snapShots_ = false;
	std::string snapShotsDirName_ = "snapShots";
};

class collapser {

public:
	// constructor

	collapser(const CollapserOpts & opts);

  CollapserOpts opts_;
private:
	template<class CLUSTER>
	void findMatch(CLUSTER &read, std::vector<CLUSTER> &comparingReads,
			const std::vector<uint64_t> & positions,
			const IterPar &runParams, size_t & amountAdded,
			aligner &alignerObj) const;

	template<class CLUSTER>
	void collapseWithParameters(std::vector<CLUSTER> &comparingReads,
			const IterPar &runParams, aligner &alignerObj) const;

	template<class CLUSTER>
	void collapseWithParameters(std::vector<CLUSTER> &comparingReads,
			std::vector<uint64_t> & positions, const IterPar &runParams,
			aligner &alignerObj) const;

public:

	template<class CLUSTER>
	std::vector<CLUSTER> collapseLowFreqOneOffs(std::vector<CLUSTER> &currentClusters,
			double lowFreqMultiplier,
			aligner &alignerObj,
			bool skipChimeras = true) const;


	template<class CLUSTER>
	std::vector<CLUSTER> runClustering(std::vector<CLUSTER> &currentClusters,
			CollapseIterations iteratorMap,
			aligner &alignerObj) const;

	template<class CLUSTER>
	void runClustering(std::vector<CLUSTER> &currentClusters,
			std::vector<uint64_t> & positons,
			CollapseIterations iteratorMap,
			aligner &alignerObj) const;



  void runFullClustering(std::vector<cluster> & clusters,
			CollapseIterations iteratorMap, CollapseIterations binIteratorMap,
			aligner & alignerObj, const std::string & mainDirectory,
			SeqIOOptions ioOpts, SeqIOOptions refOpts,
			const SnapShotsOpts & snapShotsOpts);



  template<typename READ>
	table markChimeras(std::vector<READ> &processedReads, aligner &alignerObj,
			const ChimeraOpts & chiOpts) const;


	table markChimerasTest(std::vector<cluster> &processedReads, aligner &alignerObj,
			const ChimeraOpts & chiOpts) const;


};

template <class CLUSTER>
void collapser::findMatch(CLUSTER &read,
                                 std::vector<CLUSTER> &comparingReads,
                                 const std::vector<uint64_t> & positions,
                                 const IterPar &runParams,
                                 size_t &amountAdded,
																 aligner &alignerObj) const{

	uint32_t count = 0;
  double bestScore = 0;
  bool foundMatch = false;
  uint32_t bestSearching = 0;
  uint32_t bestClusterPos = std::numeric_limits<uint32_t>::max();
  for (const auto &clusPos : positions) {
    if (comparingReads[clusPos].remove) {
      continue;
    }
  	auto & clus = comparingReads[clusPos];

  	if (clus.seqBase_.cnt_ <= runParams.smallCheckStop_) {
      continue;
    }

    if (clus.seqBase_.name_ == read.seqBase_.name_) {
      continue;
    }

  	if (clus.seqBase_.cnt_ < read.seqBase_.cnt_) {
      continue;
    }

    if (foundMatch) {
      ++bestSearching;
    }
    ++count;
    if (count > runParams.stopCheck_ || bestSearching > opts_.bestMatchOpts_.bestMatchCheck_) {
      break;
    }

    if (opts_.skipOpts_.skipOnLetterCounterDifference_) {
      double sum = read.counter_.getFracDifference(clus.counter_,
      		read.counter_.alphabet_);
      if (sum > opts_.skipOpts_.fractionDifferenceCutOff_) {
        continue;
      }
    }

    if (alignerObj.CountEndGaps()) {
      if (uAbsdiff(len(getSeqBase(read)), len(getSeqBase(clus))) > 20) {
        continue;
      }
    }

    if (opts_.skipOpts_.useReadLen_){
      if (uAbsdiff(len(getSeqBase(read)), len(getSeqBase(clus))) > opts_.skipOpts_.readLenDiff_) {
        continue;
      }
    }

//    if(5 == runParams.iterNumber_ && getSeqBase(clus).cnt_ > 1000 && getSeqBase(read).cnt_ > 1000){
//    	std::cout << std::endl << std::endl;
//    	std::cout << "clus name: "<< getSeqBase(clus).name_ << std::endl;
//    	std::cout << "read name: "<< getSeqBase(read).name_ << std::endl;
//    }
    bool matching = clus.compare(read, alignerObj, runParams, opts_);
//    if(5 == runParams.iterNumber_ && getSeqBase(clus).cnt_ > 1000 && getSeqBase(read).cnt_ > 1000){
//    	std::cout << "\tmatching:  "<< njh::colorBool(matching) << std::endl;
//    }
		if (matching) {
			foundMatch = true;
			if (opts_.bestMatchOpts_.findingBestMatch_) {
				if (opts_.alignOpts_.noAlign_) {
					alignerObj.noAlignSetAndScore(clus, read);
				} else {
					//alignerObj.alignCacheGlobal(clus, read);
					alignerObj.alignCacheGlobalDiag(clus, read);
				}
      	double currentScore = 0;
      	if(opts_.alignOpts_.eventBased_){
      		alignerObj.profilePrimerAlignment(clus, read);
      		currentScore = alignerObj.comp_.distances_.eventBasedIdentity_;
      	}else{
      		currentScore = alignerObj.parts_.score_;
      	}
        if (currentScore > bestScore) {
          bestScore = currentScore;
          bestClusterPos = clusPos;
        }
      } else {
        ++amountAdded;
        clus.addRead(read);
        read.remove = true;
        break;
      }
    }
  }
  if (foundMatch) {
    if (opts_.bestMatchOpts_.findingBestMatch_) {
    	if(bestClusterPos != std::numeric_limits<uint32_t>::max()){
    		comparingReads[bestClusterPos].addRead(read);
        read.remove = true;
        ++amountAdded;
    	}
    }
  }

}


template<class CLUSTER>
void collapser::collapseWithParameters(std::vector<CLUSTER> &comparingReads,
		const IterPar &runParams, aligner &alignerObj) const {
	std::vector<uint64_t> positions(comparingReads.size());
	njh::iota<uint64_t>(positions, 0);
	collapseWithParameters(comparingReads, positions, runParams, alignerObj);
}

template<class CLUSTER>
void collapser::collapseWithParameters(std::vector<CLUSTER> &comparingReads,
		std::vector<uint64_t> & positions, const IterPar &runParams,
		aligner &alignerObj) const {

	uint32_t sizeOfReadVector = 0;
	for (const auto & pos : positions) {
		if (!comparingReads[pos].remove) {
			++sizeOfReadVector;
		}
	}
	if (sizeOfReadVector < 2) {
		return;
	}
	if (opts_.verboseOpts_.verbose_) {
		std::cout << "Starting with " << sizeOfReadVector << " clusters"
				<< std::endl;
	}
	uint32_t clusterCounter = 0;
	size_t amountAdded = 0;
	for (const auto &reverseReadPos : iter::reversed(positions)) {
		auto & reverseRead = comparingReads[reverseReadPos];
		if (reverseRead.remove) {
			continue;
		} else {
			++clusterCounter;
		}
		if (opts_.verboseOpts_.verbose_ && clusterCounter % 100 == 0) {
			std::cout << "\r" << "Currently on cluster " << clusterCounter << " of "
					<< sizeOfReadVector;
			std::cout.flush();
		}
		findMatch(reverseRead, comparingReads, positions, runParams, amountAdded,
				alignerObj);
	}

	njh::stopWatch watch;
	watch.setLapName("updatingName");
	for (const auto & pos : positions) {
		if(!comparingReads[pos].remove){
			comparingReads[pos].updateName();
		}
	}
	watch.startNewLap("allCalculateConsensus");
	for (const auto & pos : positions) {
		if(!comparingReads[pos].remove){
			comparingReads[pos].calculateConsensus(alignerObj, true);
		}
	}
	watch.startNewLap("removeLowQualityBases");
	if (opts_.iTOpts_.removeLowQualityBases_) {
		for (const auto & pos : positions) {
			if (!comparingReads[pos].remove) {
				comparingReads[pos].seqBase_.removeLowQualityBases(
				opts_.iTOpts_.lowQualityBaseTrim_);
			}
		}
	}
	watch.startNewLap("adjustHomopolyerRuns");
	if (opts_.iTOpts_.adjustHomopolyerRuns_) {
		for (const auto & pos : positions) {
			if (!comparingReads[pos].remove) {
				comparingReads[pos].adjustHomopolyerRunQualities();
			}
		}
	}

	if (opts_.verboseOpts_.debug_) {
		watch.logLapTimes(std::cout, true, 6, true);
	}

	if (opts_.verboseOpts_.verbose_) {
		if (clusterCounter >= 100) {
			std::cout << std::endl;
		}
		std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
				<< " clusters" << std::endl;
	}
	if(opts_.clusOpts_.converge_ && amountAdded > 0){
		collapseWithParameters(comparingReads, positions, runParams, alignerObj);
	}

}

template<class CLUSTER>
std::vector<CLUSTER> collapser::runClustering(
		std::vector<CLUSTER> &currentClusters,
		CollapseIterations iteratorMap, aligner &alignerObj) const {
	for (const auto & iter : iteratorMap.iters_) {
		if (opts_.verboseOpts_.verbose_) {
			std::cout << std::endl;
			iter.second.printIterInfo(std::cout, true);
		}
		if (iter.second.errors_.hqMismatches_ == 0
				&& iter.second.errors_.lqMismatches_ == 0
				&& iter.second.errors_.largeBaseIndel_ == 0) {
			std::vector<uint64_t> positions(currentClusters.size());
			njh::iota<uint64_t>(positions, 0);
			auto byCondensed = readVec::organizeByCondensedSeqPositions(
					currentClusters, positions);
			for (auto& condensedReads : byCondensed) {
				collapseWithParameters(currentClusters, condensedReads.second,
						iter.second, alignerObj);
			}
		} else {
			collapseWithParameters(currentClusters, iter.second, alignerObj);
		}
		njh::stopWatch watch;
		watch.setLapName("sorting vector");
		readVecSorter::sortReadVector(currentClusters, "totalCount");
		if (opts_.verboseOpts_.debug_) {
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
	return readVecSplitter::splitVectorOnRemove(currentClusters).first;
}


template<class CLUSTER>
std::vector<CLUSTER> collapser::collapseLowFreqOneOffs(
		std::vector<CLUSTER> &comparingReads, double lowFreqMultiplier,
		aligner &alignerObj,
		bool skipChimeras) const {
	std::vector<uint64_t> positions(comparingReads.size());
	njh::iota<uint64_t>(positions, 0);
	uint32_t sizeOfReadVector = 0;
	for (const auto & pos : positions) {
		if (!comparingReads[pos].remove) {
			++sizeOfReadVector;
		}
	}
	if (sizeOfReadVector < 2) {
		return comparingReads;
	}
	if (opts_.verboseOpts_.verbose_) {
		std::cout << "Starting with " << sizeOfReadVector << " clusters"
				<< std::endl;
	}
	uint32_t clusterCounter = 0;
	size_t amountAdded = 0;
	for (const auto &reverseReadPos : iter::reversed(positions)) {
		auto & reverseRead = comparingReads[reverseReadPos];

		if (reverseRead.remove) {
			continue;
		} else {
			++clusterCounter;
		}
		if(skipChimeras && reverseRead.seqBase_.isChimeric()){
			continue;
		}
		if (opts_.verboseOpts_.verbose_ && clusterCounter % 100 == 0) {
			std::cout << "\r" << "Currently on cluster " << clusterCounter << " of "
					<< sizeOfReadVector;
			std::cout.flush();
		}
		uint32_t count = 0;
	  for (const auto &clusPos : positions) {
	    if (comparingReads[clusPos].remove) {
	      continue;
	    }
	  	auto & clus = comparingReads[clusPos];
	  	if (clus.seqBase_.frac_ <= reverseRead.seqBase_.frac_ * lowFreqMultiplier) {
	      continue;
	    }
	    if (clus.seqBase_.name_ == reverseRead.seqBase_.name_) {
	      continue;
	    }
	    ++count;
	    comparison comp = clus.getComparison(reverseRead, alignerObj, false);
	    //can only get here if clus.seqBase_.frac >  reverseRead.seqBase_.frac_ * lowFreqMultiplier so can just check if only diffs by 1 mismatch
//	    bool matching = ((comp.hqMismatches_ + comp.lqMismatches_ + comp.lowKmerMismatches_) <=1
	    bool matching = ((comp.hqMismatches_) <=1
	    		&& comp.largeBaseIndel_ == 0
					&& comp.twoBaseIndel_ == 0
					&& comp.oneBaseIndel_ == 0);
	    //also add if only different by one 1 base indel
	    if(!matching){
	    	if(comp.distances_.alignmentGaps_.size() == 1
	    			&& comp.distances_.alignmentGaps_.begin()->second.size_ == 1
						&& comp.largeBaseIndel_ == 0
						&& comp.twoBaseIndel_ == 0 &&
						(comp.hqMismatches_ + comp.lqMismatches_ + comp.lowKmerMismatches_) == 0){
	    		matching = true;
	    	}
	    }
			if (matching) {
        ++amountAdded;
        clus.addRead(reverseRead);
        reverseRead.remove = true;
        break;
	    }
	  }
	}

	njh::stopWatch watch;
	watch.setLapName("updatingName");
	for (const auto & pos : positions) {
		if(!comparingReads[pos].remove){
			comparingReads[pos].updateName();
		}
	}
	watch.startNewLap("allCalculateConsensus");
	for (const auto & pos : positions) {
		if(!comparingReads[pos].remove){
			comparingReads[pos].calculateConsensus(alignerObj, true);
		}
	}
	watch.startNewLap("removeLowQualityBases");
	if (opts_.iTOpts_.removeLowQualityBases_) {
		for (const auto & pos : positions) {
			if (!comparingReads[pos].remove) {
				comparingReads[pos].seqBase_.removeLowQualityBases(
				opts_.iTOpts_.lowQualityBaseTrim_);
			}
		}
	}
	watch.startNewLap("adjustHomopolyerRuns");
	if (opts_.iTOpts_.adjustHomopolyerRuns_) {
		for (const auto & pos : positions) {
			if (!comparingReads[pos].remove) {
				comparingReads[pos].adjustHomopolyerRunQualities();
			}
		}
	}

	if (opts_.verboseOpts_.debug_) {
		watch.logLapTimes(std::cout, true, 6, true);
	}

	if (opts_.verboseOpts_.verbose_) {
		if (clusterCounter >= 100) {
			std::cout << std::endl;
		}
		std::cout << "Collapsed down to " << sizeOfReadVector - amountAdded
				<< " clusters" << std::endl;
	}
	if(opts_.clusOpts_.converge_ && amountAdded > 0){
		collapseLowFreqOneOffs(comparingReads, lowFreqMultiplier, alignerObj);
	}
	return readVecSplitter::splitVectorOnRemove(comparingReads).first;
}

template<class CLUSTER>
void collapser::runClustering(std::vector<CLUSTER> &currentClusters,
		std::vector<uint64_t> & positions,
		CollapseIterations iteratorMap,
		aligner &alignerObj) const{
	{
		njh::stopWatch watch;
		auto comp =
				[&currentClusters](const uint64_t & pos1, const uint64_t & pos2) {
					return currentClusters[pos1] < currentClusters[pos2];
				};
		watch.setLapName("sortReadVector");
		njh::sort(positions, comp);
		if (opts_.verboseOpts_.debug_) {
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
	for (const auto & iter : iteratorMap.iters_) {
		if (opts_.verboseOpts_.verbose_) {
			std::cout << std::endl;
			iter.second.printIterInfo(std::cout, true);
		}
		if (iter.second.errors_.hqMismatches_ == 0
				&& iter.second.errors_.lqMismatches_ == 0
				&& iter.second.errors_.largeBaseIndel_ == 0) {
			auto byCondensed = readVec::organizeByCondensedSeqPositions(
					currentClusters, positions);
			for (auto& condensedReads : byCondensed) {
				collapseWithParameters(currentClusters, condensedReads.second,
						iter.second, alignerObj);
			}
		} else {
			collapseWithParameters(currentClusters, positions, iter.second,
					alignerObj);
		}
		njh::stopWatch watch;
		auto comp =
				[&currentClusters](const uint64_t & pos1, const uint64_t & pos2) {
					return currentClusters[pos1] < currentClusters[pos2];
				};
		watch.setLapName("sortReadVector");
		njh::sort(positions, comp);
		if (opts_.verboseOpts_.debug_) {
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
}

template<typename READ>
table collapser::markChimeras(std::vector<READ> &processedReads,
		aligner &alignerObj, const ChimeraOpts & chiOpts) const {
	table ret(VecStr { "read", "readCnt",
		"parent1", "parent1Cnt", "parent1Ratio", "parent1SeqPos",
		"parent2", "parent2Cnt", "parent2Ratio", "parent2SeqPos" });
	//assumes clusters are coming in sorted by total count
	if (processedReads.size() <= 2) {
		return ret;
	}
	if (opts_.verboseOpts_.verbose_) {
		std::cout << "Marking Chimeras" << std::endl;
		std::cout << "Initial Pass" << std::endl;
	}
	//subPos is the position of the clusters to investigate for being possibly chimeric
	for (const auto & subPos : iter::range<size_t>(2, processedReads.size())) {
		if (opts_.verboseOpts_.verbose_) {
			std::cout << subPos << ":" << processedReads.size() << "\r";
			std::cout.flush();
		}
		//skip if the clusters to investigate fall below or is equal to the run cut off
		//normally set to 1 as singlets are going to be thrown out anyways
		if (getSeqBase(processedReads[subPos]).cnt_ <= chiOpts.runCutOff_) {
			continue;
		}
		std::map<size_t, std::vector<size_t>> endChiPos;
		std::map<size_t, std::vector<size_t>> frontChiPos;
		//pos is the position of the possible parents
		for (const auto &pos : iter::range<size_t>(0, subPos)) {
			if (opts_.verboseOpts_.verbose_ && opts_.verboseOpts_.debug_) {
				std::cout << std::endl;
				std::cout << getSeqBase(processedReads[pos]).name_ << std::endl;
				std::cout << getSeqBase(processedReads[subPos]).name_ << std::endl;
				std::cout << "getSeqBase(processedReads[pos]).cnt_"
						<< getSeqBase(processedReads[pos]).cnt_ << std::endl;
				std::cout << "getSeqBase(processedReads[pos]).frac_"
						<< getSeqBase(processedReads[pos]).frac_ << std::endl;
				std::cout << "getSeqBase(processedReads[subPos]).cnt_"
						<< getSeqBase(processedReads[subPos]).cnt_ << std::endl;
				std::cout << "getSeqBase(processedReads[subPos]).frac_"
						<< getSeqBase(processedReads[subPos]).frac_ << std::endl;
			}
			//skip if the proportion of reads is less than parentFreqs
			if ((getSeqBase(processedReads[pos]).cnt_
					/ getSeqBase(processedReads[subPos]).cnt_) < chiOpts.parentFreqs_) {
				continue;
			}
			//skip if the same exact sequence
			if (getSeqBase(processedReads[pos]).seq_
					== getSeqBase(processedReads[subPos]).seq_) {
				continue;
			}
			//global align
			alignerObj.alignCacheGlobal(processedReads[pos], processedReads[subPos]);
			//profile alignment
			alignerObj.profileAlignment(processedReads[pos], processedReads[subPos],
					false, true, false);
			if (opts_.verboseOpts_.verbose_ && opts_.verboseOpts_.debug_) {
				std::cout << "alignerObj.mismatches_.size()"
						<< alignerObj.comp_.distances_.mismatches_.size() << std::endl;
				std::cout << "alignerObj.alignmentGaps_.size()"
						<< alignerObj.comp_.distances_.alignmentGaps_.size() << std::endl;
			}
			std::map<uint32_t, mismatch> savedMismatches =
					alignerObj.comp_.distances_.mismatches_;
			//remove any low quality errors
			if (!chiOpts.keepLowQaulityMismatches_) {
				std::vector<uint32_t> lowQualMismatches;
				for (const auto & mis : savedMismatches) {
					if (!mis.second.highQuality(alignerObj.qScorePars_)) {
						lowQualMismatches.emplace_back(mis.first);
					}
				}
				for (const auto & er : lowQualMismatches) {
					savedMismatches.erase(er);
				}
			}

			if (alignerObj.comp_.distances_.mismatches_.size() > 0
					|| alignerObj.comp_.twoBaseIndel_ > 0
					|| alignerObj.comp_.largeBaseIndel_ > 0) {
				//save the current mismatches and gap infos

				auto savedGapInfo = alignerObj.comp_.distances_.alignmentGaps_;
				size_t frontPos = 0;
				size_t endPos = std::numeric_limits<size_t>::max();
				if (savedMismatches.size() > 0) {
					bool passEnd = true;
					bool passFront = true;
					//check front until first mismatch
					if (savedMismatches.begin()->second.seqBasePos + 1
							<= chiOpts.overLapSizeCutoff_) {
						passFront = false;
					} else {
						alignerObj.profileAlignment(processedReads[pos],
								processedReads[subPos], false, true, false, 0,
								savedMismatches.begin()->first);
						passFront = chiOpts.chiOverlap_.passErrorProfile(alignerObj.comp_);
						if (passFront) {
							frontPos = savedMismatches.begin()->second.seqBasePos;
							//processedReads[subPos].frontChiPos.insert( {
							//	savedMismatches.begin()->second.seqBasePos, pos });
						}
					}
					//check end from last mismatch
					if (chiOpts.overLapSizeCutoff_
							>= (len(getSeqBase(processedReads[subPos]))
									- savedMismatches.rbegin()->second.seqBasePos)) {
						passEnd = false;
					} else {
						alignerObj.profileAlignment(processedReads[pos],
								processedReads[subPos], false, true, false,
								savedMismatches.rbegin()->first + 1);
						passEnd = chiOpts.chiOverlap_.passErrorProfile(alignerObj.comp_);
						if (passEnd) {
							endPos = savedMismatches.rbegin()->second.seqBasePos;

							//processedReads[subPos].endChiPos.insert( {
							//		savedMismatches.rbegin()->second.seqBasePos, pos });
						}
					}
				}

				if (savedGapInfo.size() > 0) {
					//check front until first large gap
					for (const auto & g : savedGapInfo) {
						if (g.second.size_ <= 1) {
							continue;
						}
						if (g.second.seqPos_ + 1 <= chiOpts.overLapSizeCutoff_) {
							continue;
						}
						alignerObj.profileAlignment(processedReads[pos],
								processedReads[subPos], false, true, false, 0,
								g.second.startPos_);
						if (chiOpts.chiOverlap_.passErrorProfile(alignerObj.comp_)) {
							size_t gapPos = getRealPosForAlnPos(
									alignerObj.alignObjectB_.seqBase_.seq_, g.second.startPos_);
							if (gapPos > frontPos) {
								frontPos = gapPos;
							}
							//processedReads[subPos].frontChiPos.insert(
							//		{ getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_, g.second.startPos_), pos });
						}
						break;
					}

					//check end from last large gap
					for (const auto & g : iter::reversed(savedGapInfo)) {
						if (g.second.size_ <= 1) {
							continue;
						}
						if (len(getSeqBase(processedReads[subPos]))
								- (g.second.seqPos_ + g.second.size_)
								<= chiOpts.overLapSizeCutoff_) {
							continue;
						}
						alignerObj.profileAlignment(processedReads[pos],
								processedReads[subPos], false, true, false,
								g.second.startPos_ + g.second.size_);

						if (chiOpts.chiOverlap_.passErrorProfile(alignerObj.comp_)) {
							size_t gapPos = 0;
							if (g.second.ref_) {
								gapPos = getRealPosForAlnPos(
										alignerObj.alignObjectB_.seqBase_.seq_, g.second.startPos_)
										+ g.second.size_ - 1;
								//processedReads[subPos].endChiPos.insert(
								//		{ getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
								//				g.second.startPos_) + g.second.size_ - 1, pos });
							} else {
								/**@todo something else should probably happen if this is 0 */
								if (0
										== getRealPosForAlnPos(
												alignerObj.alignObjectB_.seqBase_.seq_,
												g.second.startPos_)) {
									gapPos = getRealPosForAlnPos(
											alignerObj.alignObjectB_.seqBase_.seq_,
											g.second.startPos_);
									//processedReads[subPos].endChiPos.insert(
									//		{ getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
									//				g.second.startPos_), pos });
								} else {
									gapPos = getRealPosForAlnPos(
											alignerObj.alignObjectB_.seqBase_.seq_,
											g.second.startPos_) - 1;
									//processedReads[subPos].endChiPos.insert(
									//		{ getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
									//				g.second.startPos_) - 1, pos });
								}
							}
							if (gapPos < endPos) {
								endPos = gapPos;
							}
						}
						break;
					}
				} // pass gaps size > 0
				//check to see if we have an endPos
				if (std::numeric_limits<size_t>::max() != endPos) {
					endChiPos[endPos].emplace_back(pos);
				}
				//check to see if we have a frontPos
				if (0 != frontPos) {
					frontChiPos[frontPos].emplace_back(pos);
				}
			} // pass has errors
		} // pos loop

		struct ChiParLocs {
			ChiParLocs(size_t endReadMinReadPos, size_t endPosition,
					size_t frontReadMinReadPos, size_t frontPosition) :
					endReadMinReadPos_(endReadMinReadPos),
					endPosition_(endPosition),
					frontReadMinReadPos_(frontReadMinReadPos),
					frontPosition_(frontPosition) {

			}
			size_t endReadMinReadPos_;
			size_t endPosition_;

			size_t frontReadMinReadPos_;
			size_t frontPosition_;

		};

		if (!endChiPos.empty() && !frontChiPos.empty()) {
			if (endChiPos.begin()->first < frontChiPos.rbegin()->first) {
				std::vector<ChiParLocs> crossOvers;
				for(const auto & endChi : endChiPos){
					for(const auto & frontChi : iter::reversed(frontChiPos)){
						if(endChi.first + chiOpts.posSpacing_ < frontChi.first){
							for(const auto & endReadPos : endChi.second){
								for(const auto & frontReadPos : frontChi.second){
									if(endReadPos != frontReadPos){
										crossOvers.emplace_back(ChiParLocs{
											endReadPos, endChi.first,
											frontReadPos, frontChi.first});
									}
								}
							}
						}
					}
					if(!crossOvers.empty()){
						break;
					}
				}
				if(!crossOvers.empty()){
					getSeqBase(processedReads[subPos]).markAsChimeric();
					njh::sort(crossOvers, [](const ChiParLocs & p1, const ChiParLocs & p2){
						if(p1.endReadMinReadPos_ == p2.endReadMinReadPos_){
							return p1.frontReadMinReadPos_ < p2.frontReadMinReadPos_;
						}else{
							return p1.endReadMinReadPos_ < p2.endReadMinReadPos_;
						}
					});

					ret.content_.emplace_back(
							toVecStr(
									getSeqBase(processedReads[subPos]).name_,
									getSeqBase(processedReads[subPos]).cnt_,
									getSeqBase(processedReads[crossOvers.front().frontReadMinReadPos_]).name_,
									getSeqBase(processedReads[crossOvers.front().frontReadMinReadPos_]).cnt_,
									getSeqBase(processedReads[crossOvers.front().frontReadMinReadPos_]).cnt_/getSeqBase(processedReads[subPos]).cnt_,
									crossOvers.front().frontPosition_,
									getSeqBase(processedReads[crossOvers.front().endReadMinReadPos_]).name_,
									getSeqBase(processedReads[crossOvers.front().endReadMinReadPos_]).cnt_,
									getSeqBase(processedReads[crossOvers.front().endReadMinReadPos_]).cnt_/getSeqBase(processedReads[subPos]).cnt_,
									crossOvers.front().endPosition_));
				}
			} // is chimeric loop
		}
	} // subPos loop
	if (opts_.verboseOpts_.verbose_) {
		std::cout << std::endl;
	}
	return ret;
}

}  // namespace njhseq

