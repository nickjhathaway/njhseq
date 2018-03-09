#pragma once
//
//  collapser.h
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/1/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "bibseq/utils.h"
#include "bibseq/alignment.h"
#include "bibseq/seqToolsUtils.h"
#include "bibseq/readVectorManipulation.h"
#include "bibseq/objects/kmer/kmerCalculator.hpp"
#include "bibseq/objects/seqObjects/Clusters/cluster.hpp"
#include "bibseq/objects/collapseObjects/opts.h"
#include "bibseq/objects/dataContainers/tables/table.hpp"


namespace bibseq {



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


    bool matching = clus.compare(read, alignerObj, runParams, opts_);
    if (matching) {
      foundMatch = true;
      if (opts_.bestMatchOpts_.findingBestMatch_) {
        if(opts_.alignOpts_.noAlign_){
        	alignerObj.noAlignSetAndScore(clus, read);
        }else{
        	alignerObj.alignCacheGlobal(clus, read);
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
	bib::iota<uint64_t>(positions, 0);
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

	bib::stopWatch watch;
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
			bib::iota<uint64_t>(positions, 0);
			auto byCondensed = readVec::organizeByCondensedSeqPositions(
					currentClusters, positions);
			for (auto& condensedReads : byCondensed) {
				collapseWithParameters(currentClusters, condensedReads.second,
						iter.second, alignerObj);
			}
		} else {
			collapseWithParameters(currentClusters, iter.second, alignerObj);
		}
		bib::stopWatch watch;
		watch.setLapName("sorting vector");
		readVecSorter::sortReadVector(currentClusters, "totalCount");
		if (opts_.verboseOpts_.debug_) {
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
	return readVecSplitter::splitVectorOnRemove(currentClusters).first;
}

template<class CLUSTER>
void collapser::runClustering(std::vector<CLUSTER> &currentClusters,
		std::vector<uint64_t> & positions,
		CollapseIterations iteratorMap,
		aligner &alignerObj) const{
	{
		bib::stopWatch watch;
		auto comp =
				[&currentClusters](const uint64_t & pos1, const uint64_t & pos2) {
					return currentClusters[pos1] < currentClusters[pos2];
				};
		watch.setLapName("sortReadVector");
		bib::sort(positions, comp);
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
		bib::stopWatch watch;
		auto comp =
				[&currentClusters](const uint64_t & pos1, const uint64_t & pos2) {
					return currentClusters[pos1] < currentClusters[pos2];
				};
		watch.setLapName("sortReadVector");
		bib::sort(positions, comp);
		if (opts_.verboseOpts_.debug_) {
			watch.logLapTimes(std::cout, true, 6, true);
		}
	}
}

template<typename READ>
table collapser::markChimeras(std::vector<READ> &processedReads,
		aligner &alignerObj, const ChimeraOpts & chiOpts) const {
	table ret(VecStr { "read", "parent1", "parent1SeqPos", "parent2",
			"parent2SeqPos" });
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
		std::multimap<size_t, size_t> endChiPos;
		std::multimap<size_t, size_t> frontChiPos;
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
					endChiPos.insert( { endPos, pos });
				}
				//check to see if we have a frontPos
				if (0 != frontPos) {
					frontChiPos.insert( { frontPos, pos });
				}
			} // pass has errors
		} // pos loop

		if (!endChiPos.empty() && !frontChiPos.empty()) {
			if (endChiPos.begin()->first < frontChiPos.rbegin()->first) {
				getSeqBase(processedReads[subPos]).markAsChimeric();
				size_t endReadMinPos = endChiPos.begin()->second;
				for (const auto & end : endChiPos) {
					if (end.first != endChiPos.begin()->first) {
						break;
					} else {
						if (end.second < endReadMinPos) {
							endReadMinPos = end.second;
						}
					}
				}
				size_t frontReadMinPos = frontChiPos.rbegin()->second;
				for (const auto & front : iter::reversed(frontChiPos)) {
					if (front.first != frontChiPos.rbegin()->first) {
						break;
					} else {
						if (front.second < frontReadMinPos) {
							frontReadMinPos = front.second;
						}
					}
				}
				ret.content_.emplace_back(
						toVecStr(getSeqBase(processedReads[subPos]).name_,
								getSeqBase(processedReads[frontReadMinPos]).name_,
								frontChiPos.rbegin()->first,
								getSeqBase(processedReads[endReadMinPos]).name_,
								endChiPos.begin()->first));
			} // is chimeric loop
		}
	} // subPos loop
	if (opts_.verboseOpts_.verbose_) {
		std::cout << std::endl;
	}
	return ret;
}

}  // namespace bib

