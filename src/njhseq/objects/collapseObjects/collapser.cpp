//
//  collapser.cpp
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
#include "collapser.hpp"
#include "njhseq/objects/seqObjects/seqKmers/KmerVecUtils.hpp"
#include "njhseq/objects/helperObjects/nucCompCluster.hpp"
#include "njhseq/helpers/profiler.hpp"
#include "njhseq/objects/seqObjects/Clusters/clusterUtils.hpp"

namespace njhseq {


table collapser::markChimerasTest(std::vector<cluster> &processedReads,
                                 aligner &alignerObj,
																 const ChimeraOpts & chiOpts) const{

	table ret(VecStr{"read", "parent1", "parent1SeqPos", "parent2", "parent2SeqPos"});
	//assumes clusters are coming in sorted by total count
	if(processedReads.size() <=2){
		return ret;
	}
	if(opts_.verboseOpts_.verbose_){
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
		if(getSeqBase(processedReads[subPos]).cnt_ <= chiOpts.runCutOff_){
			continue;
		}
		std::multimap<size_t, size_t> endChiPos;
		std::multimap<size_t, size_t> frontChiPos;
		//pos is the position of the possible parents
		for (const auto &pos : iter::range<size_t>(0, subPos)) {
			if(opts_.verboseOpts_.verbose_ && opts_.verboseOpts_.debug_){
				std::cout << std::endl;
				std::cout << getSeqBase(processedReads[pos]).name_ << std::endl;
				std::cout << getSeqBase(processedReads[subPos]).name_ << std::endl;
				std::cout << "getSeqBase(processedReads[pos]).cnt_" << getSeqBase(processedReads[pos]).cnt_ << std::endl;
				std::cout << "getSeqBase(processedReads[pos]).frac_" << getSeqBase(processedReads[pos]).frac_ << std::endl;
				std::cout << "getSeqBase(processedReads[subPos]).cnt_" << getSeqBase(processedReads[subPos]).cnt_ << std::endl;
				std::cout << "getSeqBase(processedReads[subPos]).frac_" << getSeqBase(processedReads[subPos]).frac_ << std::endl;
			}
			//skip if the proportion of reads is less than parentFreqs
			if ((getSeqBase(processedReads[pos]).cnt_
					/ getSeqBase(processedReads[subPos]).cnt_) < chiOpts.parentFreqs_) {
				continue;
			}
			//skip if the same exact sequence
			if(getSeqBase(processedReads[pos]).seq_ == getSeqBase(processedReads[subPos]).seq_){
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
    	std::map<uint32_t, mismatch> savedMismatches = alignerObj.comp_.distances_.mismatches_;
    	//remove any low quality errors
    	std::vector<uint32_t> lowQualMismatches;
    	for(const auto & mis : savedMismatches){
    		if(!mis.second.highQuality(alignerObj.qScorePars_)){
    			lowQualMismatches.emplace_back(mis.first);
    		}
    	}
    	for(const auto & er : lowQualMismatches){
    		savedMismatches.erase(er);
    	}
      if (alignerObj.comp_.distances_.mismatches_.size() > 0
      		|| alignerObj.comp_.twoBaseIndel_ > 0
					|| alignerObj.comp_.largeBaseIndel_ > 0) {
      	//save the current mismatches and gap infos

      	auto savedGapInfo = alignerObj.comp_.distances_.alignmentGaps_;
      	size_t frontPos = 0;
      	size_t endPos = std::numeric_limits<size_t>::max();
      	if(savedMismatches.size() > 0){
        	bool passEnd = true;
        	bool passFront = true;
					//check front until first mismatch
					if (savedMismatches.begin()->second.seqBasePos + 1 <= chiOpts.overLapSizeCutoff_) {
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

      	if(savedGapInfo.size() > 0){
          //check front until first large gap
          for(const auto & g : savedGapInfo){
          	if(g.second.size_ <=1){
          		continue;
          	}
          	if(g.second.seqPos_ + 1 <= chiOpts.overLapSizeCutoff_){
          		continue;
          	}
  					alignerObj.profileAlignment(processedReads[pos],
  							processedReads[subPos], false, true, false,
								0, g.second.startPos_);
  					if (chiOpts.chiOverlap_.passErrorProfile(alignerObj.comp_)) {
  						size_t gapPos = getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_, g.second.startPos_);
  						if(gapPos > frontPos){
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
  					if(len(getSeqBase(processedReads[subPos])) - (g.second.seqPos_ + g.second.size_) <= chiOpts.overLapSizeCutoff_){
  						continue;
  					}
  					alignerObj.profileAlignment(processedReads[pos],
  							processedReads[subPos], false, true, false,
								g.second.startPos_ + g.second.size_);

  					if (chiOpts.chiOverlap_.passErrorProfile(alignerObj.comp_)) {
  						size_t gapPos = 0;
  						if (g.second.ref_) {
  							gapPos = getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,g.second.startPos_) + g.second.size_ - 1;
  							//processedReads[subPos].endChiPos.insert(
  							//		{ getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
  							//				g.second.startPos_) + g.second.size_ - 1, pos });
  						} else {
  							/**@todo something else should probably happen if this is 0 */
  							if (0
  									== getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
  											g.second.startPos_)) {
  								gapPos = getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_, g.second.startPos_);
  								//processedReads[subPos].endChiPos.insert(
  								//		{ getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
  								//				g.second.startPos_), pos });
  							} else {
  								gapPos = getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_, g.second.startPos_) - 1;
  								//processedReads[subPos].endChiPos.insert(
  								//		{ getRealPosForAlnPos(alignerObj.alignObjectB_.seqBase_.seq_,
  								//				g.second.startPos_) - 1, pos });
  							}
  						}
  						if(gapPos < endPos){
  							endPos = gapPos;
  						}
  					}
  					break;
  				}
      	} // pass gaps size > 0
      	//check to see if we have an endPos
      	if(std::numeric_limits<size_t>::max() != endPos){
      		endChiPos.insert({endPos, pos});
      	}
      	//check to see if we have a frontPos
				if (0 != frontPos) {
					frontChiPos.insert({frontPos, pos});
				}
      } // pass has errors
		} // pos loop


		if(!endChiPos.empty() && ! frontChiPos.empty()){
			if(endChiPos.begin()->first < frontChiPos.rbegin()-> first){
				getSeqBase(processedReads[subPos]).markAsChimeric();
				size_t endReadMinPos = endChiPos.begin()->second;
				for(const auto & end : endChiPos){
					if(end.first != endChiPos.begin()->first){
						break;
					}else{
						if(end.second < endReadMinPos){
							endReadMinPos = end.second;
						}
					}
				}
				size_t frontReadMinPos = frontChiPos.rbegin()->second;
				for(const auto & front : iter::reversed(frontChiPos)){
					if(front.first != frontChiPos.rbegin()->first){
						break;
					}else{
						if(front.second < frontReadMinPos){
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
  if(opts_.verboseOpts_.verbose_){
  	std::cout << std::endl;
  }
  return ret;
}




SnapShotsOpts::SnapShotsOpts(){

}


SnapShotsOpts::SnapShotsOpts(bool snapShots, const std::string & snapShotsDirName) :
		snapShots_(snapShots), snapShotsDirName_(snapShotsDirName) {
}


collapser::collapser(const CollapserOpts & opts) :
		opts_(opts) {
}

void collapser::runFullClustering(std::vector<cluster> & clusters,
		CollapseIterations iteratorMap, CollapseIterations binIteratorMap,
		aligner & alignerObj, const std::string & mainDirectory,
		SeqIOOptions ioOpts, SeqIOOptions refOpts,
		const SnapShotsOpts & snapShotsOpts) {
  clusterVec::allSetFractionClusters(clusters);

  if (opts_.nucCompBinOpts_.useNucComp_) {
		if (opts_.verboseOpts_.verbose_) {
			std::cout << "Binning on Nucleotide composition first"
					<< std::endl;
		}
  	uint32_t minLen = readVec::getMinLength(clusters);
  	//get the alphabet of the current cluster set and use that in the nuc comp setting
  	charCounter allCounter;
  	for(const auto & read : clusters){
  		allCounter.increaseCountByString(read.seqBase_.seq_, read.seqBase_.cnt_);
  	}
  	allCounter.resetAlphabet(false);
  	auto alphabet = allCounter.alphabet_;
  	//set the verbosity to debug level
  	bool oldVerbose = opts_.verboseOpts_.verbose_;
  	opts_.verboseOpts_.verbose_ = opts_.verboseOpts_.debug_;
  	//now cluster on nucleotide composition clusters
		if (oldVerbose) {
			std::cout << "Min Len for binning is  " << minLen << std::endl;
		}
  	for(const auto & diffCutOff : opts_.nucCompBinOpts_.diffCutOffVec_){
			std::vector<nucCompCluster> comps;
			if (opts_.nucCompBinOpts_.useMinLenNucComp_) {
				comps = clusterOnNucComp(clusters, minLen, alphabet, diffCutOff,
						opts_.nucCompBinOpts_.findBestNuc_, true,
						opts_.verboseOpts_.debug_);
			} else {
				comps = clusterOnNucComp(clusters, alphabet, diffCutOff,
						opts_.nucCompBinOpts_.findBestNuc_, true,
						opts_.verboseOpts_.debug_);
			}
			if (oldVerbose) {
				std::cout << "On Nucleotide Composition Bin difference of "
						<< diffCutOff << " in " << vectorToString(opts_.nucCompBinOpts_.diffCutOffVec_, ", ")
						<< std::endl;
			}
			//now run clustering on the nucleotide clusters
    	uint32_t compCount = 0;
    	for(auto & comp : comps){
    		++compCount;
    		if(oldVerbose){
      		std::cout << "On " << compCount << " of " << comps.size() << "\r";
      		std::cout.flush();
    		}
    		if(comp.readPositions_.size() < 2){
    			continue;
    		}
    		double currentReadCnt = 0;
    		for(const auto & pos : comp.readPositions_){
    			currentReadCnt += clusters[pos].seqBase_.cnt_;
    		}
    		//set the run cut off of low kmer frequency to the current count
    		auto oldRunCutoff = alignerObj.kMaps_.runCutOff_ ;
    		uint32_t currentRunCutoff = processRunCutoff(opts_.kmerOpts_.runCutOffString_,currentReadCnt);
    		alignerObj.kMaps_.runCutOff_ = currentRunCutoff;
    		//binIteratorMap.writePars(std::cout);
    		runClustering(clusters, comp.readPositions_, binIteratorMap, alignerObj);
    		alignerObj.kMaps_.runCutOff_ = oldRunCutoff;
    	}
    	if(oldVerbose){
    		std::cout << std::endl;
    	}
  	}
  	opts_.verboseOpts_.verbose_ = oldVerbose;
  	//remove all the clusters marked removed
  	clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
  } else if (opts_.kmerBinOpts_.useKmerBinning_){
		if (opts_.verboseOpts_.verbose_) {
			std::cout << "Binning on Kmer composition first"
					<< std::endl;
		}
  	//set verbosity at debug level
  	bool oldVerbose = opts_.verboseOpts_.verbose_;
  	opts_.verboseOpts_.verbose_ = opts_.verboseOpts_.debug_;
  	//create vector of seqWtihKmerInfo reads to use to clusters by kmer distance
  	std::vector<std::unique_ptr<seqWithKmerInfo>> reads;
  	reads.reserve(clusters.size());
  	for(const auto & read : clusters){
  		reads.emplace_back(std::make_unique<seqWithKmerInfo>(read.seqBase_));
  	}
  	allSetKmers(reads, opts_.kmerBinOpts_.kCompareLen_, false);
  	std::vector<kmerClusterPos> kClusters;
  	//now cluster reads based on kmers
  	for(const auto & readPos : iter::range(reads.size())){
  		if(opts_.verboseOpts_.debug_ && readPos % 50 == 0){
  			std::cout << "Currently on " << readPos << " of "
  					<< reads.size() << std::endl;
  		}
  		bool foundMatch = false;
  		for(auto & kClus : kClusters){
  			foundMatch = kClus.compareRead(reads[readPos], readPos, opts_.kmerBinOpts_.kmerCutOff_,
  					false);
  			if(foundMatch){
  				break;
  			}
  		}
  		if(!foundMatch){
  			kClusters.emplace_back(kmerClusterPos(reads[readPos], readPos));
  		}
  	}
  	//run clustering on the different kmer clustering
  	uint32_t clusterCount = 0;
  	for(auto & kClus : kClusters){
  		++clusterCount;
  		if(clusterCount % 50 == 0){
  			std::cout << "On " << clusterCount << " of " << kClusters.size() << std::endl;
  		}
  		if(kClus.readPositions_.size() < 2){
  			continue;
  		}
  		double currentReadCnt = 0;
  		for(const auto & pos : kClus.readPositions_){
  			currentReadCnt += clusters[pos].seqBase_.cnt_;
  		}
  		auto oldRunCutoff = alignerObj.kMaps_.runCutOff_ ;
  		uint32_t currentRunCutoff = processRunCutoff(opts_.kmerOpts_.runCutOffString_, currentReadCnt );
  		alignerObj.kMaps_.runCutOff_ = currentRunCutoff;
  		runClustering(clusters, kClus.readPositions_, binIteratorMap, alignerObj);
  		alignerObj.kMaps_.runCutOff_ = oldRunCutoff;
  	}
  	opts_.verboseOpts_.verbose_ = oldVerbose;
  	//remove the clusters marked removed
  	clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
  }

  readVecSorter::sortReadVector(clusters, "totalCount");
  //create snapshots directory if taking snap shots
  std::string snapShotsDirectoryName = "";
	if (snapShotsOpts.snapShots_) {
		snapShotsDirectoryName = njh::files::makeDir(mainDirectory,
				njh::files::MkdirPar(snapShotsOpts.snapShotsDirName_, false)).string();
	}
  for (const auto & iter : iteratorMap.iters_) {
  	if(opts_.verboseOpts_.verbose_){
  		std::cout << std::endl;
  	}
    if(opts_.verboseOpts_.verbose_){
    	iter.second.printIterInfo(std::cout, true);
    }
    // if collapse is only on indels, create cluerters based on condensed seq and run clustering on those
    if (!iteratorMap.onPerId_ && iter.second.errors_.hqMismatches_ == 0 &&
    		iter.second.errors_.lqMismatches_ == 0 &&
				iter.second.errors_.largeBaseIndel_ == 0) {
			bool oldVerbose = opts_.verboseOpts_.verbose_;
			uint32_t startingSizeOfReadVector = 0;
			for (const auto & read : clusters) {
				if (!read.remove) {
					++startingSizeOfReadVector;
				}
			}
			if (opts_.verboseOpts_.verbose_) {
				std::cout << "Starting with " << startingSizeOfReadVector << " clusters"
						<< std::endl;
			}
			opts_.verboseOpts_.verbose_ = opts_.verboseOpts_.debug_;

			std::vector<uint64_t> positions(clusters.size(), 0);
			njh::iota<uint64_t>(positions, 0);
			auto byCondensed = readVec::organizeByCondensedSeqPositions(
					clusters, positions);
			for (auto& condensedReads : byCondensed) {
				collapseWithParameters(clusters, condensedReads.second,
						iter.second, alignerObj);
			}
			opts_.verboseOpts_.verbose_ = oldVerbose;
      int stopSizeOfReadVector = readVec::getReadVectorSize(clusters);
      if(opts_.verboseOpts_.verbose_){
      	std::cout << "Collapsed down to " << stopSizeOfReadVector << std::endl;
      }
    } else {
//    	if(5 == iter.second.iterNumber_){
//    		std::cout << clusters[0].seqBase_.name_ << std::endl;
//    		std::cout << clusters[1].seqBase_.name_ << std::endl;
//    	}
    	collapseWithParameters(clusters, iter.second, alignerObj);
//    	if(5 == iter.second.iterNumber_){
//    		std::cout << clusters[0].seqBase_.name_ << std::endl;
//    		std::cout << clusters[1].seqBase_.name_ << std::endl;
//    	}
    }
		njh::stopWatch watch;
		watch.setLapName("sorting vector");
		readVecSorter::sortReadVector(clusters, "totalCount");
		if (opts_.verboseOpts_.debug_) {
			watch.logLapTimes(std::cout, true, 6, true);
		}
    //write out current iteration if taking snapshots
    if (snapShotsOpts.snapShots_) {

    	std::string iterName = njh::leftPadNumStr<uint32_t>(iter.first, iteratorMap.iters_.size());

      bfs::path iterDir =
          njh::files::makeDir(snapShotsDirectoryName, njh::files::MkdirPar(iterName, false));
      std::vector<cluster> currentClusters =
          readVecSplitter::splitVectorOnRemove(clusters).first;
      std::string seqName = bfs::basename(ioOpts.firstName_);
      renameReadNames(currentClusters, seqName, true, false);
      SeqOutput writer(SeqIOOptions(snapShotsDirectoryName + iterName,
      		ioOpts.outFormat_,ioOpts.out_));
      writer.openWrite(currentClusters);
      clusterVec::allWriteClustersInDir(currentClusters, iterDir.string(),ioOpts );
      if (refOpts.firstName_ == "") {
        profiler::getFractionInfoCluster(currentClusters, snapShotsDirectoryName,
                                  iterName + ".tab.txt");
      } else {
        profiler::getFractionInfoCluster(
        		currentClusters,
						snapShotsDirectoryName,
						iterName + ".tab.txt",
						refOpts.firstName_.string(),
						alignerObj, false);
      }
    }
  }
  clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
}




}  // namespace njhseq



