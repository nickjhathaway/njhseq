//
//  collapser.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 1/1/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

#include "collapser.hpp"
#include "bibseq/objects/helperObjects/nucCompCluster.hpp"
#include "bibseq/helpers/profiler.hpp"
#include "bibseq/objects/seqObjects/Clusters/clusterUtils.hpp"

namespace bibseq {




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
  }

  readVecSorter::sortReadVector(clusters, "totalCount");
  //create snapshots directory if taking snap shots
  std::string snapShotsDirectoryName = "";
	if (snapShotsOpts.snapShots_) {
		snapShotsDirectoryName = bib::files::makeDir(mainDirectory,
				bib::files::MkdirPar(snapShotsOpts.snapShotsDirName_, false));
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
			bib::iota<uint64_t>(positions, 0);
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
    	collapseWithParameters(clusters, iter.second, alignerObj);
    }
		bib::stopWatch watch;
		watch.setLapName("sorting vector");
		readVecSorter::sortReadVector(clusters, "totalCount");
		if (opts_.verboseOpts_.debug_) {
			watch.logLapTimes(std::cout, true, 6, true);
		}
    //write out current iteration if taking snapshots
    if (snapShotsOpts.snapShots_) {
      std::string iterDir =
          bib::files::makeDir(snapShotsDirectoryName, bib::files::MkdirPar(std::to_string(iter.first), false));
      std::vector<cluster> currentClusters =
          readVecSplitter::splitVectorOnRemove(clusters).first;
      std::string seqName = bib::files::getFileName(ioOpts.firstName_);
      renameReadNames(currentClusters, seqName, true, false);
      SeqOutput writer(SeqIOOptions(snapShotsDirectoryName + std::to_string(iter.first),
      		ioOpts.outFormat_,ioOpts.out_));
      writer.openWrite(currentClusters);
      clusterVec::allWriteClustersInDir(currentClusters, iterDir,ioOpts );
      if (refOpts.firstName_ == "") {
        profiler::getFractionInfoCluster(currentClusters, snapShotsDirectoryName,
                                  std::to_string(iter.first) + ".tab.txt");
      } else {
        profiler::getFractionInfoCluster(currentClusters, snapShotsDirectoryName,
                                  std::to_string(iter.first) + ".tab.txt",
																	refOpts.firstName_, alignerObj, false,
                                  true);
      }
    }
  }
  clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
}




}  // namespace bibseq



