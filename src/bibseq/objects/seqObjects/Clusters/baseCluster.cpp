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

#include "baseCluster.hpp"
#include "bibseq/helpers/consensusHelper.hpp"
#include "bibseq/readVectorManipulation/readVectorHelpers.h"
#include "bibseq/IO/SeqIO/SeqOutput.hpp"
#include "bibseq/objects/dataContainers/graphs/ConBasePathGraph.hpp"


namespace bibseq {


baseCluster::baseCluster() : readObject() {
  firstReadCount_ = 0;
  needToCalculateConsensus_ = true;
}

baseCluster::baseCluster(const seqInfo& firstRead) : readObject(firstRead) {
  
	firstReadName_ = firstRead.name_;
  firstReadCount_ = firstRead.cnt_;
  reads_.emplace_back(std::make_shared<readObject>(firstRead));
  needToCalculateConsensus_ = true;
  remove = false;
  updateName();
}


void baseCluster::addRead(const baseCluster& otherCluster) {
  seqBase_.cnt_ += otherCluster.seqBase_.cnt_;
  seqBase_.frac_ = ((seqBase_.frac_ * reads_.size()) + otherCluster.seqBase_.frac_) /
                   (reads_.size() + 1);
  if (seqBase_.cnt_ / 2.0 >= firstReadCount_) {
    needToCalculateConsensus_ = true;
  }
  // needToCalculateConsensus = true;
  reads_.insert(reads_.end(), otherCluster.reads_.begin(), otherCluster.reads_.end());
}
void baseCluster::calculateConsensusToCurrent(aligner& alignerObj, bool setToConsensus){
	// if the cluster is only one read, no need to create consensus
	if (reads_.size() <= 1) {
		return;
	}
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus_) {
    return;
  }
  calculateConsensusTo(seqBase_, alignerObj, setToConsensus);

}











void baseCluster::calculateConsensusTo(const seqInfo & seqBase,
		aligner& alignerObj, bool setToConsensus) {

	// create the map for letter counters for each position
	std::map<uint32_t, charCounter> counters;
	// create a map in case of insertions
	std::map<uint32_t, std::map<uint32_t, charCounter>> insertions;
	std::map<int32_t, charCounter> beginningGap;
	auto getSeqBase =
			[](const std::shared_ptr<readObject> & read) ->const seqInfo& {return read->seqBase_;};

	consensusHelper::increaseCounters(seqBase_, reads_, getSeqBase, alignerObj, counters,
			insertions, beginningGap);

	calcConsensusInfo_ = seqBase_;
	for (auto & counter : counters) {
		counter.second.resetAlphabet(true);
		counter.second.setFractions();
	}
	for (auto & counter : beginningGap) {
		counter.second.resetAlphabet(true);
		counter.second.setFractions();
	}
	for (auto & counter : insertions) {
		for (auto & subCounter : counter.second) {
			subCounter.second.resetAlphabet(true);
			subCounter.second.setFractions();
		}
	}

	//if the count is just 2 then just a majority rules consensus
	if(seqBase_.cnt_ > 2){

		//find out if there are several locations with larger contention of majority rules
		//consensus and therefore could lead to bad consensus building
		double contentionCutOff = 0.25;
		uint32_t countAbovepCutOff = 0;
		std::vector<uint32_t> importantPositions;
		for(const auto & counter : counters){
			uint32_t count = 0;
			for(const auto base : counter.second.alphabet_){
				if(counter.second.fractions_[base] > contentionCutOff){
					++count;
					if(count >=2){
						break;
					}
				}
			}
			if(count >=2){
				++countAbovepCutOff;
				importantPositions.emplace_back(counter.first);
			}
		}
		/* for debugging
		if(seqBase_.name_ == "lib03_Pf3D7-AMA1-var0_seq.0001_2" ||
				firstReadName_ == "lib03_Pf3D7-AMA1-var0_seq.0001_2"){
			std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
			std::cout << seqBase_.name_ << std::endl;
			std::cout << "countAbovepCutOff: " << countAbovepCutOff << std::endl;
			bfs::path counterFnp = firstReadName_ + "_baseCounts.tab.txt";
			counterFnp = bib::files::findNonexitantFile(counterFnp.string());
			std::ofstream outFile(counterFnp.string());

			outFile << "pos\tbase\tcount\tfrac" << "\n";
			for (const auto & counter : counters) {
				for (const auto base : counter.second.alphabet_) {
					outFile << counter.first << "\t" << base << "\t"
							<< counter.second.chars_[base] << "\t"
							<< counter.second.fractions_[base] << "\n";
				}
			}
			bfs::path insert_counterFnp = firstReadName_ + "_insertBaseCounts.tab.txt";
			insert_counterFnp = bib::files::findNonexitantFile(insert_counterFnp.string());
			std::ofstream outInsertFile(insert_counterFnp.string());
			outFile << "pos\tinsertPos\tbase\tcount\tfrac" << "\n";
			for (const auto & insert : insertions) {
				for (const auto & counter : insert.second) {
					for (const auto base : counter.second.alphabet_) {
						outInsertFile << insert.first << "\t"
								<< counter.first << "\t" << base << "\t"
								<< counter.second.chars_[base] << "\t"
								<< counter.second.fractions_[base] << "\n";
					}
				}
			}
		}*/
		//if there are several points of contention
		if(countAbovepCutOff >= 2){

			//for debugging;
			/*
			std::ofstream outFile(firstReadName_ + "_baseCounts.tab.txt");
			outFile << "pos\tbase\tcount\tfrac" << "\n";
			for(const auto & counter : counters){
				for(const auto base : counter.second.alphabet_){
					outFile << counter.first
							<< "\t" << base
							<< "\t" << counter.second.chars_[base]
							<< "\t" << counter.second.fractions_[base] <<"\n";
				}
			}*/

			ConBasePathGraph graph;
			for(const auto pos : importantPositions){
				auto & counter = counters.at(pos);
				for(const auto base : counter.alphabet_){
					if(counter.chars_[base] > 0){
						graph.addNode(ConBasePathGraph::ConPath::PosBase{pos, base},counter.chars_[base], counter.fractions_[base]);
					}
				}
			}

			//count up the pathways that seqs take and pick the path most traveled as the consensus path
			//ConPathCounter pathCounter;
			for(const auto & seq : reads_){
				ConBasePathGraph::ConPath currentPath(seq->seqBase_.cnt_);
				alignerObj.alignCacheGlobal(seqBase_, seq);
				for(const auto pos : importantPositions){
					currentPath.addPosBase(pos, alignerObj.alignObjectB_.seqBase_.seq_[alignerObj.getAlignPosForSeqAPos(pos)]);
				}

				bib::sort(importantPositions);
				for(const auto pos : iter::range(importantPositions.size() - 1)){
					auto headPos = importantPositions[pos];
					auto headBase = alignerObj.alignObjectB_.seqBase_.seq_[alignerObj.getAlignPosForSeqAPos(headPos)];
					auto tailPos = importantPositions[pos + 1];
					auto tailBase = alignerObj.alignObjectB_.seqBase_.seq_[alignerObj.getAlignPosForSeqAPos(tailPos)];
					graph.addEdge(ConBasePathGraph::ConPath::PosBase{headPos,headBase}.getUid(),ConBasePathGraph::ConPath::PosBase{tailPos,tailBase}.getUid(),seq->seqBase_.cnt_ );
				}

				//pathCounter.addConPath(currentPath);
			}
			/*
			std::ofstream outPathFile(firstReadName_ + "_basePaths.tab.txt");
			outPathFile << "Name: " << seqBase_.name_ << std::endl;
			outPathFile << "cnt: " << seqBase_.cnt_ << std::endl;
			outPathFile << "frac: " << seqBase_.frac_ << std::endl;
			graph.writePaths(outPathFile);
			*/
			auto paths = graph.getPaths();
			std::vector<ConBasePathGraph::ConPath> bestPaths;
			double bestCount = 0;

			for (const auto & path : paths) {
				if(path.count_ > bestCount){
					bestCount = path.count_;
					bestPaths.clear();
					bestPaths.emplace_back(path);
				}else if (path.count_ == bestCount){
					bestPaths.emplace_back(path);
				}
			}

			ConBasePathGraph::ConPath bestPath(0);
			if(bestPaths.size() > 1){

				double bestImprovement = std::numeric_limits<double>::lowest();
				std::vector<ConBasePathGraph::ConPath> secondaryBestPaths;
				for (const auto & path : bestPaths) {
					//outPathFile << path.getUid() << " : " << path.count_ << std::endl;
					double avgOppositeQual = 0;
					uint32_t baseCount = 0;
					for (const auto & base : path.bases_) {
						if('-' != base.base_){
							++baseCount;
							avgOppositeQual += counters.at(base.pos_).qualities_[base.base_]/static_cast<double>(counters.at(base.pos_).chars_[base.base_]) ;
						}
						//outPathFile << counters.at(base.pos_).qualities_[base.base_]/static_cast<double>(counters.at(base.pos_).chars_[base.base_]) << " ";
					}
					avgOppositeQual /=baseCount;
					//outPathFile << avgOppositeQual;
					//outPathFile << std::endl;
					double avgPresentQual = 0;

					for (const auto & base : path.bases_) {
						avgPresentQual+=seqBase_.qual_[base.pos_];
						//outPathFile << seqBase_.qual_[base.pos_] << " ";
					}
					avgPresentQual /= path.bases_.size();
					//outPathFile << avgPresentQual;
					//outPathFile << std::endl;
					//outPathFile << avgOppositeQual - avgPresentQual << std::endl;
					double qualImprovement = avgOppositeQual - avgPresentQual;

					if(qualImprovement > bestImprovement){
						bestImprovement = qualImprovement;
						secondaryBestPaths.clear();
						secondaryBestPaths.emplace_back(path);
					}else if(qualImprovement == bestImprovement){
						secondaryBestPaths.emplace_back(path);
					}
				}
				bestPath = secondaryBestPaths.front();
			}else if(bestPaths.size() == 1){
				bestPath = bestPaths.front();
			}else{
				std::stringstream ss;
				ss << __FILE__ << ":" << __LINE__ << " - " << __PRETTY_FUNCTION__ << "\n";
				ss << "Error, bestPaths should contain at least 1 path" << "\n";
				throw std::runtime_error{ss.str()};
			}

			//auto bestPath = pathCounter.getBestPath();
			//if there is just one supporting reads as the best path just do a majority's rule's consensus
			if(bestPath.count_ > 1){
				/*if("lib07_Minor.00_seq.0001_54-2_t1" == firstReadName_){
					std::cout << "lib07_Minor.00_seq.0001_54-2_t1" << std::endl;
					std::cout << "BestPath " << std::endl;
					std::cout << bestPath.toJson() << std::endl;
					std::cout << "All Paths" << std::endl;
					std::cout << pathCounter.toJson() << std::endl;
				}*/
				for(const auto & pb : bestPath.bases_){
					uint32_t otherCount = 0;
					auto & counter = counters.at(pb.pos_);
					for(const auto base : counter.alphabet_){
						if(base!= pb.base_){
							otherCount += counter.chars_[base];
							counter.chars_[base] = 0;
						}
					}
					//add the other's count so this becomes the majority
					//add qualities as well so that quality doesn't artificially drop
					double avgQual = counter.qualities_[pb.base_]/static_cast<double>(counter.chars_[pb.base_]);
					counter.qualities_[pb.base_] += std::round(avgQual * otherCount);
					counter.chars_[pb.base_] += otherCount;
					counter.setFractions();
				}
			}
		}
	}


	consensusHelper::genConsensusFromCounters(calcConsensusInfo_, counters, insertions,
			beginningGap);

	if (setToConsensus) {
		if (seqBase_.seq_ != calcConsensusInfo_.seq_) {
			seqBase_.seq_ = calcConsensusInfo_.seq_;
			setLetterCount();
		} else {
			seqBase_.seq_ = calcConsensusInfo_.seq_;
		}
		seqBase_.qual_ = calcConsensusInfo_.qual_;
	}
	previousErrorChecks_.clear();
	needToCalculateConsensus_ = false;
}

void baseCluster::calculateConsensus(aligner& alignerObj, bool setToConsensus) {
	// if the cluster is only one read, no need to create consensus
	if (reads_.size() <= 1) {
		return;
	}
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus_) {
    return;
  }
  auto averageSize = getAverageReadLength();
  seqInfo longestcluster;
  if (uAbsdiff(this->seqBase_.seq_.size(), averageSize) >
      0.1 * averageSize) {
    uint64_t biggest = 0;
    for (const auto& read : reads_) {
      if (read->seqBase_.seq_.length() > biggest) {
        longestcluster = read->seqBase_;
        biggest = read->seqBase_.seq_.length();
      }
    }
  } else {
    longestcluster = seqInfo(seqBase_.name_, seqBase_.seq_, seqBase_.qual_);
  }

  calculateConsensusTo(longestcluster, alignerObj, setToConsensus);

}


//////align the current clusters to the curent consensus
std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
baseCluster::calculateAlignmentsToConsensus(aligner& alignObj) {
  std::vector<baseReadObject> withConAlignments;
  std::vector<baseReadObject> withOutConAlignments;

  for (const auto & read : reads_) {
    alignObj.alignCacheGlobal(*this, read);

    baseReadObject tempConsensus = alignObj.alignObjectA_;
    tempConsensus.seqBase_.name_ = seqBase_.name_;
    tempConsensus.seqBase_.name_.append("_consensus");
    baseReadObject tempRead = alignObj.alignObjectB_;
    tempRead.seqBase_.name_ = read->seqBase_.name_;

    withConAlignments.push_back(tempRead);
    withConAlignments.push_back(tempConsensus);
    withOutConAlignments.push_back(tempRead);
  }
  return {withConAlignments, withOutConAlignments};
}
// write out the alignments to the consensus
void baseCluster::writeOutAlignments(const std::string& directoryName,
                                     aligner& alignObj) {
  std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
      alignments = calculateAlignmentsToConsensus(alignObj);
  std::stringstream alignMentsFileName;
  std::stringstream alignmentsOnlyFileName;
  alignMentsFileName << directoryName << "align_" << seqBase_.name_ << ".fasta";
  alignmentsOnlyFileName << directoryName << "noConAlign_" << seqBase_.name_
                         << ".fasta";

  SeqIOOptions alignMentsFileNameOpts = SeqIOOptions::genFastaOut(alignMentsFileName.str());
  SeqOutput firstWriter(alignMentsFileNameOpts);
  firstWriter.openWrite(alignments.first);
  SeqIOOptions alignmentsOnlyFileNameOpts = SeqIOOptions::genFastaOut(alignmentsOnlyFileName.str());;
  SeqOutput secondWriter(alignmentsOnlyFileNameOpts);
  secondWriter.openWrite(alignments.second);
}

void baseCluster::writeClustersInDir(const std::string& workingDir,
		const SeqIOOptions & ioOptions) const {
	SeqIOOptions options(workingDir + seqBase_.name_,ioOptions.outFormat_, ioOptions.out_);
	writeClusters(options);
}

/// output the clusters currently clustered to the this cluster
void baseCluster::writeClusters(const SeqIOOptions & ioOptions) const {
	SeqOutput writer(ioOptions);
	writer.openWrite(reads_);
}

VecStr baseCluster::getReadNames() const { return readVec::getNames(reads_); }

double baseCluster::getAverageReadLength() const {
  if (reads_.size() == 1) {
    return (double)seqBase_.seq_.length();
  } else {

    int sumOfLength = 0;
    int numberOfReads = 0;
    double averageReadLength = 0.0;
    for (const auto& read : reads_) {
      sumOfLength += read->seqBase_.seq_.length() * read->seqBase_.cnt_;
      numberOfReads += read->seqBase_.cnt_;
    }
    averageReadLength = (double)sumOfLength / numberOfReads;
    return averageReadLength;
  }
}



readObject baseCluster::createRead() const {
	return readObject(seqBase_);
}

std::string toSlimJsonErrors(const comparison & comp){
	Json::Value ret;
	ret["highQualityMatches_"] = comp.highQualityMatches_;
	ret["hqMismatches_"] = comp.hqMismatches_;
	ret["largeBaseIndel_"] = comp.largeBaseIndel_;
	ret["lowKmerMismatches_"] = comp.lowKmerMismatches_;
	ret["lowQualityMatches_"] = comp.lowQualityMatches_;
	ret["lqMismatches_"] = comp.lqMismatches_;
	ret["oneBaseIndel_"] = comp.oneBaseIndel_;
	ret["twoBaseIndel_"] = comp.twoBaseIndel_;
	Json::FastWriter jWriter;
	return jWriter.write(ret);
}


bool baseCluster::compare(baseCluster & read, aligner & alignerObj,
		const IterPar & runParams, const CollapserOpts & collapserOptsObj) {
	bool ret = false;
	if (previousErrorChecks_.find(read.firstReadName_) != previousErrorChecks_.end()
			&& read.previousErrorChecks_.find(firstReadName_)
					!= read.previousErrorChecks_.end()) {
		ret = runParams.passErrorCheck(
				read.previousErrorChecks_.at(firstReadName_));
		//if(read.previousErrorChecks_.at(firstReadName_).distances_.eventBasedIdentity_ > .95 and read.seqBase_.cnt_ == 1 and !ret){
		//	std::cout << toSlimJsonErrors(read.previousErrorChecks_.at(firstReadName_)) << std::endl;
		//}
	} else {
    if(collapserOptsObj.alignOpts_.noAlign_){
    	alignerObj.noAlignSetAndScore(seqBase_, read.seqBase_);
    }else{
    	alignerObj.alignCacheGlobal(seqBase_, read.seqBase_);
    }
		comparison currentProfile = alignerObj.compareAlignment(seqBase_, read.seqBase_,
				 collapserOptsObj.kmerOpts_.checkKmers_);
		if (currentProfile.distances_.query_.coverage_ < 0.50
				|| currentProfile.distances_.ref_.coverage_< 0.50) {
			ret = false;
		} else {
			read.previousErrorChecks_[firstReadName_] = currentProfile;
			previousErrorChecks_[read.firstReadName_] = currentProfile;
			//std::stringstream ss;
			//ss << bib::json::toJson(currentProfile) << std::endl;;
			//std::cout << bib::removeAllWhitespace(ss.str()) << std::endl;;
			//ss.str("");
			//ss << bib::json::toJson(runParams.errors_) << std::endl;
			//std::cout << bib::removeAllWhitespace(ss.str()) << std::endl;;
			ret = runParams.passErrorCheck(currentProfile);
			//if (currentProfile.distances_.eventBasedIdentity_ > .95
			//		&& read.seqBase_.cnt_ == 1 and !ret) {
			//	std::cout << toSlimJsonErrors(currentProfile) << std::endl;
			//}
			//std::cout << ret << std::endl;
		}
	}
	return ret;
}

bool baseCluster::isClusterCompletelyChimeric() {
  for (const auto &read : reads_) {
    if (read->seqBase_.name_.find("CHI") == std::string::npos) {
      return false;
    }
  }
  return true;
}

bool baseCluster::isClusterAtLeastHalfChimeric() {
	size_t chiCount = 0;
	uint32_t chiReadCnt = 0;
	double total = readVec::getTotalReadCount(reads_);
	for (const auto &read : reads_) {
		if (read->seqBase_.name_.find("CHI") != std::string::npos) {
			++chiCount;
			chiReadCnt += read->seqBase_.cnt_;
		}
	}
	if (chiReadCnt >= total / 2.0) {
		return true;
	}
	return false;
}

bool baseCluster::isClusterAtLeastChimericCutOff(double cutOff) {
	size_t chiCount = 0;
	uint32_t chiReadCnt = 0;
	double total = readVec::getTotalReadCount(reads_);
	for (const auto &read : reads_) {
		if (read->seqBase_.name_.find("CHI") != std::string::npos) {
			++chiCount;
			chiReadCnt += read->seqBase_.cnt_;
		}
	}
	if (chiReadCnt/total >= cutOff) {
		return true;
	}
	return false;
}


void baseCluster::removeReads(std::vector<uint32_t> readPositions){
	std::sort(readPositions.rbegin(), readPositions.rend());
	for(const auto & pos : readPositions){
		removeRead(pos);
	}
}

void baseCluster::removeRead(uint32_t readPos) {
	if (readPos >= reads_.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error readPos: " << readPos
				<< " out of range of reads size: " << reads_.size() << "\n";
		throw std::out_of_range { ss.str() };
	}
	seqBase_.cnt_ -= reads_[readPos]->seqBase_.cnt_;
	if (reads_[readPos]->seqBase_.name_ == firstReadName_) {
		if (0 != readPos) {
			firstReadName_ = reads_[0]->seqBase_.name_;
			firstReadCount_ = reads_[0]->seqBase_.cnt_;
		} else if (reads_.size() > 1) {
			firstReadName_ = reads_[1]->seqBase_.name_;
			firstReadCount_ = reads_[1]->seqBase_.cnt_;
		}else{
			firstReadName_ = "";
		}
	}
	reads_.erase(reads_.begin() + readPos);
	needToCalculateConsensus_ = true;
	updateName();
}


void baseCluster::removeRead(const std::string & stubName){
	uint32_t readPos = std::numeric_limits<uint32_t>::max();
	for(const auto & pos : iter::range(reads_.size())){
		if(reads_[pos]->getStubName(true) == stubName){
			readPos = pos;
			break;
		}
	}
	removeRead(readPos);
}

}  // namespace bib
