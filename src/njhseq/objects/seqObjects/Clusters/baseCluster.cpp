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

#include "baseCluster.hpp"
#include "njhseq/helpers/consensusHelper.hpp"
#include "njhseq/readVectorManipulation/readVectorHelpers.h"
#include "njhseq/IO/SeqIO/SeqOutput.hpp"
#include "njhseq/objects/dataContainers/graphs/ConBasePathGraph.hpp"


namespace njhseq {


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
void baseCluster::calculateConsensusToCurrent(aligner& alignerObj, calculateConsensusPars conPars){
	// if the cluster is only one read, no need to create consensus
	if (reads_.size() <= 1) {
		return;
	}
  //check to see if a consensus needs to be built
  if (!needToCalculateConsensus_) {
    return;
  }
  calculateConsensusTo(seqBase_, alignerObj, conPars);

}











bool baseCluster::calculateConsensusTo(
		const seqInfo & calcToThisInput,
		aligner& alignerObj,
		calculateConsensusPars conPars) {
	bool matchCurrentSeq = false;
	uint32_t convergeCount = 0;
	seqInfo calcToThisCopy = calcToThisInput;

//	for(const auto & clus : reads_){
//		std::cout << clus->seqBase_.name_ << std::endl;
//		std::cout << "\t" << clus->seqBase_.cnt_ << std::endl;
//	}

	while(!matchCurrentSeq && convergeCount <= conPars.convergeAttempts){
		++convergeCount;
		//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//	std::cout << njh::bashCT::boldRed(seqBase.name_) << std::endl;
		//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
		// create the map for letter counters for each position
		std::map<uint32_t, charCounter> counters;
		// create a map in case of insertions
		std::map<uint32_t, std::map<uint32_t, charCounter>> insertions;
		std::map<int32_t, charCounter> beginningGap;
		auto getSeqBase = [](const std::shared_ptr<readObject> & read) ->const seqInfo& {return read->seqBase_;};
//		std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
		//consensusHelper::increaseCounters(seqBase_, reads_, getSeqBase, alignerObj, counters, insertions, beginningGap);

		if(noWeightConsensus_){

			consensusHelper::increaseCountersNoReadWeights(calcToThisCopy, reads_, getSeqBase, alignerObj, counters, insertions, beginningGap);
		}else{
			consensusHelper::increaseCounters(calcToThisCopy, reads_, getSeqBase, alignerObj, counters, insertions, beginningGap);
		}

//		std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
		if(noWeightConsensus_){
			calcConsensusInfo_.cnt_ = reads_.size();
			calcConsensusInfo_.frac_ = seqBase_.frac_; //??? not sure about this
		}else{
			calcConsensusInfo_.cnt_ = seqBase_.cnt_;
			calcConsensusInfo_.frac_ = seqBase_.frac_;
		}

		calcConsensusInfo_.name_ = calcToThisCopy.name_;
		calcConsensusInfo_.seq_ = calcToThisCopy.seq_;
		calcConsensusInfo_.qual_ = calcToThisCopy.qual_;

		//calcConsensusInfo_ = calcToThisCopy;

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		seqBase_.outPutFastq(std::cout);
//		calcToThisCopy.outPutFastq(std::cout);

//		std::cout << "count: " << seqBase_.cnt_ << std::endl;
	//	consensusHelper::increaseCounters(seqBase, reads_, getSeqBase, alignerObj, counters, insertions, beginningGap);
	//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//	calcConsensusInfo_ = seqBase;
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
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//if the count is just 2 then just a majority rules consensus
		if(reads_.size() > 2){
		//if(seqBase_.cnt_ > 2){
			/*//for debugging
			bool print = false;
			if(njh::containsSubString(seqBase_.name_, "lib1_Minor.00_seq.0001_5")){
				print = true;
			}*/

			//find out if there are several locations with larger contention of majority rules
			//consensus and therefore could lead to bad consensus building
			double contentionCutOff = 0.25;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			while(importantPositions.size() > 20){
				contentionCutOff += .10;
				countAbovepCutOff = 0;
				importantPositions.clear();
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
			}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "countAbovepCutOff: " << countAbovepCutOff << std::endl;

			//for debugging
			/*if(print){
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				std::cout << seqBase_.name_ << std::endl;
				std::cout << "countAbovepCutOff: " << countAbovepCutOff << std::endl;
				bfs::path counterFnp = firstReadName_ + "_baseCounts.tab.txt";
				counterFnp = njh::files::findNonexitantFile(counterFnp.string());
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
				insert_counterFnp = njh::files::findNonexitantFile(insert_counterFnp.string());
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
					{
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				std::cout << seqBase_.name_ << std::endl;
				std::cout << firstReadName_ << std::endl;
				std::cout << "countAbovepCutOff: " << countAbovepCutOff << std::endl;
				bfs::path counterFnp = firstReadName_ + "_baseCounts.tab.txt";
				counterFnp = njh::files::findNonexitantFile(counterFnp.string());
				std::ofstream outFile(njh::replaceString(counterFnp.string(), "/", "_"));

				outFile << "pos\tbase\tcount\tfrac" << std::endl;
				for (const auto & counter : counters) {
					if(njh::in(counter.first, importantPositions)){
						for (const auto base : counter.second.alphabet_) {
							if (counter.second.fractions_[base] > contentionCutOff) {
								outFile << counter.first << "\t" << base << "\t"
																<< counter.second.chars_[base] << "\t"
																<< counter.second.fractions_[base] << std::endl;;
							}
						}
					}
				}
				bfs::path insert_counterFnp = firstReadName_ + "_insertBaseCounts.tab.txt";
				insert_counterFnp = njh::files::findNonexitantFile(insert_counterFnp.string());
				std::ofstream outInsertFile(njh::replaceString(insert_counterFnp.string(), "/", "_"));
				outFile << "pos\tinsertPos\tbase\tcount\tfrac" << std::endl;
				for (const auto & insert : insertions) {
					for (const auto & counter : insert.second) {
						for (const auto base : counter.second.alphabet_) {
							if(njh::in(insert.first, importantPositions)){
								outInsertFile << insert.first << "\t"
										<< counter.first << "\t" << base << "\t"
										<< counter.second.chars_[base] << "\t"
										<< counter.second.fractions_[base] << std::endl;;
							}
						}
					}
				}
				outFile.close();
				outInsertFile.close();
			}
			}*/


			//if there are several points of contention
			if(countAbovepCutOff >= 2){
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				//for debugging;
				/*
				if(print){
					//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				}

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
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				ConBasePathGraph graph;
				for(const auto pos : importantPositions){
					auto & counter = counters.at(pos);
					for(const auto base : counter.alphabet_){
						if (counter.fractions_[base] > contentionCutOff) {
							graph.addNode(ConBasePathGraph::ConPath::PosBase { pos, base },
									counter.chars_[base], counter.fractions_[base]);
						}
						/*if(counter.chars_[base] > 0){
							graph.addNode(ConBasePathGraph::ConPath::PosBase{pos, base},counter.chars_[base], counter.fractions_[base]);
						}*/
					}
				}
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				//count up the pathways that seqs take and pick the path most traveled as the consensus path
				for(const auto & seq : reads_){
					//alignerObj.alignCacheGlobal(seqBase_, seq);
					//alignerObj.alignCacheGlobalDiag(seqBase_, seq);
					alignerObj.alignCacheGlobalDiag(calcToThisCopy, seq);

					njh::sort(importantPositions);
					for(const auto pos : iter::range(importantPositions.size() - 1)){
						auto headPos = importantPositions[pos];
						auto headBase = alignerObj.alignObjectB_.seqBase_.seq_[alignerObj.getAlignPosForSeqAPos(headPos)];
						auto tailPos = importantPositions[pos + 1];
						auto tailBase = alignerObj.alignObjectB_.seqBase_.seq_[alignerObj.getAlignPosForSeqAPos(tailPos)];
						if(tailPos <= headPos){
							std::stringstream ss;
							ss << "tailPos: " << tailPos << " is less than headPos: " << headPos << std::endl;
							throw std::runtime_error{ss.str()};
						}
						if(counters.at(headPos).fractions_[headBase] > contentionCutOff &&
								counters.at(tailPos).fractions_[tailBase] > contentionCutOff ){
							if(noWeightConsensus_){
								graph.addEdge(ConBasePathGraph::ConPath::PosBase{headPos,headBase}.getUid(),ConBasePathGraph::ConPath::PosBase{tailPos,tailBase}.getUid(),1 );
								//here
							}else{
								graph.addEdge(ConBasePathGraph::ConPath::PosBase{headPos,headBase}.getUid(),ConBasePathGraph::ConPath::PosBase{tailPos,tailBase}.getUid(),seq->seqBase_.cnt_ );
							}
						}
						//graph.addEdge(ConBasePathGraph::ConPath::PosBase{headPos,headBase}.getUid(),ConBasePathGraph::ConPath::PosBase{tailPos,tailBase}.getUid(),seq->seqBase_.cnt_ );
					}
				}
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				//for debugging;
				/*
				if(print){
					std::ofstream outPathFile(firstReadName_ + "_basePaths.tab.txt");
					outPathFile << "Name: " << seqBase_.name_ << std::endl;
					outPathFile << "cnt: " << seqBase_.cnt_ << std::endl;
					outPathFile << "frac: " << seqBase_.frac_ << std::endl;
					graph.writePaths(outPathFile);
				}

				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				{
					std::ofstream outJson("conPaths.json");
					outJson << graph.createSankeyOutput() << std::endl;;
				}*/
				auto paths = graph.getPaths();
				std::vector<ConBasePathGraph::ConPath> bestPaths;
				double bestCount = 0;
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				for (const auto & path : paths) {
					if(path.count_ > bestCount){
						bestCount = path.count_;
						bestPaths.clear();
						bestPaths.emplace_back(path);
					}else if (path.count_ == bestCount){
						bestPaths.emplace_back(path);
					}
				}
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				ConBasePathGraph::ConPath bestPath(0);
				if(bestPaths.size() > 1){
					//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
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
							//avgPresentQual+=seqBase_.qual_[base.pos_];
							avgPresentQual+=calcToThisCopy.qual_[base.pos_];
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
				//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
				//auto bestPath = pathCounter.getBestPath();
				//if there is just one supporting reads as the best path just do a majority's rule's consensus
				if(bestPath.count_ > 1){
					//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
					/*if("lib07_Minor.00_seq.0001_54-2_t1" == firstReadName_){
						std::cout << "lib07_Minor.00_seq.0001_54-2_t1" << std::endl;
						std::cout << "BestPath " << std::endl;
						std::cout << bestPath.toJson() << std::endl;
						std::cout << "All Paths" << std::endl;
						std::cout << pathCounter.toJson() << std::endl;
					}*/
					for(const auto & pb : bestPath.bases_){
						double otherCount = 0;
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
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "noWeightConsensus_: " << njh::colorBool(noWeightConsensus_) << std::endl;
//		std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
//    calcConsensusInfo_.outPutSeqAnsi(std::cout);
		consensusHelper::genConsensusFromCounters(calcConsensusInfo_, counters, insertions, beginningGap);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		matchCurrentSeq = seqBase_.seq_ == calcConsensusInfo_.seq_;
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(calcConsensusInfo_.seq_.size() >= alignerObj.parts_.maxSize_){
			alignerObj.parts_.setMaxSize(calcConsensusInfo_.seq_.size() + 20 );
		}
		if (conPars.setToConsensus) {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			if (seqBase_.seq_ != calcConsensusInfo_.seq_) {
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				seqBase_.seq_ = calcConsensusInfo_.seq_;
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				setLetterCount();
			}
			seqBase_.qual_ = calcConsensusInfo_.qual_;
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << calcToThisCopy.name_ << std::endl;
//		std::cout << seqBase_.name_ << std::endl;
//		std::cout << "\tconvergeCount  : " << convergeCount << std::endl;
//		std::cout << "\tmatchCurrentSeq: " << njh::colorBool(matchCurrentSeq) << std::endl;
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(!matchCurrentSeq && conPars.convergeConsensus && convergeCount <= conPars.convergeAttempts){
//			calcConsensusInfo_.outPutFastq(std::cout);
//			seqBase_.outPutFastq(std::cout);
			calcToThisCopy.seq_ = calcConsensusInfo_.seq_;
			calcToThisCopy.qual_ = calcConsensusInfo_.qual_;

//			calcToThisCopy.outPutFastq(std::cout);
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		{
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			auto debugOpts = SeqIOOptions::genFastqOut("test.fastq");
//			debugOpts.out_.append_ = true;
//			SeqOutput writter(debugOpts);
//			auto outDebug = calcConsensusInfo_;
//			outDebug.name_ = njh::pasteAsStr(convergeCount);
//			writter.openWrite(outDebug);
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}


		//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
		if(!conPars.convergeConsensus){
			break;
		}
	}
	previousErrorChecks_.clear();
	needToCalculateConsensus_ = false;
	if(reads_.size() > 2){
		//exit(1);
	}
	return matchCurrentSeq;
}

void baseCluster::calculateConsensus(aligner& alignerObj, calculateConsensusPars conPars) {
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

  calculateConsensusTo(longestcluster, alignerObj, conPars);

}


//////align the current clusters to the curent consensus
std::pair<std::vector<baseReadObject>, std::vector<baseReadObject>>
baseCluster::calculateAlignmentsToConsensus(aligner& alignObj) {
  std::vector<baseReadObject> withConAlignments;
  std::vector<baseReadObject> withOutConAlignments;

  for (const auto & read : reads_) {
    //alignObj.alignCacheGlobal(*this, read);
  		alignObj.alignCacheGlobalDiag(*this, read);

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
	return njh::json::writeAsOneLine(ret);
}


comparison baseCluster::getComparison(baseCluster & read, aligner & alignerObj, bool checkKmers) const{
	alignerObj.alignCacheGlobalDiag(seqBase_, read.seqBase_);
	alignerObj.compareAlignment(seqBase_, read.seqBase_, checkKmers);
	return alignerObj.comp_;
}

bool baseCluster::compare(baseCluster & read, aligner & alignerObj,
		const IterPar & runParams, const CollapserOpts & collapserOptsObj) {
	bool ret = false;
	if (previousErrorChecks_.find(read.firstReadName_) != previousErrorChecks_.end()
			&& read.previousErrorChecks_.find(firstReadName_)
			!= read.previousErrorChecks_.end()) {
//
//		if("[hap=Chi.0[PCRRound=1];idNum=03;mutated=false]" == firstReadName_ && read.seqBase_.cnt_ > 1000){
//			std::cout << njh::bashCT::blue;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "runParams: " << runParams.iterNumber_ << std::endl;
//			std::cout << "\tseqBase_.name_     : " << seqBase_.name_ << std::endl;
//			std::cout << "\tfirstReadName_     : " << firstReadName_ << std::endl;
//			std::cout << "\tread.seqBase_.name_: " << read.seqBase_.name_ << std::endl;
//			std::cout << "\tread.firstReadName_: " << read.firstReadName_ << std::endl;
//			std::cout << "\t" << toSlimJsonErrors(read.previousErrorChecks_.at(firstReadName_)) << std::endl;
//			std::cout << "\tqueryName_:" << read.previousErrorChecks_.at(firstReadName_).queryName_ << std::endl;
//			std::cout << "\trefName_  :" << read.previousErrorChecks_.at(firstReadName_).refName_ << std::endl;
//
//			std::cout << njh::bashCT::reset;
//		}
//		if("[hap=Chi.1[PCRRound=1];idNum=07;mutated=false]" == firstReadName_ && read.seqBase_.cnt_ > 1000){
//			std::cout << njh::bashCT::cyan;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "runParams: " << runParams.iterNumber_ << std::endl;
//			std::cout << "\tseqBase_.name_     : " << seqBase_.name_ << std::endl;
//			std::cout << "\tfirstReadName_     : " << firstReadName_ << std::endl;
//			std::cout << "\tread.seqBase_.name_: " << read.seqBase_.name_ << std::endl;
//			std::cout << "\tread.firstReadName_: " << read.firstReadName_ << std::endl;
//			std::cout << "\t" << toSlimJsonErrors(read.previousErrorChecks_.at(firstReadName_)) << std::endl;
//			std::cout << "\tqueryName_:" << read.previousErrorChecks_.at(firstReadName_).queryName_ << std::endl;
//			std::cout << "\trefName_  :" << read.previousErrorChecks_.at(firstReadName_).refName_ << std::endl;
//			std::cout << njh::bashCT::reset;
//		}
		ret = runParams.passErrorCheck(read.previousErrorChecks_.at(firstReadName_));
		//if(read.previousErrorChecks_.at(firstReadName_).distances_.eventBasedIdentity_ > .95 and read.seqBase_.cnt_ == 1 and !ret){
		//	std::cout << toSlimJsonErrors(read.previousErrorChecks_.at(firstReadName_)) << std::endl;
		//}
	} else {
    if(collapserOptsObj.alignOpts_.noAlign_){
    	alignerObj.noAlignSetAndScore(seqBase_, read.seqBase_);
    }else{
    //	alignerObj.alignCacheGlobal(seqBase_, read.seqBase_);
    	alignerObj.alignCacheGlobalDiag(seqBase_, read.seqBase_);
    }
		comparison currentProfile = alignerObj.compareAlignment(seqBase_, read.seqBase_,
				 collapserOptsObj.kmerOpts_.checkKmers_);
//		if("[hap=Chi.0[PCRRound=1];idNum=03;mutated=false]" == firstReadName_ && read.seqBase_.cnt_ > 1000){
//			std::cout << njh::bashCT::green;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "runParams: " << runParams.iterNumber_ << std::endl;
//			std::cout << "\tseqBase_.name_     : " << seqBase_.name_ << std::endl;
//			std::cout << "\tfirstReadName_     : " << firstReadName_ << std::endl;
//			std::cout << "\tread.seqBase_.name_: " << read.seqBase_.name_ << std::endl;
//			std::cout << "\tread.firstReadName_: " << read.firstReadName_ << std::endl;
//			std::cout << "\t" << toSlimJsonErrors(currentProfile) << std::endl;
//			std::cout << njh::bashCT::reset;
//		}
//		if("[hap=Chi.1[PCRRound=1];idNum=07;mutated=false]" == firstReadName_ && read.seqBase_.cnt_ > 1000){
//			std::cout << njh::bashCT::purple;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "runParams: " << runParams.iterNumber_ << std::endl;
//			std::cout << "\tseqBase_.name_     : " << seqBase_.name_ << std::endl;
//			std::cout << "\tfirstReadName_     : " << firstReadName_ << std::endl;
//			std::cout << "\tread.seqBase_.name_: " << read.seqBase_.name_ << std::endl;
//			std::cout << "\tread.firstReadName_: " << read.firstReadName_ << std::endl;
//			std::cout << "\t" << toSlimJsonErrors(currentProfile) << std::endl;
//			std::cout << njh::bashCT::reset;
//		}
		if (currentProfile.distances_.query_.coverage_ < 0.50
				|| currentProfile.distances_.ref_.coverage_< 0.50) {
			ret = false;
		} else {
			read.previousErrorChecks_[firstReadName_] = currentProfile;
			previousErrorChecks_[read.firstReadName_] = currentProfile;
			//std::stringstream ss;
			//ss << njh::json::toJson(currentProfile) << std::endl;;
			//std::cout << njh::removeAllWhitespace(ss.str()) << std::endl;;
			//ss.str("");
			//ss << njh::json::toJson(runParams.errors_) << std::endl;
			//std::cout << njh::removeAllWhitespace(ss.str()) << std::endl;;
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
    //if (read->seqBase_.name_.find("CHI") == std::string::npos) {
  	if(read->seqBase_.isChimeric()){
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
		//if (read->seqBase_.name_.find("CHI") != std::string::npos) {
		if(read->seqBase_.isChimeric()){
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
		//if (read->seqBase_.name_.find("CHI") != std::string::npos) {
		if(read->seqBase_.isChimeric()){
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
	for(const auto pos : iter::range(reads_.size())){
		if(reads_[pos]->getStubName(true) == stubName){
			readPos = pos;
			break;
		}
	}
	removeRead(readPos);
}

}  // namespace njh
