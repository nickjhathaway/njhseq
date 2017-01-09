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



class ConPath {
public:
	struct PosBase {
		uint32_t pos_;
		char base_;
		std::string getUid() const {
			return estd::to_string(pos_) + estd::to_string(base_);
		}

		Json::Value toJson() const {
			Json::Value ret;
			ret["class"] = bib::getTypeName(*this);
			ret["pos_"] = bib::json::toJson(pos_);
			ret["base_"] = bib::json::toJson(base_);
			return ret;
		}
	};

	ConPath(double count):count_(count){}

	double count_ = 1;
	std::vector<PosBase> bases_;
	void addPosBase(uint32_t pos, char base){
		addPosBase({pos, base});
	}
	void addPosBase(const PosBase & pb){
		bases_.emplace_back(pb);
	}

	std::string getUid()const{
		std::string ret = "";
		for(const auto & pb : bases_){
			ret += pb.getUid();
		}
		return ret;
	}

	Json::Value toJson() const {
		Json::Value ret;
		ret["class"] = bib::getTypeName(*this);
		ret["count_"] = bib::json::toJson(count_);
		ret["bases_"] = bib::json::toJson(bases_);
		return ret;
	}
};

class ConPathCounter {
public:

	std::unordered_map<std::string, ConPath> paths_;

	void addConPath(const ConPath & path) {
		auto search = paths_.find(path.getUid());
		if (paths_.end() != search) {
			search->second.count_ += path.count_;
		} else {
			paths_.emplace(path.getUid(), path);
		}
	}

	ConPath getBestPath() const {
		ConPath ret(0);
		for (const auto & path : paths_) {
			if (path.second.count_ > ret.count_) {
				ret = path.second;
			}
		}
		return ret;
	}

	Json::Value toJson() const {
		Json::Value ret;
		ret["class"] = bib::getTypeName(*this);
		ret["paths_"] = bib::json::toJson(paths_);
		return ret;
	}

};



class ConBaseOverlapGraph {
public:

	class edge;
	class node {
	public:

		node(const ConPath::PosBase & val,
				double cnt, double frac) :
				val_(val),
				cnt_(cnt),
				frac_(frac){
		}
		ConPath::PosBase val_;
		double cnt_;
		double frac_;

		std::vector<std::shared_ptr<edge>> headEdges_;
		std::vector<std::shared_ptr<edge>> tailEdges_;

		uint32_t visitCount_ = 0;

		void resetVisitCount();

		void addHead(const std::shared_ptr<edge> & e) {
			headEdges_.push_back(e);
		}

		void addTail(const std::shared_ptr<edge> & e) {
			tailEdges_.push_back(e);
		}

		bool headless() const {
			return headEdges_.empty();
		}

		bool tailless() const {
			return tailEdges_.empty();
		}

		void addToWritingPath(std::ostream & out,
				std::string currentPath) {
			++visitCount_;
			currentPath += val_.getUid();
			if (tailless()) {
				out << currentPath << std::endl;
			}
			for (const auto & tail : tailEdges_) {
				tail->tail_.lock()->addToWritingPath(out, currentPath + " - " + estd::to_string(tail->cnt_) + " > ");
			}
		};
		void addToPath(std::vector<ConPath> & paths,
				ConPath currentPath) {
			++visitCount_;
			currentPath.addPosBase(val_.pos_, val_.base_);
			if (tailless()) {
				paths.emplace_back(currentPath);
			}
			for (const auto & tail : tailEdges_) {
				auto tempPath = currentPath;
				tempPath.count_ += tail->cnt_;
				tail->tail_.lock()->addToPath(paths, tempPath);
			}
		};
	};
	class edge {
	public:
		edge(const std::shared_ptr<node> & head,
				const std::shared_ptr<node> & tail,
				double cnt) :
				head_(head), tail_(tail),
				cnt_(cnt){

		};
		std::weak_ptr<node> head_;
		std::weak_ptr<node> tail_;
		double cnt_;
	};
	std::unordered_map<std::string, std::shared_ptr<node>> nodes_;

	void addNode(const ConPath::PosBase & n, double cnt, double frac) {
		if (bib::has(nodes_, n.getUid())) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, alreay contains node with name "
					<< n.getUid() << "\n";
			throw std::runtime_error { ss.str() };
		}
		nodes_.emplace(n.getUid(), std::make_shared<node>(n, cnt, frac));
	}

	void addEdge(const std::string & head,
			const std::string & tailName,
			double cnt) {
		/**@todo add way to check if edge already exists */
		auto headNode = nodes_.at(head);
		bool foundEdge = false;
		for(const auto & tail : headNode->tailEdges_){
			auto tailNode = tail->tail_.lock();
			if(tailName == tailNode->val_.getUid() ){
				foundEdge = true;
				tail->cnt_ += cnt;
				break;
			}
		}
		if(!foundEdge){
			auto tailNode = nodes_.at(tailName);
			std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode, cnt);
			headNode->addTail(e);
			tailNode->addHead(e);
		}
	}

	void writePaths(std::ostream & out) const {
		/**@todo add way to check if there are any cycles or no headless nodes */
		for (const auto & n : nodes_) {
			if (n.second->headless()) {
				n.second->addToWritingPath(out, "");
			}
		}
		VecStr notVisitedNodes;
		for (const auto & n : nodes_) {
			if (n.second->visitCount_ == 0) {
				notVisitedNodes.emplace_back(n.first);
			}
		}
		if (!notVisitedNodes.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ": Error, the following nodes weren't visited:\n";
			ss << bib::conToStr(notVisitedNodes, ",") << "\n";
			throw std::runtime_error { ss.str() };
		}
	}

	std::vector<ConPath> getPaths() const {
		std::vector<ConPath> ret;
		/**@todo add way to check if there are any cycles or no headless nodes */
		for (const auto & n : nodes_) {
			if (n.second->headless()) {
				n.second->addToPath(ret, ConPath{0});
			}
		}
		VecStr notVisitedNodes;
		for (const auto & n : nodes_) {
			if (n.second->visitCount_ == 0) {
				notVisitedNodes.emplace_back(n.first);
			}
		}
		if (!notVisitedNodes.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ": Error, the following nodes weren't visited:\n";
			ss << bib::conToStr(notVisitedNodes, ",") << "\n";
			throw std::runtime_error { ss.str() };
		}
		return ret;
	}

	Json::Value createSankeyOutput() const {
		Json::Value ret;
		std::vector<std::shared_ptr<node>> nodesVec;
		std::unordered_map<std::string, uint32_t> nodePosition;
		for (const auto & n : nodes_) {
			nodePosition[n.first] = nodesVec.size();
			nodesVec.push_back(n.second);
		}
		auto &nodes = ret["nodes"];
		auto &links = ret["links"];

		for (const auto & n : nodesVec) {
			Json::Value nodeJson;
			nodeJson["name"] = n->val_.getUid();
			nodeJson["cnt"] = n->cnt_;
			nodeJson["frac"] = n->frac_;
			nodes.append(nodeJson);
			double totalTail = 0;
			for (const auto & tl : n->tailEdges_) {
				totalTail += tl->cnt_;
			}
			for (const auto & tl : n->tailEdges_) {
				Json::Value linkJsons;
				linkJsons["source"] = nodePosition[tl->head_.lock()->val_.getUid()];
				linkJsons["target"] = nodePosition[tl->tail_.lock()->val_.getUid()];
				;
				linkJsons["value"] = tl->cnt_ / totalTail;
				links.append(linkJsons);
			}
		}
		return ret;
	}


};



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
			ConBaseOverlapGraph graph;
			for(const auto pos : importantPositions){
				auto & counter = counters.at(pos);
				for(const auto base : counter.alphabet_){
					if(counter.chars_[base] > 0){
						graph.addNode(ConPath::PosBase{pos, base},counter.chars_[base], counter.fractions_[base]);
					}
				}
			}
			//count up the pathways that seqs take and pick the path most traveled as the consensus path
			//ConPathCounter pathCounter;
			for(const auto & seq : reads_){
				ConPath currentPath(seq->seqBase_.cnt_);
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
					graph.addEdge(ConPath::PosBase{headPos,headBase}.getUid(),ConPath::PosBase{tailPos,tailBase}.getUid(),seq->seqBase_.cnt_ );
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
			std::vector<ConPath> bestPaths;
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
			ConPath bestPath(0);
			if(bestPaths.size() > 1){
				double bestImprovement = std::numeric_limits<double>::lowest();
				std::vector<ConPath> secondaryBestPaths;
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

}  // namespace bib
