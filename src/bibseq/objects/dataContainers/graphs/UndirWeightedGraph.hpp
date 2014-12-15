//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#pragma once
/*
 * UndirWeightedGraph.hpp
 *
 *  Created on: Dec 10, 2014
 *      Author: nickhathaway
 */


#include "bibseq/common.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
#include <bibcpp/jsonUtils.h>
#include <bibcpp/graphics.h>


namespace bibseq {


template<typename DIST, typename VALUE>
class njhUndirWeightedGraph {
public:
	class node;

	class edge{
	public:

		edge(const DIST &dist,
				const std::shared_ptr<node> & node1,
				const std::shared_ptr<node> & node2): dist_(dist),
		on_(true), visited_(false){
			nodeToNode_[node1->name_] = node2;
			nodeToNode_[node2->name_] = node1;
		}

		DIST dist_;
		std::map<std::string, std::weak_ptr<node>> nodeToNode_;
		bool on_;
		bool visited_;
		bool best_ = false;
		uint32_t visitedAmount_ = 0;

		bool operator <(const edge & otherEdge){
			return dist_ < otherEdge.dist_;
		}
		bool operator >(const edge & otherEdge){
			return dist_ > otherEdge.dist_;
		}
	}; // class edge

	class node {
	public:

		node(const std::string & name,
				const VALUE & value):name_(name), value_(value),
				on_(true), visited_(false){

		}
		std::string name_;
		VALUE value_;
		std::vector<std::shared_ptr<edge>> edges_;
		bool on_;
		bool visited_;
		uint32_t visitedAmount_ = 0;
		uint32_t group_ = std::numeric_limits<uint32_t>::max();


		void sortEdges(){
			auto greaterComp = [](const std::shared_ptr<edge> & e1,
					const std::shared_ptr<edge> & e2){
				return *e1 > *e2;
			};
			sortEdges(greaterComp);
		}
		template<typename COMP>
		void sortEdges(COMP comp){
			std::sort(edges_.begin(), edges_.end(), comp);
		}

		void printEdges(std::ostream & out, bool printVisited){
			for(const auto & e : edges_){
				if(e->on_){
					if(e->visited_ && printVisited){
						out << name_ << " -- " << e->nodeToNode_[name_].lock()->name_  << " " << e->dist_ << "\n";
					}else{
						out << name_ << " -- " << e->nodeToNode_[name_].lock()->name_  << " " << e->dist_ << "\n";
					}
				}
			}
		}

		void printFirstEdge(std::ostream & out, bool visitEdge){
			if(!edges_.empty()){
				for(const auto & e : edges_){
					if(e->on_){
						if(!e->visited_){
							out << name_ << " -- " << e->nodeToNode_[name_].lock()->name_
									<< " " << e->dist_ << "\n";
						}
						if(visitEdge){
							e->visited_ = true;
						}
						break;
					}
				}
			}
		}

		void printLastEdge(std::ostream & out, bool visitEdge){
			if(!edges_.empty()){
				for(const auto & e : iter::reverse(edges_)){
					if(e->on_){
						if(!e->visited_){
							out << name_ << " -- " << e->nodeToNode_[name_].lock()->name_
									<< " " << e->dist_ << "\n";
						}
						if(visitEdge){
							e->visited_ = true;
						}
						break;
					}
				}
			}
		}

		void visitAndGroupCons(uint32_t currentGroup){
			if(!visited_){
				group_ = currentGroup;
				visited_ = true;
				++visitedAmount_;
				for(auto & e : edges_){
					if(e->on_ && !e->visited_){
						e->visited_ = true;
						e->nodeToNode_[name_].lock()->visitAndGroupCons(currentGroup);
					}
				}
			}
		}

		void determineHighestDistance(bool doTies){
			if(!edges_.empty()){
				DIST best = edges_.front()->dist_;
				std::vector<uint64_t> bestIndex = {std::numeric_limits<uint64_t>::max()};
				for(const auto & e : iter::enumerate(edges_)){
					if(e.element->on_){
						if(e.element->dist_ > best){
							best = e.element->dist_;
							bestIndex.clear();
							bestIndex = {e.index};
						} else if (doTies && e.element->dist_  == best){
							bestIndex.emplace_back(e.index);
						}
					}
				}
				if(bestIndex.front() != std::numeric_limits<uint64_t>::max()){
					if(doTies){
						for(const auto & b : bestIndex){
							edges_[b]->best_ = true;
						}
					}else{
						edges_[bestIndex.front()]->best_ = true;
					}
				}
			}
		}

		void determineLowestDistance(bool doTies){
			if(!edges_.empty()){
				DIST best = edges_.front()->dist_;
				std::vector<uint64_t> bestIndex = {std::numeric_limits<uint64_t>::max()};
				for(const auto & e : iter::enumerate(edges_)){
					if(e.element->on_){
						if(e.element->dist_ < best){
							best = e.element->dist_;
							bestIndex.clear();
							bestIndex = {e.index};
						} else if (doTies && e.element->dist_  == best){
							bestIndex.emplace_back(e.index);
						}
					}
				}
				if(bestIndex.front() != std::numeric_limits<uint64_t>::max()){
					if(doTies){
						for(const auto & b : bestIndex){
							edges_[b]->best_ = true;
						}
					}else{
						edges_[bestIndex.front()]->best_ = true;
					}
				}
			}
		}


	}; // class node



	std::vector<std::shared_ptr<node>> nodes_;
	std::vector<std::shared_ptr<edge>> edges_;
	std::unordered_map<std::string, uint64_t> nameToNodePos_;
	uint32_t numberOfGroups_ = 1;

	void addNode(const std::string & name, const VALUE & value){
		nameToNodePos_[name] = nodes_.size();
		nodes_.emplace_back(std::make_shared<node>(name, value));
	}

	void addEdge(const std::string & name1,
			const std::string & name2, const DIST & dist){
		edges_.emplace_back(std::make_shared<edge>(dist,
				nodes_[nameToNodePos_[name1]], nodes_[nameToNodePos_[name2]]));
		nodes_[nameToNodePos_[name1]]->edges_.emplace_back(edges_.back());
		nodes_[nameToNodePos_[name2]]->edges_.emplace_back(edges_.back());
	}


	void turnOffAllCons(){
		for(const auto & e : edges_){
			e->on_ = false;
		}
	}

	void turnOnAllCons(){
		for(const auto & e : edges_){
			e->on_ = true;
		}
	}

	void resetBestAndVistEdges(){
		for(const auto & e : edges_){
			e->visited_ = false;
			e->visitedAmount_ = 0;
			e->best_ = false;
		}
	}

	void resetVisitedEdges(){
		for(const auto & e : edges_){
			e->visited_ = false;
			e->visitedAmount_ = 0;
		}
	}

	void resetBestEdges(){
		for(const auto & e : edges_){
			e->best_ = false;
		}
	}

	void turnOffAllNodes(){
		for(const auto & n : nodes_){
			n->on_ = false;
		}
	}

	void turnOnAllNodes(){
		for(const auto & n : nodes_){
			n->on_ = true;
		}
	}

	void resetVisitedNodes(){
		for(const auto & n : nodes_){
			n->visited_ = false;
			n->visitedAmount_ = 0;
		}
	}


	void turnOnEdgesUnder(const DIST & cutOff){
		for(const auto & c : edges_){
			if(c->dist_ < cutOff){
				c->on_ = true;
			}else{
				c->on_ = false;
			}
		}
	}

	void turnOffEdgesUnder(const DIST & cutOff){
		for(const auto & c : edges_){
			if(c->dist_ < cutOff){
				c->on_ = false;
			}else{
				c->on_ = true;
			}
		}
	}

	void turnOnEdgesAbove(const DIST & cutOff){
		for(const auto & c : edges_){
			if(c->dist_ > cutOff){
				c->on_ = true;
			}else{
				c->on_ = false;
			}
		}
	}

	void turnOffEdgesAbove(const DIST & cutOff){
		for(const auto & c : edges_){
			if(c->dist_ > cutOff){
				c->on_ = false;
			}else{
				c->on_ = true;
			}
		}
	}

	void allSortEdges(){
		for(const auto & n : nodes_){
			n->sortEdges();
		}
	}
	template<typename COMP>
	void allSortEdges(COMP comp){
		for(const auto & n : nodes_){
			n->sortEdges(comp);
		}
	}

	void printAdj(std::ostream & out){
		for(const auto & e : edges_){
			if(e->on_){
				out << e->nodeToNode_.begin()->first << " -- "
						<< e->nodeToNode_.begin()->second.lock()->name_ << " "
						<< e->dist_ << "\n";
			}
		}
	}

	void allPrintEdges(std::ostream & out, bool printVisited){
		for(const auto & n : nodes_){
			n->printEdges(out, printVisited);
		}
	}

	void allPrintFirstEdge(std::ostream & out, bool visitEdge){
		for(const auto & n : nodes_){
			n->printFirstEdge(out, visitEdge);
		}
	}

	void allPrintLastEdge(std::ostream & out, bool visitEdge){
		for(const auto & n : nodes_){
			n->printLastEdge(out, visitEdge);
		}
	}

	void determineGroups(){
		numberOfGroups_ = 0;
		for(auto & n : nodes_){
			if(!n->visited_){
				n->visitAndGroupCons(numberOfGroups_);
				++numberOfGroups_;
			}
		}
	}



	void allDetermineHigestBest(bool doTies){
		for(const auto & n : nodes_){
			n->determineHighestDistance(doTies);
		}
		for(const auto & e : edges_){
			if(e->on_){
				if(!e->best_){
					e->on_ = false;
				}
			}
		}
	}

	void allDetermineLowestBest(bool doTies){
		for(const auto & n : nodes_){
			n->determineLowestDistance(doTies);
		}
		for(const auto & e : edges_){
			if(e->on_){
				if(!e->best_){
					e->on_ = false;
				}
			}
		}
	}

	std::map<uint32_t, std::uint32_t> getGroupCounts(){
	  std::map<uint32_t, uint32_t> groupCounts;
	  for(const auto & n : nodes_){
	  	++groupCounts[n->group_];
	  }
	  return groupCounts;
	}

	std::map<uint32_t, std::vector<VALUE>> getGroupValues(){
		std::map<uint32_t, std::vector<VALUE>> ret;
	  for(const auto & n : nodes_){
	  	ret[n->group_].emplace_back(n->value_);
	  }
	  return ret;
	}

	Json::Value toJson(uint32_t groupSizeCutOff){
	  Json::Value graphJson;
	  auto & nodes = graphJson["nodes"];
	  auto & links = graphJson["links"];
	  uint32_t nCount = 0;
	  std::unordered_map<uint32_t, uint32_t> groupCounts;
	  std::unordered_map<std::string, uint64_t> nameToNewPos;
	  for(const auto & n : nodes_){
	  	++groupCounts[n->group_];
	  }
	  uint64_t pos = 0;
	  for(const auto & n : nodes_){
	  	if(groupCounts[n->group_] >= groupSizeCutOff){
	  		nameToNewPos[n->name_] = pos;
	  		++pos;
		  	//std::cout << n->name_ << " : " << n->group_  << " : " << n->value_ << std::endl;
		  	nodes[nCount]["name"] = bib::json::toJson(n->name_);
		  	nodes[nCount]["group"] = bib::json::toJson(n->group_);
		  	++nCount;
	  	}
	  }
	  uint32_t lCount=0;
		for(const auto & e : edges_){
			if(e->on_){
				if(groupCounts[e->nodeToNode_.begin()->second.lock()->group_] >= groupSizeCutOff){
					links[lCount]["source"] = bib::json::toJson(nameToNewPos[nodes_[nameToNodePos_[e->nodeToNode_.begin()->first]]->name_]);
					links[lCount]["target"] = bib::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
					links[lCount]["value"] = bib::json::toJson(e->dist_);
					links[lCount]["on"] = bib::json::toJson(e->on_);
					++lCount;
				}
			}
		}
		return graphJson;
	}


	Json::Value toJsonBestHigest(DIST cutOff,
			uint32_t groupSizeCutOff){
		turnOffEdgesUnder(cutOff);
		resetBestAndVistEdges();
		resetVisitedNodes();
		bool doTies = false;
		allDetermineHigestBest(doTies);
		determineGroups();
		return toJson(groupSizeCutOff);
	}

	Json::Value toJsonBestLowest(DIST cutOff,
			uint32_t groupSizeCutOff){
		turnOffEdgesAbove(cutOff);
		resetBestAndVistEdges();
		resetVisitedNodes();
		bool doTies = false;
		allDetermineLowestBest(doTies);
		determineGroups();
		return toJson(groupSizeCutOff);
	}
}; //class njhUndirWeightedGraph

inline std::vector<bib::color> getColsBetweenExcludeClosest(bib::color first,
																						bib::color last,
																						uint32_t num){
	if(std::abs(first.hue_ - last.hue_) > 180){
		if(first.hue_ < last.hue_){
			first.hue_ += 360;
		}else{
			last.hue_ += 360;
		}
	}
	auto ret = bib::getColsBetweenInc(first.hue_, last.hue_, first.lum_, last.lum_, first.lSat_, last.lSat_, num + 2);
	return std::vector<bib::color>(ret.begin() + 1, ret.end() - 1);
}
template<typename DIST>
class readDistGraph : public njhUndirWeightedGraph<DIST, std::shared_ptr<seqInfo>>{
public:

	/**@b Construct with a distance matrix and a vector of reads that were used to create the distances
	 *
	 * @param distances The distance matrix, the matrix should at least have the diagonal values (each row has as many elements as it's row position
	 * @param reads the reads the distance graph is describing
	 * @todo do a size check for the distance matrix
	 */
	template<typename T>
	readDistGraph(const std::vector<std::vector<DIST>> & distances,
			const std::vector<T> & reads){
	  for(const auto & pos : iter::range(reads.size())){
	  	this->addNode(reads[pos].seqBase_.name_,
	  			std::make_shared<seqInfo>(reads[pos].seqBase_));
	  }
	  for(const auto & pos : iter::range(distances.size())){
	  	for(const auto & subPos : iter::range<uint64_t>(distances[pos].size())){
	  		this->addEdge(reads[pos].seqBase_.name_,
	  				reads[subPos].seqBase_.name_,
	  				distances[pos][subPos]);
	  	}
	  }
	}

	/**@b Construct with a distance matrix and a vector of reads that were used to create the distances, for shared pointer vecs
	 *
	 * @param distances The distance matrix, the matrix should at least have the diagonal values (each row has as many elements as it's row position
	 * @param reads the reads the distance graph is describing
	 * @todo do a size check for the distance matrix
	 */
	template<typename T>
	readDistGraph(const std::vector<std::vector<DIST>> & distances,
			const std::vector<std::shared_ptr<T>> & reads){
	  for(const auto & pos : iter::range(reads.size())){
	  	this->addNode(reads[pos]->seqBase_.name_,
	  			std::make_shared<seqInfo>(reads[pos]->seqBase_));
	  }
	  for(const auto & pos : iter::range(distances.size())){
	  	for(const auto & subPos : iter::range<uint64_t>(distances[pos].size())){
	  		this->addEdge(reads[pos]->seqBase_.name_,
	  				reads[subPos]->seqBase_.name_,
	  				distances[pos][subPos]);
	  	}
	  }
	}


	Json::Value toJsonMismatchGraphAll(bib::color backgroundColor,
			std::unordered_map<std::string, bib::color> nameColors ){
	  Json::Value graphJson;
	  graphJson["backgroundColor"] = "#" +  backgroundColor.hexStr_ ;
	  auto & nodes = graphJson["nodes"];
	  auto & links = graphJson["links"];
	  uint32_t nCount = 0;
	  std::unordered_map<std::string, uint64_t> nameToNewPos;
		uint32_t mismatchesAllowed = 1;
		this->turnOffEdgesAbove(mismatchesAllowed);
		this->determineGroups();
		while(this->numberOfGroups_ > 1){
			++mismatchesAllowed;
			this->resetVisitedNodes();
			this->resetVisitedEdges();
			this->turnOffEdgesAbove(mismatchesAllowed);
			this->determineGroups();
		}

		bib::randomGenerator gen;
	  uint64_t pos = 0;
	  double minReadCnt = std::numeric_limits<double>::max();
	  double maxReadCnt = std::numeric_limits<double>::min();
	  for(const auto & n : this->nodes_){
	  	if(n->value_->cnt_ < minReadCnt){
	  		minReadCnt = n->value_->cnt_;
	  	}
	  	if(n->value_->cnt_ > maxReadCnt){
	  		maxReadCnt = n->value_->cnt_;
	  	}
	  }
	  scale<double> cntScale({minReadCnt, maxReadCnt},{50.0, 1000.0});
	  for(const auto & n : this->nodes_){
  		nameToNewPos[n->name_] = pos;
  		++pos;
	  	//std::cout << n->name_ << " : " << n->group_  << " : " << n->value_ << std::endl;
	  	nodes[nCount]["name"] = bib::json::toJson(n->name_);
	  	nodes[nCount]["group"] = bib::json::toJson(n->group_);
	  	nodes[nCount]["color"] = "#" + nameColors[n->name_].hexStr_;
	  	nodes[nCount]["size"] = cntScale.get(n->value_->cnt_);
	  	++nCount;
	  }
	  uint32_t lCount=0;
		for(const auto & e : this->edges_){
			if(e->on_){
				if(e->dist_ == 0){
					links[lCount]["source"] = bib::json::toJson(nameToNewPos[this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_]);
					links[lCount]["target"] = bib::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
					links[lCount]["value"] = bib::json::toJson(1);
					links[lCount]["on"] = bib::json::toJson(e->on_);
					auto lColor = getColsBetweenExcludeClosest(nameColors[this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_],
							nameColors[e->nodeToNode_.begin()->second.lock()->name_], 1);
					links[lCount]["color"] = "#" + lColor.front().hexStr_;
					++lCount;
				}else{
					std::string lastName = this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_;
					auto lColors = getColsBetweenExcludeClosest(nameColors[this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_],
							nameColors[e->nodeToNode_.begin()->second.lock()->name_], e->dist_ + 1);
					for(const auto & mis : iter::range(e->dist_)){
						std::string newName = this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_
								+ estd::to_string(mis) + e->nodeToNode_.begin()->second.lock()->name_;
			  		nameToNewPos[newName] = pos;
			  		++pos;
				  	nodes[nCount]["name"] = bib::json::toJson(newName);
				  	nodes[nCount]["group"] = bib::json::toJson(e->nodeToNode_.begin()->second.lock()->group_);
				  	nodes[nCount]["color"] = "red";
				  	nodes[nCount]["size"] = 10;
				  	++nCount;
						links[lCount]["source"] = bib::json::toJson(nameToNewPos[lastName]);
						links[lCount]["target"] = bib::json::toJson(nameToNewPos[newName]);
						links[lCount]["value"] = bib::json::toJson(1);
						links[lCount]["on"] = bib::json::toJson(true);
						links[lCount]["color"] = "#" + lColors[mis].hexStr_;
						++lCount;
						lastName = newName;
					}

					links[lCount]["source"] = bib::json::toJson(nameToNewPos[lastName]);
					links[lCount]["target"] = bib::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
					links[lCount]["value"] = bib::json::toJson(1);
					links[lCount]["on"] = bib::json::toJson(true);
					links[lCount]["color"] = "#" + lColors[e->dist_].hexStr_;
					++lCount;
				}
			}
		}
		return graphJson;
	}
	Json::Value toJsonMismatchGraph(uint32_t groupCutOff,
				uint32_t mismatchesAllowed,
				bib::color backgroundColor, double hueStart, double hueStop,
		    double lumStart, double lumStop,
		    double satStart, double satStop){
		Json::Value graphJson;
		graphJson["backgroundColor"] = "#" +  backgroundColor.hexStr_ ;
		auto & nodes = graphJson["nodes"];
		auto & links = graphJson["links"];
		uint32_t nCount = 0;
		std::unordered_map<std::string, uint64_t> nameToNewPos;
		this->turnOffEdgesAbove(mismatchesAllowed);
		this->determineGroups();

		bib::randomGenerator gen;
		uint64_t pos = 0;
		double minReadCnt = std::numeric_limits<double>::max();
		double maxReadCnt = std::numeric_limits<double>::min();
		std::unordered_map<uint32_t, bib::color> groupColors;
	  std::unordered_map<uint32_t, uint32_t> groupCounts;
	  uint32_t numOfCutOffGroups = 0;
	  for(const auto & n : this->nodes_){
	  	++groupCounts[n->group_];
	  }
	  std::vector<uint32_t> groups;
		for(const auto & g : groupCounts){
			if(g.second >= groupCutOff){
				++numOfCutOffGroups;
				groups.emplace_back(g.first);
			}
		}
		//printOutMapContents(groupCounts,"\t", std::cout);
		//printVector(groups);
		auto gColors = bib::getColsBetweenInc(hueStart, hueStop,
	  		lumStart, lumStop,
	  		satStart, satStop,
	  		groups.size());

		for(const auto & pos : iter::range(groups.size())){
			groupColors[groups[pos]] = gColors[pos];
		}

		for(const auto & n : this->nodes_){
			if(n->value_->cnt_ < minReadCnt){
				minReadCnt = n->value_->cnt_;
			}
			if(n->value_->cnt_ > maxReadCnt){
				maxReadCnt = n->value_->cnt_;
			}
		}
		scale<double> cntScale({minReadCnt, maxReadCnt},{50.0, 1000.0});
		for(const auto & n : this->nodes_){
	  	if(groupCounts[n->group_] >= groupCutOff){
				nameToNewPos[n->name_] = pos;
				++pos;
				//std::cout << n->name_ << " : " << n->group_  << " : " << n->value_ << std::endl;
				nodes[nCount]["name"] = bib::json::toJson(n->name_);
				nodes[nCount]["group"] = bib::json::toJson(n->group_);
				nodes[nCount]["color"] = "#" + groupColors[n->group_].hexStr_;
				nodes[nCount]["size"] = cntScale.get(n->value_->cnt_);
				++nCount;
	  	}
		}
		uint32_t lCount=0;
		for(const auto & e : this->edges_){
			if(e->on_ && groupCounts[e->nodeToNode_.begin()->second.lock()->group_] >= groupCutOff){
				if(e->dist_ == 0){
					links[lCount]["source"] = bib::json::toJson(nameToNewPos[this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_]);
					links[lCount]["target"] = bib::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
					links[lCount]["value"] = bib::json::toJson(1);
					links[lCount]["on"] = bib::json::toJson(e->on_);
					auto lColor = getColsBetweenExcludeClosest(groupColors[this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->group_],
							groupColors[e->nodeToNode_.begin()->second.lock()->group_], 1);
					links[lCount]["color"] = "#" + lColor.front().hexStr_;
					++lCount;
				}else{
					std::string lastName = this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_;
					auto lColors = getColsBetweenExcludeClosest(groupColors[this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->group_],
							groupColors[e->nodeToNode_.begin()->second.lock()->group_], e->dist_ + 1);
					for(const auto & mis : iter::range(e->dist_)){
						std::string newName = this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_
								+ estd::to_string(mis) + e->nodeToNode_.begin()->second.lock()->name_;
						nameToNewPos[newName] = pos;
						++pos;
						nodes[nCount]["name"] = bib::json::toJson(newName);
						nodes[nCount]["group"] = bib::json::toJson(e->nodeToNode_.begin()->second.lock()->group_);
						nodes[nCount]["color"] = "red";
						nodes[nCount]["size"] = 10;
						++nCount;
						links[lCount]["source"] = bib::json::toJson(nameToNewPos[lastName]);
						links[lCount]["target"] = bib::json::toJson(nameToNewPos[newName]);
						links[lCount]["value"] = bib::json::toJson(1);
						links[lCount]["on"] = bib::json::toJson(true);
						links[lCount]["color"] = "#" + lColors[mis].hexStr_;
						++lCount;
						lastName = newName;
					}

					links[lCount]["source"] = bib::json::toJson(nameToNewPos[lastName]);
					links[lCount]["target"] = bib::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
					links[lCount]["value"] = bib::json::toJson(1);
					links[lCount]["on"] = bib::json::toJson(true);
					links[lCount]["color"] = "#" + lColors[e->dist_].hexStr_;
					++lCount;
				}
			}
		}
		return graphJson;
	}
};


}  // namespace bibseq
