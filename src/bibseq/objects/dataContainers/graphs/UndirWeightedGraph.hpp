#pragma once
/*
 * UndirWeightedGraph.hpp
 *
 *  Created on: Dec 10, 2014
 *      Author: nickhathaway
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
#include "bibseq/objects/dataContainers/graphs/graphsCommon.hpp"



namespace bibseq {


/**@brief
 *
 * @todo need to go through and make sure checks for if edges and nodes are on for all functions
 */
template<typename DIST, typename VALUE>
class njhUndirWeightedGraph {
public:
	struct dbscanPars{
		DIST eps_;
		uint32_t minEpNeighbors_;
	};

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

		//operator only works if a operator< has been defined for DIST
		bool operator <(const edge & otherEdge){
			return dist_ < otherEdge.dist_;
		}

		//operator only works if a operator> has been defined for DIST
		bool operator >(const edge & otherEdge){
			return dist_ > otherEdge.dist_;
		}

		//checks to see if both nodes are on
		bool nodesOn() const {
			for (const auto & n : nodeToNode_) {
				if (!n.second.lock()->on_) {
					return false;
				}
			}
			return true;
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

		void turnOff(){
			on_ = false;
			for(auto & e : edges_){
				e->on_ = false;
			}
		}

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

		void printEdges(std::ostream & out, bool printVisited) {
			for (const auto & e : edges_) {
				if (e->on_) {
					if (e->visited_ && printVisited) {
						out << name_ << " -- " << e->nodeToNode_[name_].lock()->name_ << " "
								<< e->dist_ << "\n";
					} else {
						out << name_ << " -- " << e->nodeToNode_[name_].lock()->name_ << " "
								<< e->dist_ << "\n";
					}
				}
			}
		}

		DIST getAverageDist() {
			DIST sum = 0;
			double edgeNum = 0;
			for (const auto & e : edges_) {
				if (e->on_) {
					++edgeNum;
					sum += e->dist_;
				}
			}
			return sum / edgeNum;
		}

		DIST getAverageDist(std::function<DIST(std::vector<DIST>)> avgFunc){
			std::vector<DIST> dists;
			for(const auto & e : edges_){
				if(e->on_){
					dists.emplace_back(e->dist_);
				}
			}
			return avgFunc(dists);
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
				for(const auto & e : iter::reversed(edges_)){
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
					if(e->nodeToNode_[name_].lock()->on_ && e->on_ && !e->visited_){
						e->visited_ = true;
						e->nodeToNode_[name_].lock()->visitAndGroupCons(currentGroup);
					}
				}
			}
		}

		std::vector<std::shared_ptr<edge>> getNeighbors() const {
			return edges_;
		}

		std::vector<std::shared_ptr<edge>> getNeighborsEps(double eps) const {
			std::vector<std::shared_ptr<edge>> neighbors;
			for (auto & e : edges_) {
				if (e->dist_ <= eps) {
					neighbors.push_back(e);
				}
			}
			return neighbors;
		}

		void dbscanSpread(uint32_t currentGroup, const dbscanPars & pars) {
			//should never have already been visited, maybe throw if it has?
			visited_ = true;
			++visitedAmount_;
			auto neighbors = getNeighborsEps(pars.eps_);
			if (neighbors.size() >= pars.minEpNeighbors_) {
				for (auto & neigh : neighbors) {
					neigh->on_ = true;
					auto next = neigh->nodeToNode_[name_].lock();
					if(std::numeric_limits<uint32_t>::max() == next->group_){
						next->group_ = currentGroup;
						next->on_ = true;
						if(!next->visited_){
							next->dbscanSpread(currentGroup, pars);
						}
					}
				}
			}
		}

		void dbscanExpand(uint32_t currentGroup, const dbscanPars & pars) {
			//should never have already been visited, maybe throw if it has?
			visited_ = true;
			++visitedAmount_;
			std::vector<std::shared_ptr<edge>> neighbors;
			for (auto & e : edges_) {
				if (e->dist_ <= pars.eps_) {
					neighbors.push_back(e);
				}
			}
			if (neighbors.size() >= pars.minEpNeighbors_) {
				group_ = currentGroup;
				on_ = true;
				for (auto & neigh : neighbors) {
					neigh->on_ = true;
					auto next = neigh->nodeToNode_[name_].lock();
					next->group_ = currentGroup;
					next->on_ = true;
					if (!next->visited_) {
						next->dbscanSpread(currentGroup, pars);
					}
				}
			} else {
				//mark self as noise, leave neighbors alone
				on_ = false;
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

		uint32_t numConnections() const{
			uint32_t ret = 0;
			for(auto & e : edges_){
				if(e->on_){
					++ret;
				}
			}
			return ret;
		}

		void reset() {
			on_ = true;
			visited_ = false;
			visitedAmount_ = 0;
			group_ = std::numeric_limits<uint32_t>::max();
		}

	}; // class node



	std::vector<std::shared_ptr<node>> nodes_;
	std::vector<std::shared_ptr<edge>> edges_;
	std::unordered_map<std::string, uint64_t> nameToNodePos_;
	uint32_t numberOfGroups_ = 1;

	void addNode(const std::string & uid, const VALUE & value){
		if(bib::in(uid, nameToNodePos_)){
			throw std::runtime_error{std::string(__PRETTY_FUNCTION__) + ": already contains node with uid: " + uid + ", can't have duplicate names"};
		}
		nameToNodePos_[uid] = nodes_.size();
		nodes_.emplace_back(std::make_shared<node>(uid, value));
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

	void resetAllNodes(){
		for(const auto & n : nodes_){
			n->reset();
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

	template<typename COMP>
	void turnOnEdgesWtihComp(const DIST & cutOff, COMP comp) {
		for (const auto & c : edges_) {
			if(c->nodesOn()){
				if (comp(c->dist_, cutOff)) {
					c->on_ = true;
				} else {
					c->on_ = false;
				}
			}
		}
	}
	template<typename COMP>
	void turnOffEdgesWithComp(const DIST & cutOff, COMP comp) {
		for (const auto & c : edges_) {
			if(c->nodesOn()){
				if (comp(c->dist_, cutOff)) {
					c->on_ = false;
				} else {
					c->on_ = true;
				}
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

	void printAdjByGroup(std::ostream & out){
		std::unordered_map<uint32_t, VecStr> groupOut;
		for(const auto & e : edges_){
			if(e->on_){
				std::stringstream ss;
				ss << e->nodeToNode_.begin()->second.lock()->group_ << " "
						<< e->nodeToNode_.begin()->first << " -- "
						<< e->nodeToNode_.begin()->second.lock()->name_ << " "
						<< e->dist_;
				groupOut[e->nodeToNode_.begin()->second.lock()->group_].emplace_back(ss.str());
			}
		}
		for(const auto & g : groupOut){
			out << "Group: " << g.first << std::endl;
			for(const auto & go : g.second){
				out << "\t" << go << std::endl;
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
			if(n->on_ && !n->visited_){
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

	void dbscan(const dbscanPars & pars) {
		//reset node's visited and group values
		this->resetAllNodes();
		//turn off whole graph
		this->turnOffAllCons();
		this->turnOffAllNodes();
		//set number of groups to be 0
		numberOfGroups_ = 0;
		for (auto & n : nodes_) {
			//if the node has not be visited by an expand or spread try to expand it
			if (!n->visited_) {
				n->dbscanExpand(numberOfGroups_, pars);
				//if it was assigned a group and expanded, increase group number
				if (std::numeric_limits<uint32_t>::max() != n->group_) {
					++numberOfGroups_;
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
	  uint32_t groupsAboveCutOff =0;
	  std::set<uint32_t> largeGroups;
	  for(const auto & n : nodes_){
	  	if(groupCounts[n->group_] >= groupSizeCutOff){
	  		largeGroups.emplace(n->group_);
	  	}
	  }
	  groupsAboveCutOff = largeGroups.size();
	  std::vector<bib::color> groupColors;
	  if(groupsAboveCutOff > 1){
	  	groupColors = bib::getColsBetweenInc(0, 360 - 360.0/groupsAboveCutOff, 0.5, 0.5, 1,1, groupsAboveCutOff);
	  }else{
	  	groupColors = {bib::color{"#FF0000"}};
	  }
	  std::unordered_map<uint32_t, bib::color> gColors;
	  for(auto e : iter::enumerate(largeGroups)){
	  	gColors[e.element] = groupColors[e.index];
	  }
	  uint64_t pos = 0;
	  for(const auto & n : nodes_){
	  	if(groupCounts[n->group_] >= groupSizeCutOff){
	  		nameToNewPos[n->name_] = pos;
	  		++pos;
		  	//std::cout << n->name_ << " : " << n->group_  << " : " << n->value_ << std::endl;
		  	nodes[nCount]["name"] = bib::json::toJson(n->name_);
		  	nodes[nCount]["group"] = bib::json::toJson(n->group_);
		  	nodes[nCount]["color"] = bib::json::toJson("#" + gColors[n->group_].hexStr_);
		  	nodes[nCount]["size"]  = 30;
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
					links[lCount]["color"] = bib::json::toJson("#" + gColors[e->nodeToNode_.begin()->second.lock()->group_].hexStr_);
					++lCount;
				}
			}
		}
		return graphJson;
	}

	Json::Value toJson(uint32_t groupSizeCutOff, std::unordered_map<std::string, std::string> nameToColor){
	  Json::Value graphJson;
	  auto & nodes = graphJson["nodes"];
	  auto & links = graphJson["links"];
	  uint32_t nCount = 0;
	  std::unordered_map<uint32_t, uint32_t> groupCounts;
	  std::unordered_map<std::string, uint64_t> nameToNewPos;

	  for(const auto & n : nodes_){
	  	++groupCounts[n->group_];
	  }
	  uint32_t groupsAboveCutOff =0;
	  std::set<uint32_t> largeGroups;
	  for(const auto & n : nodes_){
	  	if(groupCounts[n->group_] >= groupSizeCutOff){
	  		largeGroups.emplace(n->group_);
	  	}
	  }
	  groupsAboveCutOff = largeGroups.size();
	  std::vector<bib::color> groupColors;

	  if(nameToColor.empty()){
	  	groupColors = bib::getColsBetweenInc(0, 360 - 360.0/groupsAboveCutOff, 0.5, 0.5, 1,1, groupsAboveCutOff);
		  std::unordered_map<uint32_t, bib::color> gColors;
		  for(auto e : iter::enumerate(largeGroups)){
		  	gColors[e.element] = groupColors[e.index];
		  }
	  }
	  uint64_t pos = 0;
	  for(const auto & n : nodes_){
	  	if(groupCounts[n->group_] >= groupSizeCutOff){
	  		nameToNewPos[n->name_] = pos;
	  		++pos;
		  	//std::cout << n->name_ << " : " << n->group_  << " : " << n->value_ << std::endl;
		  	nodes[nCount]["name"] = bib::json::toJson(n->name_);
		  	nodes[nCount]["group"] = bib::json::toJson(n->group_);
		  	if(nameToColor.empty()){
		  		nodes[nCount]["color"] = bib::json::toJson(groupColors[n->group_].getHexStr());
		  	}else{
		  		nodes[nCount]["color"] = bib::json::toJson(nameToColor[n->name_]);
		  	}
		  	nodes[nCount]["size"]  = 30 ;
		  	++nCount;
	  	}
	  }
	  uint32_t lCount = 0;
		for(const auto & e : edges_){
			if(e->on_){
				if(groupCounts[e->nodeToNode_.begin()->second.lock()->group_] >= groupSizeCutOff){
					links[lCount]["source"] = bib::json::toJson(nameToNewPos[nodes_[nameToNodePos_[e->nodeToNode_.begin()->first]]->name_]);
					links[lCount]["target"] = bib::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
					links[lCount]["value"] = bib::json::toJson(e->dist_);
					links[lCount]["on"] = bib::json::toJson(e->on_);
			  	if(nameToColor.empty()){
			  		links[lCount]["color"] = bib::json::toJson(groupColors[e->nodeToNode_.begin()->second.lock()->group_].getHexStr());
			  	}else{
			  		links[lCount]["color"] = bib::json::toJson(nameToColor[e->nodeToNode_.begin()->second.lock()->name_]);
			  	}
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








}  // namespace bibseq
