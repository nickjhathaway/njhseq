#pragma once
/*
 * UndirWeightedGraph.hpp
 *
 *  Created on: Dec 10, 2014
 *      Author: nickhathaway
 */
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
#include <boost/math/distributions/normal.hpp>
#undef BOOST_HAS_THREADS
#include <boost/math/statistics/t_test.hpp>
#include "njhseq/objects/dataContainers/graphs/graphsCommon.hpp"
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/concurrency/PairwisePairFactory.hpp"
#include "njhseq/concurrency/AllByAllPairFactory.hpp"


namespace njhseq {
/**@brief
 *
 * @todo need to go through and make sure checks for if edges and nodes are on for all functions
 */
template<typename DIST, typename VALUE>
class njhUndirWeightedGraph {
public:
  struct dbscanPars {
    dbscanPars(const DIST &eps, const uint32_t minEpNeighbors) : eps_(eps), minEpNeighbors_(minEpNeighbors) {

    }

    dbscanPars() = default;

    DIST eps_;
    uint32_t minEpNeighbors_{};
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
				on_(true), visited_(false), corePoint_(false){

		}
		std::string name_;
		VALUE value_;
		std::vector<std::shared_ptr<edge>> edges_;
		bool on_;
		bool visited_;
		bool corePoint_;
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
			//only count neighbors from current group or unassigned neighbors
			std::vector<std::shared_ptr<edge>> filtNeighbors;
			for(const auto & neigh : neighbors){
				auto next = neigh->nodeToNode_[name_].lock();
				if(std::numeric_limits<uint32_t>::max() == next->group_ || currentGroup == next->group_){
					filtNeighbors.emplace_back(neigh);
				}
			}
			//based off the r implementation i think the origin point counts in the neighbor counts
			if (filtNeighbors.size() + 1 >= pars.minEpNeighbors_) {
				corePoint_ = true;
				for (auto & neigh : filtNeighbors) {
					auto next = neigh->nodeToNode_[name_].lock();
					if(std::numeric_limits<uint32_t>::max() == next->group_){
						neigh->on_ = true;
						next->group_ = currentGroup;
						next->on_ = true;
						if(!next->visited_){
							next->dbscanSpread(currentGroup, pars);
						}
					}else if(group_ == next->group_){
						//can turn on edges that connect same gorup
						neigh->on_ = true;
					}
				}
			}
		}

		void dbscanExpand(uint32_t currentGroup, const dbscanPars & pars) {
			//should never have already been visited, maybe throw if it has?
			visited_ = true;
			++visitedAmount_;
			auto neighbors = getNeighborsEps(pars.eps_);
			std::vector<std::shared_ptr<edge>> filtNeighbors;
			for(const auto & neigh : neighbors){
				auto next = neigh->nodeToNode_[name_].lock();
				if(std::numeric_limits<uint32_t>::max() == next->group_ || currentGroup == next->group_){
					filtNeighbors.emplace_back(neigh);
				}
			}
			//based off the r implementation i think the origin point counts in the neighbor counts
			if (filtNeighbors.size() + 1 >= pars.minEpNeighbors_) {
				group_ = currentGroup;
				on_ = true;
				corePoint_ = true;
				for (auto & neigh : filtNeighbors) {
					//if neighbor is unassigned add
					auto next = neigh->nodeToNode_[name_].lock();
					if(std::numeric_limits<uint32_t>::max() == next->group_){
						neigh->on_ = true;

						next->group_ = currentGroup;
						next->on_ = true;
						if (!next->visited_) {
							next->dbscanSpread(currentGroup, pars);
						}
					}else if(group_ == next->group_){
						//can turn on edges that connect same gorup
						neigh->on_ = true;
					}
				}
			} else {
				//mark self as noise, leave neighbors alone
				on_ = false;
			}
		}

		void determineHighestDistance(bool doTies){
			if(!edges_.empty() && numConnections()> 0){
				DIST best = edges_.front()->dist_;
				std::vector<uint64_t> bestIndex = {0};
				for(const auto ePos : iter::range(edges_.size())){
					const auto & e = edges_[ePos];
					if(e->on_){
						best = e->dist_;
						bestIndex.clear();
						bestIndex.emplace_back(ePos);
						break;
					}
				}
				for(const auto e : iter::enumerate(edges_)){
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
				if(doTies){
					for(const auto & b : bestIndex){
						edges_[b]->best_ = true;
					}
				}else{
					edges_[bestIndex.front()]->best_ = true;
				}
			}
		}

		void determineLowestDistance(bool doTies){
			if(!edges_.empty() && numConnections()> 0){
				DIST best = edges_.front()->dist_;
				std::vector<uint64_t> bestIndex = {0};
				for(const auto ePos : iter::range(edges_.size())){
					const auto & e = edges_[ePos];
					if(e->on_){
						best = e->dist_;
						bestIndex.clear();
						bestIndex.emplace_back(ePos);
						break;
					}
				}
				for(const auto e : iter::enumerate(edges_)){
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
				if(doTies){
					for(const auto & b : bestIndex){
						edges_[b]->best_ = true;
					}
				}else{
					edges_[bestIndex.front()]->best_ = true;
				}
			}
		}

		void determineBestDistanceWithComp(bool doTies,
				std::function<bool(const DIST&,const DIST&)> distCompFunc,
				std::function<bool(const DIST&,const DIST&)> equalCompFunc){
			if(!edges_.empty() && numConnections()> 0){
				DIST best = edges_.front()->dist_;
				std::vector<uint64_t> bestIndex = {0};
				for(const auto ePos : iter::range(edges_.size())){
					const auto & e = edges_[ePos];
					if(e->on_){
						best = e->dist_;
						bestIndex.clear();
						bestIndex.emplace_back(ePos);
						break;
					}
				}
				for(const auto e : iter::enumerate(edges_)){
					if(e.element->on_){
						if(distCompFunc(e.element->dist_,best)){
							best = e.element->dist_;
							bestIndex.clear();
							bestIndex = {e.index};
						} else if (doTies && equalCompFunc(e.element->dist_, best)){
							bestIndex.emplace_back(e.index);
						}
					}
				}
				if(doTies){
					for(const auto & b : bestIndex){
						edges_[b]->best_ = true;
					}
				}else{
					edges_[bestIndex.front()]->best_ = true;
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

		void removeOffEdges(){
			edges_.erase(std::remove_if(edges_.begin(), edges_.end(), [](const std::shared_ptr<edge> & e){
				return !e->on_;
			}),edges_.end());
		}

	}; // class node



	std::vector<std::shared_ptr<node>> nodes_;
	std::vector<std::shared_ptr<edge>> edges_;
	std::unordered_map<std::string, uint64_t> nameToNodePos_;
	uint32_t numberOfGroups_ = 1;

	void resetNodePositionMap(){
		nameToNodePos_.clear();
		for(const auto pos : iter::range(nodes_.size())){
			nameToNodePos_[nodes_[pos]->name_] = pos;
		}
	}

	void removeOffNodes(){
		std::vector<uint32_t> toRemove;
		for(const auto nodePos : iter::range(nodes_.size())){
			const auto & n = nodes_[nodePos];
			if(!n->on_){
				for(const auto & edge : n->edges_){
					//also turn off edges so they can be removed as well
					edge->on_ = false;
				}
				toRemove.emplace_back(nodePos);
			}
		}
		if(!toRemove.empty()){
			std::sort(toRemove.rbegin(), toRemove.rend());
			for(const auto & remove : toRemove){
				nodes_.erase(nodes_.begin() + remove);
			}
			removeOffEdges();
			resetNodePositionMap();
		}
	}

	void removeOffEdges(){
		for(const auto & n : nodes_){
			n->removeOffEdges();
		}

		edges_.erase(std::remove_if(edges_.begin(), edges_.end(), [](const std::shared_ptr<edge> & e){
			return !e->on_;
		}), edges_.end());
	}

	void addNode(const std::string & uid, const VALUE & value){
		if(njh::in(uid, nameToNodePos_)){
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

	void allDetermineBestDistanceWithComp(bool doTies, std::function<bool(const DIST&,const DIST&)> distCompFunc,
			std::function<bool(const DIST&,const DIST&)> equalCompFunc){
		for(const auto & n : nodes_){
			n->determineBestDistanceWithComp(doTies, distCompFunc, equalCompFunc);
		}
		for(const auto & e : edges_){
			if(e->on_){
				if(!e->best_){
					e->on_ = false;
				}
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
		//std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;

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
				//std::cout << n->value_->name_ << std::endl;
				n->dbscanExpand(numberOfGroups_, pars);
				//if it was assigned a group and expanded, increase group number
				if (std::numeric_limits<uint32_t>::max() != n->group_) {
					++numberOfGroups_;
				}
			}
		}
	}

	void assignNoiseNodesAGroup() {
		for (auto & n : nodes_) {
			if (std::numeric_limits<uint32_t>::max() == n->group_) {
				n->group_ = numberOfGroups_;
				++numberOfGroups_;
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


	Json::Value toJson(uint32_t groupSizeCutOff,
			std::unordered_map<std::string, std::string> nameToColor = {}){
	  Json::Value graphJson;
	  auto & nodes = graphJson["nodes"];
	  auto & links = graphJson["links"];
	  uint32_t nCount = 0;
	  std::unordered_map<uint32_t, uint32_t> groupCounts;
	  std::unordered_map<std::string, uint64_t> nameToNewPos;

	  for(const auto & n : nodes_){
	  	++groupCounts[n->group_];
	  }
	  std::set<uint32_t> largeGroups;
	  for(const auto & n : nodes_){
	  	if(groupCounts[n->group_] >= groupSizeCutOff){
	  		largeGroups.emplace(n->group_);
	  	}
	  }
	  std::unordered_map<uint32_t, njh::color> gColors;
	  uint32_t groupsAboveCutOff = largeGroups.size();
	  std::vector<njh::color> groupColors;
	  if(groupsAboveCutOff > 1){
	  	groupColors = njh::getColsBetweenInc(0, 360 - 360.0/groupsAboveCutOff, 0.5, 0.5, 1,1, groupsAboveCutOff);
	  }else{
	  	groupColors = {njh::color{"#FF0000"}};
	  }
	  for(const auto e : iter::enumerate(largeGroups)){
	  	gColors[e.element] = groupColors[e.index];
	  }
	  uint64_t pos = 0;
	  for(const auto & n : nodes_){
	  	if(groupCounts[n->group_] >= groupSizeCutOff){
	  		nameToNewPos[n->name_] = pos;
	  		++pos;
		  	//std::cout << n->name_ << " : " << n->group_  << " : " << n->value_ << std::endl;
		  	nodes[nCount]["name"] = njh::json::toJson(n->name_);
		  	nodes[nCount]["group"] = njh::json::toJson(n->group_);
		  	nodes[nCount]["corePoint"] = njh::json::toJson(n->corePoint_);
		  	if(nameToColor.empty() || !njh::in(n->name_, nameToColor)){
		  		nodes[nCount]["color"] = njh::json::toJson(gColors[n->group_].getHexStr());
		  	}else{
		  		nodes[nCount]["color"] = njh::json::toJson(nameToColor[n->name_]);
		  	}
		  	nodes[nCount]["size"]  = 30 ;
		  	++nCount;
	  	}
	  }
	  uint32_t lCount = 0;
		for(const auto & e : edges_){
			if(e->on_){
				if(groupCounts[e->nodeToNode_.begin()->second.lock()->group_] >= groupSizeCutOff){
					links[lCount]["source"] = njh::json::toJson(nameToNewPos[nodes_[nameToNodePos_[e->nodeToNode_.begin()->first]]->name_]);
					links[lCount]["target"] = njh::json::toJson(nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_]);
					links[lCount]["value"] = njh::json::toJson(e->dist_);
					links[lCount]["on"] = njh::json::toJson(e->on_);
			  	if(nameToColor.empty() || !njh::in(e->nodeToNode_.begin()->second.lock()->name_, nameToColor)){
			  		links[lCount]["color"] = njh::json::toJson(gColors[e->nodeToNode_.begin()->second.lock()->group_].getHexStr());
			  	}else{
			  		links[lCount]["color"] = njh::json::toJson(nameToColor[e->nodeToNode_.begin()->second.lock()->name_]);
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


  struct HDBScanInputPars {
    uint32_t proposedClusters{std::numeric_limits<uint32_t>::max()};
    bool debug{false};
    bool verbose{false};
    bool HDBScountZeroNeighbors{true};
    uint32_t numThreads{1};
    bool HDBSCountSingletGroups{false};
    bool HDBSredetermineMaxEps{false};
    double HDBSmaxInitialEps{std::numeric_limits<double>::max()};
    uint32_t groupDiffToReCalc = 20;
    //dbscanPars initialDbPars{2, 0.20};
  };

  struct HDBScanResults{
    table hdbsRunInfo{VecStr{"Run", "nonSingletGroupCounts","totalGroupCount","lowestCentroidDist"}};
    table nearestNeibhorDists{VecStr{"name", "neighbor", "nonZeroNeighborPos", "dist"}};
    table tTests{VecStr{"group", "groupSize", "t-statistic", "p-value", "diffMean", "diffSD", "groupMean", "groupSD"}};

    table initialGroupNames{VecStr{"name", "group"}};
    table initial_centroidDistances;
    table final_centroidDistances;
  };

  HDBScanResults runHDBScan(const HDBScanInputPars & pars, const std::vector<std::vector<DIST>> & dist){
    HDBScanResults ret;

    uint32_t minimum_minPts = 3;
    uint32_t proposedClusters = std::numeric_limits<uint32_t>::max() == pars.proposedClusters ?  nodes_.size()/minimum_minPts : pars.proposedClusters;

    if(std::numeric_limits<uint32_t>::max() == proposedClusters){
      proposedClusters = nodes_.size()/minimum_minPts;
    }
    if(pars.verbose){
      std::cout << "using: " << proposedClusters << " proposed clusters" << std::endl;
    }
    //HDBS

    uint32_t maximum_minPts = nodes_.size()/proposedClusters;
    if(0 == maximum_minPts){
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "maximum_minPts can't be 0" << "\n";
      throw std::runtime_error { ss.str() };
    }
    if(pars.verbose){
      std::cout << "using " << maximum_minPts << " for maximum min pts" << std::endl;
    }
    DIST minEps = std::numeric_limits<DIST>::max();
    std::vector<DIST> neighbor2Dists;
//    std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
    for(const auto & n : nodes_){
      //sort so edges are sorted from lowest dist to highest dist
      n->sortEdges([](const std::shared_ptr<edge> & e1,
                      const std::shared_ptr<edge> & e2){
        return *e1 < *e2;
      });

      uint32_t total = 0;
      uint32_t count = 0;
      for(const auto & e : n->edges_){
        ++total;
        if(pars.HDBScountZeroNeighbors || e->dist_ != 0){
          ++count;
          ret.nearestNeibhorDists.addRow(n->value_->name_, total, count, e->dist_);
          if((minimum_minPts - 1) == count){
            neighbor2Dists.emplace_back(e->dist_);
            if(e->dist_ < minEps){
              minEps = e->dist_;
            }
            break;
          }
        }
      }
    }
//    std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
    njh::sort(neighbor2Dists);
    std::vector<DIST> eps;
    auto sq = static_cast<uint32_t>(std::round(std::sqrt(nodes_.size())));
    for(uint32_t i = sq; i < nodes_.size(); i += sq){
      eps.emplace_back(neighbor2Dists[i]);
    }
//    std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
    double maxEps = std::numeric_limits<DIST>::lowest();
    for(const auto & n : nodes_){
      //since all nodes are connected no need to check size edge
//      std::cout << "n->edges_.size(): " << n->edges_.size() << std::endl;
      if(n->edges_.size() >= maximum_minPts && n->edges_[maximum_minPts -1]->dist_ > maxEps){
        maxEps = n->edges_[maximum_minPts - 1]->dist_;
      }
    }
//    exit(1);
//    std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
    //std::cout << __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ << std::endl;
    if(pars.debug){
      std::cout << "Eps: " << std::endl;
      for(const auto & epEnum : iter::enumerate(eps)){
        std::cout << "\t"<< epEnum.index << ": " << epEnum.second << std::endl;
      }
      std::cout << "Max Eps: " << maxEps << std::endl;
      std::cout << "maximum_minPts: " << maximum_minPts << std::endl;
      std::cout << "medianNearestNeighbor: " << vectorMedianRef(neighbor2Dists) << std::endl;
      std::cout << "meanNearestNeighbor: " << vectorMean(neighbor2Dists) << std::endl;
      std::cout << "  sdNearestNeighbor: " << vectorStandardDeviationPop(neighbor2Dists) << std::endl;
    }
    //reset node's visited and group values
    resetAllNodes();
    //turn off whole graph
    turnOffAllCons();
    turnOffAllNodes();
    //set number of groups to be 0
    numberOfGroups_ = 0;
    for(const auto & epEnum : iter::enumerate(eps)){
      if(epEnum.second > pars.HDBSmaxInitialEps){
        break;
      }
      dbscanPars currentPars;
      currentPars.eps_ = epEnum.second;
      currentPars.minEpNeighbors_ = minimum_minPts;
      for (auto & n : nodes_) {
        //if the node has not be visited by an expand or spread try to expand it
        if (!n->visited_) {
          //std::cout << n->value_->name_ << std::endl;
          n->dbscanExpand(numberOfGroups_, currentPars);
          //if it was assigned a group and expanded, increase group number
          if (std::numeric_limits<uint32_t>::max() != n->group_) {
            ++numberOfGroups_;
          }
        }
      }
      if(epEnum.index + 1 != eps.size() && eps[epEnum.first + 1] < pars.HDBSmaxInitialEps ) {
        //reset unclustered nodes back on and unvisted for the next eps
        for(auto & n : nodes_){
          if(!n->on_){
            n->visitedAmount_ = 0;
            n->visited_ = false;
            n->on_ = true;
          }
        }
      }
    }
    assignNoiseNodesAGroup();
    if(pars.HDBSredetermineMaxEps){
      std::vector<DIST> notSameGroupDists;
      for(const auto & n : nodes_){
        for(const auto & e : n->edges_){
          if(e->nodeToNode_[n->name_].lock()->group_ != n->group_){
            notSameGroupDists.emplace_back(e->dist_);
          }
        }
      }
      auto newHDBSmaxInitialEps = vectorMean(notSameGroupDists) - 2 * vectorStandardDeviationPop(notSameGroupDists);
      if(newHDBSmaxInitialEps < 0){
        newHDBSmaxInitialEps = vectorMean(notSameGroupDists);
      }
      if(pars.debug){
        std::cout << "medianNotSameGroupDists: " << vectorMedianRef(notSameGroupDists) << std::endl;
        std::cout << "meanNotSameGroupDists: " << vectorMean(notSameGroupDists) << std::endl;
        std::cout << "  sdNotSameGroupDists: " << vectorStandardDeviationPop(notSameGroupDists) << std::endl;
        std::cout << std::endl;
        std::cout << "max eps based on not same group =" << newHDBSmaxInitialEps << std::endl;
      }
      //reset node's visited and group values
      resetAllNodes();
      //turn off whole graph
      turnOffAllCons();
      turnOffAllNodes();
      //set number of groups to be 0
      numberOfGroups_ = 0;
      for(const auto & epEnum : iter::enumerate(eps)){
        if(epEnum.second > newHDBSmaxInitialEps){
          break;
        }
        dbscanPars currentPars;
        currentPars.eps_ = epEnum.second;
        currentPars.minEpNeighbors_ = minimum_minPts;
        for (auto & n : nodes_) {
          //if the node has not be visited by an expand or spread try to expand it
          if (!n->visited_) {
            //std::cout << n->value_->name_ << std::endl;
            n->dbscanExpand(numberOfGroups_, currentPars);
            //if it was assigned a group and expanded, increase group number
            if (std::numeric_limits<uint32_t>::max() != n->group_) {
              ++numberOfGroups_;
            }
          }
        }
        if(epEnum.index + 1 != eps.size() && eps[epEnum.first + 1] < newHDBSmaxInitialEps ) {
          //reset unclustered nodes back on and unvisted for the next eps
          for(auto & n : nodes_){
            if(!n->on_){
              n->visitedAmount_ = 0;
              n->visited_ = false;
              n->on_ = true;
            }
          }
        }
      }
      assignNoiseNodesAGroup();
    }
    std::vector<double> differentInitialGroupDists;

    double lowestCentroidDistInitial = std::numeric_limits<DIST>::max();
    std::vector<std::vector<DIST>> centroidDistances;
    for(const auto pos : iter::range(numberOfGroups_)) {
      centroidDistances.emplace_back(std::vector<DIST>(pos + 1));
    }
    //set the zeros
    for(const auto pos : iter::range(numberOfGroups_)) {
      centroidDistances[pos][pos] = 0;
    }
    {
      //initial centroid info
      //compute centroid distances
      std::map<uint32_t, std::vector<uint32_t>> groupNodes;
      for(const auto & node : nodes_){
        groupNodes[node->group_].emplace_back(nameToNodePos_[node->name_]);
      }
      auto groups = getVectorOfMapKeys(groupNodes);
      PairwisePairFactory pFac(groups.size());
      PairwisePairFactory::PairwisePair pair;
      while(pFac.setNextPair(pair)){
        //group1 == row_
        //group2 == col_
        uint32_t group1 = groups[pair.row_];
        uint32_t group2 = groups[pair.col_];
        uint32_t group1Size = groupNodes[group1].size();
        uint32_t group2Size = groupNodes[group2].size();
        double sumOfSquaresAll = 0;
        double sumOfSquaresGroup1 = 0;
        double sumOfSquaresGroup2 = 0;
        //group 1
        if(group1Size > 1){
          PairwisePairFactory group1_pFac(group1Size);
          PairwisePairFactory::PairwisePair group1_pair;
          while(group1_pFac.setNextPair(group1_pair)){

            uint32_t node1 = groupNodes[group1][group1_pair.col_];
            uint32_t node2 = groupNodes[group1][group1_pair.row_];
            uint32_t distRow = std::max(node1,  node2);
            uint32_t distCol = std::min(node1,  node2);
            double squareDist = std::pow(dist[distRow][distCol], 2.0);
            sumOfSquaresGroup1 += squareDist;
            sumOfSquaresAll += squareDist;
          }
        }
        //group 2
        if(group2Size > 1){
          PairwisePairFactory group2_pFac(group2Size);
          PairwisePairFactory::PairwisePair group2_pair;
          while(group2_pFac.setNextPair(group2_pair)){
            uint32_t node1 = groupNodes[group2][group2_pair.col_];
            uint32_t node2 = groupNodes[group2][group2_pair.row_];
            uint32_t distRow = std::max(node1,  node2);
            uint32_t distCol = std::min(node1,  node2);
            double squareDist = std::pow(dist[distRow][distCol], 2.0);
            sumOfSquaresGroup2 += squareDist;
            sumOfSquaresAll += squareDist;
          }
        }
        //between group
        for(const auto & group1Node : groupNodes[group1]){
          for(const auto & group2Node : groupNodes[group2]){
            uint32_t distRow = std::max(group1Node,  group2Node);
            uint32_t distCol = std::min(group1Node,  group2Node);
            double squareDist = std::pow(dist[distRow][distCol], 2.0);
            differentInitialGroupDists.emplace_back(dist[distRow][distCol]);
            sumOfSquaresAll += squareDist;
          }
        }
        double distBetweenCentroids = (sumOfSquaresAll - (group1Size + group2Size) * (sumOfSquaresGroup1 / group1Size + sumOfSquaresGroup2 / group2Size)) / (group1Size * group2Size);
        centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)] = distBetweenCentroids;
        if(distBetweenCentroids < lowestCentroidDistInitial) {
          lowestCentroidDistInitial = distBetweenCentroids;
        }
      }
      if(pars.debug){
        double differentInitialGroupDists_mean = vectorMean(differentInitialGroupDists);
        double differentInitialGroupDists_sd = vectorStandardDeviationPop(differentInitialGroupDists);
        std::cout << "differentInitialGroupDists_mean: " << differentInitialGroupDists_mean << std::endl;
        std::cout << "differentInitialGroupDists_sd: " << differentInitialGroupDists_sd << std::endl;

        boost::math::normal_distribution diffDistr(differentInitialGroupDists_mean, differentInitialGroupDists_sd);

        for(const auto & group : groups){
          uint32_t groupSize = groupNodes[group].size();
          if(groupSize > 1){
            std::vector<double> sameGroupDist;
            PairwisePairFactory group_pFac(groupSize);
            PairwisePairFactory::PairwisePair group_pair;
            while(group_pFac.setNextPair(group_pair)){

              uint32_t node1 = groupNodes[group][group_pair.col_];
              uint32_t node2 = groupNodes[group][group_pair.row_];
              uint32_t distRow = std::max(node1,  node2);
              uint32_t distCol = std::min(node1,  node2);
              sameGroupDist.emplace_back(dist[distRow][distCol]);

            }

            auto [t, p] = boost::math::statistics::two_sample_t_test(differentInitialGroupDists, sameGroupDist);

            ret.tTests.addRow(
                group, groupSize, t, p, differentInitialGroupDists_mean, differentInitialGroupDists_sd,
                vectorMean(sameGroupDist), vectorStandardDeviationSamp(sameGroupDist)
            );
          } else {
            ret.tTests.addRow(
                group, groupSize, "NA", "NA", differentInitialGroupDists_mean, differentInitialGroupDists_sd, "NA", "NA"
            );
          }
        }
      }
      {
        //initial infos
        for(const auto & n : nodes_){
          ret.initialGroupNames.addRow(n->name_, n->group_);
        }
        //write out distance matrix
        ret.initial_centroidDistances = table(toVecStr(njh::getVecOfMapKeys(groupNodes)));
        for (const auto & rowIdx : iter::range(centroidDistances.size())) {
          auto outVec = centroidDistances[rowIdx];
          for(const auto & addRowIdx : iter::range(rowIdx + 1, centroidDistances.size())){
            outVec.emplace_back(centroidDistances[addRowIdx][rowIdx]);
          }
          ret.initial_centroidDistances.addRow(toVecStr(outVec));
        }
      }
    }
    uint32_t numberOfNonSingletClusters = 0;
    {
      auto groupCounts = getGroupCounts();
      for(const auto & count : groupCounts){
        if(pars.HDBSCountSingletGroups || count.second > 1){
          ++numberOfNonSingletClusters;
        }
      }
    }
    ret.hdbsRunInfo.addRow(
        "0", numberOfNonSingletClusters, numberOfGroups_, lowestCentroidDistInitial
    );
    uint32_t runCount = 0;
    if(pars.verbose){
      std::cout << "Initial Number of Clusters: " << numberOfNonSingletClusters << std::endl;
      std::cout << "Proposed Clusters: " << proposedClusters << std::endl;
    }
    auto groupCountsSinceLastReCalc = getGroupCounts();
    while(numberOfNonSingletClusters > proposedClusters){
      ++runCount;
      if(pars.verbose){
        std::cout << "Run: " << runCount << std::endl;
        std::cout << "Current Number of Clusters: " << numberOfNonSingletClusters << std::endl;
        std::cout << "Proposed Clusters: " << proposedClusters << std::endl;
      }
      double lowestCentroidDist = std::numeric_limits<double>::max();
      std::set<uint32_t> groupsToCollapse;
      //get lowest centroid distances
      std::map<uint32_t, std::vector<uint32_t>> groupNodes;
      for(const auto & node : nodes_){
        groupNodes[node->group_].emplace_back(nameToNodePos_[node->name_]);
      }
      auto groups = getVectorOfMapKeys(groupNodes);
      PairwisePairFactory pFac(groups.size());
      PairwisePairFactory::PairwisePair pair;
      while(pFac.setNextPair(pair)){
        uint32_t group1 = groups[pair.row_];
        uint32_t group2 = groups[pair.col_];
        double distBetweenCentroids = centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)];
        if(distBetweenCentroids < lowestCentroidDist){
          lowestCentroidDist = distBetweenCentroids;
          groupsToCollapse.clear();
          groupsToCollapse.emplace(groups[pair.row_]);
          groupsToCollapse.emplace(groups[pair.col_]);
        }else if(distBetweenCentroids == lowestCentroidDist){
          groupsToCollapse.emplace(groups[pair.row_]);
          groupsToCollapse.emplace(groups[pair.col_]);
        }
      }
      if(pars.verbose){
        std::cout << "lowestCentroidDist: " << lowestCentroidDist << std::endl;
      }
      //collapse group by setting group value
      for(const auto & group : groupsToCollapse){
        for(const auto & nodeIdx : groupNodes[group]){
          nodes_[nodeIdx]->group_ = *groupsToCollapse.begin();
        }
        if(*groupsToCollapse.begin() != group){
          //add this group's nodes to the collapsed to group
          njh::addOtherVec(groupNodes[*groupsToCollapse.begin()], groupNodes[group]);
        }
      }
      numberOfNonSingletClusters = 0;
      auto currentGroupCounts = getGroupCounts();
      numberOfGroups_ = currentGroupCounts.size();

      for(const auto & count : currentGroupCounts){
        if(pars.HDBSCountSingletGroups || count.second > 1){
          ++numberOfNonSingletClusters;
        }
      }
      //re-compute centroid distances now that clusters have been collapsed
      uint32_t modifiedGroup = *groupsToCollapse.begin();
      if(pars.debug){
        std::cout << "\tgroupCountsSinceLastReCalc[modifiedGroup]: " << groupCountsSinceLastReCalc[modifiedGroup] << std::endl;
        std::cout << "\tpars.groupDiffToReCalc: " << pars.groupDiffToReCalc << std::endl;
        std::cout << "\tcurrentGroupCounts[modifiedGroup]: " << currentGroupCounts[modifiedGroup] << std::endl;
        std::cout << "\tgroupCountsSinceLastReCalc[modifiedGroup] + pars.groupDiffToReCalc > currentGroupCounts[modifiedGroup]: " << njh::colorBool(groupCountsSinceLastReCalc[modifiedGroup] + pars.groupDiffToReCalc > currentGroupCounts[modifiedGroup]) << std::endl;
        std::cout << "\tcurrentGroupCounts[modifiedGroup] - groupCountsSinceLastReCalc[modifiedGroup] > pars.groupDiffToReCalc: " << njh::colorBool(currentGroupCounts[modifiedGroup] - groupCountsSinceLastReCalc[modifiedGroup] > pars.groupDiffToReCalc) << std::endl;
        std::cout << "\tmodifiedGroup: " << modifiedGroup << std::endl;
      }

      if(currentGroupCounts[modifiedGroup] - groupCountsSinceLastReCalc[modifiedGroup] > pars.groupDiffToReCalc){
        groupCountsSinceLastReCalc = currentGroupCounts;
        njh::concurrent::LockableQueue<uint32_t> groupsQueue(groups);
        uint32_t otherGroup = std::numeric_limits<uint32_t>::max();
        //std::vector<double> allCurrentDists;
        while(groupsQueue.getVal(otherGroup)){
          if(!njh::in(otherGroup, groupsToCollapse)){
            uint32_t group1 = modifiedGroup;
            uint32_t group2 = otherGroup;
            uint32_t group1Size = groupNodes[group1].size();
            uint32_t group2Size = groupNodes[group2].size();
            double sumOfSquaresAll = 0;
            double sumOfSquaresGroup1 = 0;
            double sumOfSquaresGroup2 = 0;

            std::mutex sumsMut;
            //group 1
            if(group1Size > 1){
              PairwisePairFactory group1_pFac(group1Size);
              std::function<void()> computeSumsOfSqaures = [&group1_pFac,&sumsMut,&sumOfSquaresGroup1,&sumOfSquaresAll,&group1,&dist,&groupNodes](){
                PairwisePairFactory::PairwisePairVec group1_pair_vec;
                while(group1_pFac.setNextPairs(group1_pair_vec, 200)){
                  double current_sumOfSquaresAll = 0;
                  for(const auto & group1_pair : group1_pair_vec.pairs_){
                    uint32_t node1 = groupNodes[group1][group1_pair.col_];
                    uint32_t node2 = groupNodes[group1][group1_pair.row_];
                    uint32_t distRow = std::max(node1,  node2);
                    uint32_t distCol = std::min(node1,  node2);
                    double squareDist = std::pow(dist[distRow][distCol], 2.0);
                    current_sumOfSquaresAll += squareDist;
                  }
                  //allCurrentDists.emplace_back(dist[distRow][distCol]);
                  {
                    std::lock_guard<std::mutex> lock(sumsMut);
                    sumOfSquaresGroup1 += current_sumOfSquaresAll;
                    sumOfSquaresAll += current_sumOfSquaresAll;
                  }
                }
              };
              njh::concurrent::runVoidFunctionThreaded(computeSumsOfSqaures, pars.numThreads);
            }
            //group 2
            if(group2Size > 1){
              PairwisePairFactory group2_pFac(group2Size);
              std::function<void()> computeSumsOfSqaures = [&group2_pFac,&sumsMut,&sumOfSquaresGroup2,&sumOfSquaresAll,&group2,&dist,&groupNodes](){
                PairwisePairFactory::PairwisePairVec group2_pair_vec;
                while(group2_pFac.setNextPairs(group2_pair_vec, 200)){
                  double current_sumOfSquaresAll = 0;
                  for(const auto & group2_pair : group2_pair_vec.pairs_){
                    uint32_t node1 = groupNodes[group2][group2_pair.col_];
                    uint32_t node2 = groupNodes[group2][group2_pair.row_];
                    uint32_t distRow = std::max(node1,  node2);
                    uint32_t distCol = std::min(node1,  node2);
                    double squareDist = std::pow(dist[distRow][distCol], 2.0);
                    current_sumOfSquaresAll += squareDist;
                  }
                  //allCurrentDists.emplace_back(dist[distRow][distCol]);
                  {
                    std::lock_guard<std::mutex> lock(sumsMut);
                    sumOfSquaresGroup2 += current_sumOfSquaresAll;
                    sumOfSquaresAll += current_sumOfSquaresAll;
                  }
                }
              };
              njh::concurrent::runVoidFunctionThreaded(computeSumsOfSqaures, pars.numThreads);
            }
            //between group
            AllByAllPairFactory allFac(groupNodes[group1].size(), groupNodes[group2].size());
            std::function<void()> computeSumsOfSqaures = [&allFac,&sumsMut,&sumOfSquaresAll,&dist,&groupNodes,&group1,&group2](){
              AllByAllPairFactory::AllByAllPairVec allPairVec;
              while(allFac.setNextPairs(allPairVec, 200)){
                double current_sumOfSquaresAll = 0;
                for(const auto & allPair : allPairVec.pairs_){
                  auto group1Node = groupNodes[group1][allPair.row_];
                  auto group2Node = groupNodes[group2][allPair.col_];
                  uint32_t distRow = std::max(group1Node,  group2Node);
                  uint32_t distCol = std::min(group1Node,  group2Node);
                  double squareDist = std::pow(dist[distRow][distCol], 2.0);
                  current_sumOfSquaresAll += squareDist;
                }
                //allCurrentDists.emplace_back(dist[distRow][distCol]);
                {
                  std::lock_guard<std::mutex> lock(sumsMut);
                  sumOfSquaresAll += current_sumOfSquaresAll;
                }
              }
            };
            njh::concurrent::runVoidFunctionThreaded(computeSumsOfSqaures, pars.numThreads);
            double distBetweenCentroids = (sumOfSquaresAll - (group1Size + group2Size) * (sumOfSquaresGroup1 / group1Size + sumOfSquaresGroup2 / group2Size)) / (group1Size * group2Size);
            centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)] = distBetweenCentroids;
          }
        }
//        {
//          auto [t, p] = boost::math::statistics::two_sample_t_test(differentInitialGroupDists, allCurrentDists);
//          if(pars.verbose){
//            auto diffDists_mean = vectorMean(differentInitialGroupDists);
//            auto diffDists_sd = vectorStandardDeviationSamp(differentInitialGroupDists);
//            boost::math::normal_distribution diffDistr(diffDists_mean, diffDists_sd);
//            auto allCurrentDists_mean = vectorMean(allCurrentDists);
//            auto allCurrentDists_sd = vectorStandardDeviationSamp(allCurrentDists);
//            boost::math::normal_distribution currentDistr(allCurrentDists_mean, allCurrentDists_sd);
//
//            std::cout << "differentInitialGroupDists.size(): " << differentInitialGroupDists.size() << std::endl;
//            std::cout << "allCurrentDists.size(): " << allCurrentDists.size() << std::endl;
//
//            std::cout << "\t" << vectorMean(differentInitialGroupDists) << "\t" <<  vectorStandardDeviationSamp(differentInitialGroupDists) << std::endl
//                    << "\t" << vectorMean(allCurrentDists) << "\t" <<  vectorStandardDeviationSamp(allCurrentDists) << std::endl;
//            std::cout << "t-statistic: " << t << ", p-value: " << p << ": " << njh::colorBool(p < 0.01) << std::endl;
//            std::cout << "in different overlap: " << cdf(diffDistr,allCurrentDists_mean + allCurrentDists_sd*2) << ": " << njh::colorBool(cdf(diffDistr,allCurrentDists_mean + allCurrentDists_sd*2) < 0.01) << std::endl;
//            std::cout << "in current overlap  : " << 1 - cdf(currentDistr, diffDists_mean - 2 * diffDists_sd) << std::endl;
//          }
//        }
      }

      ret.hdbsRunInfo.addRow(runCount, numberOfNonSingletClusters, numberOfGroups_, lowestCentroidDist);
    }
    //regroup based on counts;
    auto groupCounts = getGroupCounts();
    std::vector<uint32_t> allGroups = getVectorOfMapKeys(groupCounts);
    njh::sort(allGroups,[&groupCounts](uint32_t g1, uint32_t g2){
      return groupCounts[g1] > groupCounts[g2];
    });
    std::unordered_map<uint32_t, uint32_t> reGroupingKey;
    for(const auto & idx : iter::enumerate(allGroups)){
      reGroupingKey[idx.element] = idx.index;
    }
    for(auto & n : nodes_){
      n->group_ = reGroupingKey[n->group_];
    }
    //final centroid
    {
      //compute centroid distances
      std::map<uint32_t, std::vector<uint32_t>> groupNodes;
      for(const auto & node : nodes_){
        groupNodes[node->group_].emplace_back(nameToNodePos_[node->name_]);
      }
      auto groups = getVectorOfMapKeys(groupNodes);
      //
      {
        ret.final_centroidDistances = table(toVecStr(njh::getVecOfMapKeys(groupNodes)));
        for(const auto & group1Idx : iter::range(groups.size())){
          std::vector<DIST> outVec;
          uint32_t group1 = groups[group1Idx];
          for(const auto & group2Idx : iter::range(0UL, group1Idx + 1)){
            uint32_t group2 = groups[group2Idx];
            outVec.emplace_back(centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)] );
          }
          for(const auto & group2Idx : iter::range(group1Idx + 1, groupNodes.size())){
            uint32_t group2 = groups[group2Idx];
            outVec.emplace_back(centroidDistances[std::max(group1,  group2)][std::min(group1,  group2)] );
          }
          ret.final_centroidDistances.addRow(toVecStr(outVec));
        }
      }
    }
    return ret;
  }


}; //class njhUndirWeightedGraph






}  // namespace njhseq
