/*
 * ConBasePathGraph.cpp
 *
 *  Created on: Jan 14, 2017
 *      Author: nick
 */
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "ConBasePathGraph.hpp"

namespace bibseq {

std::string ConBasePathGraph::ConPath::PosBase::getUid() const {
	return estd::to_string(pos_) + estd::to_string(base_);
}

Json::Value ConBasePathGraph::ConPath::PosBase::toJson() const {
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["pos_"] = bib::json::toJson(pos_);
	ret["base_"] = bib::json::toJson(base_);
	return ret;
}


ConBasePathGraph::ConPath::ConPath(double count):count_(count){}

void ConBasePathGraph::ConPath::addPosBase(uint32_t pos, char base){
	addPosBase({pos, base});
}
void ConBasePathGraph::ConPath::addPosBase(const PosBase & pb){
	bases_.emplace_back(pb);
}

std::string ConBasePathGraph::ConPath::getUid()const{
	std::string ret = "";
	for(const auto & pb : bases_){
		ret += pb.getUid();
	}
	return ret;
}

Json::Value ConBasePathGraph::ConPath::toJson() const {
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["count_"] = bib::json::toJson(count_);
	ret["bases_"] = bib::json::toJson(bases_);
	return ret;
}


ConBasePathGraph::node::node(const ConBasePathGraph::ConPath::PosBase & val,
		double cnt, double frac) :
		val_(val),
		cnt_(cnt),
		frac_(frac){
}

Json::Value ConBasePathGraph::node::toJson() const {
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["val_"] = bib::json::toJson(val_);
	ret["cnt_"] = bib::json::toJson(cnt_);
	ret["frac_"] = bib::json::toJson(frac_);
	ret["headEdges_"] = bib::json::toJson(headEdges_);
	ret["tailEdges_"] = bib::json::toJson(tailEdges_);
	ret["visitCount_"] = bib::json::toJson(visitCount_);
	return ret;
}

Json::Value ConBasePathGraph::edge::toJson() const {
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["head_"] = bib::json::toJson(head_.lock()->val_.getUid());
	ret["tail_"] = bib::json::toJson(tail_.lock()->val_.getUid());
	ret["cnt_"] = bib::json::toJson(cnt_);
	return ret;
}

void ConBasePathGraph::node::resetVisitCount(){
	visitCount_ = 0;
}

void ConBasePathGraph::node::addHead(const std::shared_ptr<edge> & e) {
	headEdges_.push_back(e);
}

void ConBasePathGraph::node::addTail(const std::shared_ptr<edge> & e) {
	tailEdges_.push_back(e);
}

bool ConBasePathGraph::node::headless() const {
	return headEdges_.empty();
}

bool ConBasePathGraph::node::tailless() const {
	return tailEdges_.empty();
}

void ConBasePathGraph::node::addToWritingPath(std::ostream & out,
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

void ConBasePathGraph::node::addToPath(std::vector<ConPath> & paths,
		ConPath currentPath) {
/*
 * 	if(visitCount_ > 0){
		std::cout << "paths" << std::endl;
		for(const auto & p : paths){
			std::cout << p.getUid() << std::endl;
		}
		std::cout << "current path: " << std::endl;
		std::cout << currentPath.getUid() << std::endl;
		throw std::runtime_error{"visited before"};
	}
 */
	/*
	std::cout << "paths" << std::endl;
	for (const auto & p : paths) {
		std::cout << p.getUid() << std::endl;
	}*/
	++visitCount_;
	/*
	if (currentPath.bases_.size() > 0) {
		if (val_.pos_ <= currentPath.bases_.back().pos_ ) {
			std::cout << "paths" << std::endl;
			for (const auto & p : paths) {
				std::cout << p.getUid() << std::endl;
			}
			std::cout << "current path: " << std::endl;
			std::cout << currentPath.getUid() << std::endl;
			throw std::runtime_error { "smaller pos" };
		}
	}
	*/
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

ConBasePathGraph::edge::edge(const std::shared_ptr<node> & head,
		const std::shared_ptr<node> & tail,
		double cnt) :
		head_(head), tail_(tail),
		cnt_(cnt){

};


void ConBasePathGraph::addNode(const ConPath::PosBase & n, double cnt, double frac) {
	if (bib::has(nodes_, n.getUid())) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, already contains node with name "
				<< n.getUid() << "\n";
		throw std::runtime_error { ss.str() };
	}
	nodes_.emplace(n.getUid(), std::make_shared<node>(n, cnt, frac));
}

void ConBasePathGraph::addEdge(const std::string & head,
		const std::string & tailName,
		double cnt) {
	//get head node
	//see if tail node already exists and if it does, just add count to the edge
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

void ConBasePathGraph::writePaths(std::ostream & out) const {
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

std::vector<ConBasePathGraph::ConPath> ConBasePathGraph::getPaths() const {

	std::vector<ConPath> ret;
	/**@todo add way to check if there are any cycles or no headless nodes */
	for (const auto & n : nodes_) {
		if (n.second->headless()) {
			//std::cout << __FILE__ << " : " << __LINE__  << " : " << __PRETTY_FUNCTION__ << std::endl;
			//std::cout << n.second->val_.toJson() << std::endl;
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

Json::Value ConBasePathGraph::createSankeyOutput() const {
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
		nodeJson["pos"] = n->val_.pos_;
		nodeJson["base"] = n->val_.base_;
		nodeJson["cnt"] = n->cnt_;
		nodeJson["frac"] = n->frac_;
		nodes.append(nodeJson);
		double totalTail = 0;
		for (const auto & tl : n->tailEdges_) {
			totalTail += tl->cnt_;
		}
		for (const auto & tl : n->tailEdges_) {
			Json::Value linkJsons;
			//linkJsons["source"] = tl->head_.lock()->val_.getUid();
			//linkJsons["target"] = tl->tail_.lock()->val_.getUid();
			linkJsons["source"] = nodePosition[tl->head_.lock()->val_.getUid()];
			linkJsons["target"] = nodePosition[tl->tail_.lock()->val_.getUid()];
			linkJsons["value"] = tl->cnt_ / totalTail;
			links.append(linkJsons);
		}
	}
	return ret;
}


}  // namespace bibseq

