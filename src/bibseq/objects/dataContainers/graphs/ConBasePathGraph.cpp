/*
 * ConBasePathGraph.cpp
 *
 *  Created on: Jan 14, 2017
 *      Author: nick
 */


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
			linkJsons["value"] = tl->cnt_ / totalTail;
			links.append(linkJsons);
		}
	}
	return ret;
}


}  // namespace bibseq

