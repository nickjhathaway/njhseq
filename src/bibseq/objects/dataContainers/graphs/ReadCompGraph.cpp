/*
 * ReadCompGraph.cpp
 *
 *  Created on: Dec 1, 2015
 *      Author: nick
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


#include "ReadCompGraph.hpp"
#include "bibseq/IO/OutputStream.hpp"
#include "bibseq/concurrency/PairwisePairFactory.hpp"


namespace bibseq {


std::map<uint32_t, std::vector<char>> ReadCompGraph::getVariantSnpLociMap(
		const std::string & name, VecStr names, uint32_t expand) const {
	std::map<uint32_t, std::set<char>> ret;
	auto & n = nodes_[nameToNodePos_.at(name)];
	for (const auto & e : n->edges_) {
		if (bib::in(e->nodeToNode_.at(name).lock()->name_, names)) {
			if (e->dist_.refName_ == name) {
				for (const auto & m : e->dist_.distances_.mismatches_) {
					ret[m.second.refBasePos].insert(m.second.seqBase);
					if (expand > 0) {
						for (const auto pos : iter::range<uint32_t>(1, expand + 1)) {
							if (pos <= m.second.refBasePos) {
								ret[m.second.refBasePos - pos].insert(m.second.seqBase);
							}
							if (pos + m.second.refBasePos < n->value_->seq_.size()) {
								ret[m.second.refBasePos + pos].insert(m.second.seqBase);
							}
						}
					}
				}
			} else {
				for (const auto & m : e->dist_.distances_.mismatches_) {
					ret[m.second.seqBasePos].insert(m.second.refBase);
					if (expand > 0) {
						for (const auto pos : iter::range<uint32_t>(1, expand + 1)) {
							if (pos <= m.second.seqBasePos) {
								ret[m.second.seqBasePos - pos].insert(m.second.refBase);
							}
							if (pos + m.second.seqBasePos < n->value_->seq_.size()) {
								ret[m.second.seqBasePos + pos].insert(m.second.refBase);
							}
						}
					}
				}
			}
		}
	}
	std::map<uint32_t, std::vector<char>> realRet;
	for (const auto & l : ret) {
		realRet[l.first] = std::vector<char> { l.second.begin(), l.second.end() };
	}
	return realRet;
}

std::map<uint32_t, std::vector<gap>> ReadCompGraph::getVariantIndelLociMap(
		const std::string & name, VecStr names, uint32_t expand) const{
	std::map<uint32_t, std::vector<gap>> ret;
	auto & n = nodes_[nameToNodePos_.at(name)];
	for (const auto & e : n->edges_) {
		if (bib::in(e->nodeToNode_.at(name).lock()->name_, names)) {
			if (e->dist_.refName_ == name) {
				for (const auto & g : e->dist_.distances_.alignmentGaps_) {
					auto mainSearch = ret.find(g.second.refPos_);
					if(ret.end() == mainSearch ||
							!bib::contains(mainSearch->second, g.second, [](const gap & g1, const gap & g2){ return g1.compare(g2);})){
						ret[g.second.refPos_].emplace_back(g.second);
					}
					if (expand > 0) {
						for (const auto pos : iter::range<uint32_t>(1, expand + 1)) {
							if (pos <= g.second.refPos_) {
								auto expandPos = g.second.refPos_ - pos;
								auto expandSearch = ret.find(expandPos);
								if (ret.end() == expandSearch
										|| !bib::contains(expandSearch->second, g.second,
												[](const gap & g1, const gap & g2) {return g1.compare(g2);})) {
									ret[expandPos].emplace_back(g.second);
								}
							}
							if (pos + g.second.refPos_ < n->value_->seq_.size()) {
								auto expandPos = g.second.refPos_ + pos;
								auto expandSearch = ret.find(expandPos);
								if (ret.end() == expandSearch
										|| !bib::contains(expandSearch->second, g.second,
												[](const gap & g1, const gap & g2) {return g1.compare(g2);})) {
									ret[expandPos].emplace_back(g.second);
								}
							}
						}
					}
				}
			} else {
				for (const auto & g : e->dist_.distances_.alignmentGaps_) {
					auto gapCopy = g.second;
					gapCopy.switchSeqAndRef();
					auto mainSearch = ret.find(gapCopy.refPos_);
					if(ret.end() == mainSearch ||
							!bib::contains(mainSearch->second, gapCopy, [](const gap & g1, const gap & g2){ return g1.compare(g2);})){
						ret[gapCopy.refPos_].emplace_back(gapCopy);
					}
					if (expand > 0) {
						for (const auto pos : iter::range<uint32_t>(1, expand + 1)) {
							if (pos <= gapCopy.refPos_) {
								auto expandPos = gapCopy.refPos_ - pos;
								auto expandSearch = ret.find(expandPos);
								if (ret.end() == expandSearch
										|| !bib::contains(expandSearch->second, gapCopy,
												[](const gap & g1, const gap & g2) {return g1.compare(g2);})) {
									ret[expandPos].emplace_back(gapCopy);
								}
							}
							if (pos + gapCopy.refPos_ < n->value_->seq_.size()) {
								auto expandPos = gapCopy.refPos_ + pos;
								auto expandSearch = ret.find(expandPos);
								if (ret.end() == expandSearch
										|| !bib::contains(expandSearch->second, gapCopy,
												[](const gap & g1, const gap & g2) {return g1.compare(g2);})) {
									ret[expandPos].emplace_back(gapCopy);
								}
							}
						}
					}
				}
			}
		}
	}
	return ret;
}

comparison ReadCompGraph::setMinimumConnections(
		std::function<void(comparison &)> modFunc,
		std::function<bool(const comparison &, const comparison &)> compFunc) {
	comparison eventsAllowed;
	modFunc(eventsAllowed);
	this->turnOffEdgesWithComp(eventsAllowed, compFunc);
	this->determineGroups();
	while (this->numberOfGroups_ > 1) {
		modFunc(eventsAllowed);
		this->resetVisitedNodes();
		this->resetVisitedEdges();
		this->turnOffEdgesWithComp(eventsAllowed, compFunc);
		this->determineGroups();
	}
	return eventsAllowed;
}

comparison ReadCompGraph::setMinimumEventConnections() {
	auto modFunc =
			[](comparison & comp)->void {++comp.distances_.overLappingEvents_;};
	auto compFunc =
			[](const comparison & observed, const comparison & allowed)->bool {
				return observed.distances_.getNumOfEvents(false) > allowed.distances_.overLappingEvents_;};
	return setMinimumConnections(modFunc, compFunc);
}

comparison ReadCompGraph::setMinimumHqMismatchConnections() {
	auto modFunc = [](comparison & comp)->void {++comp.hqMismatches_;};
	auto compFunc =
			[](const comparison & observed, const comparison & allowed)->bool {
				return observed.hqMismatches_ > allowed.hqMismatches_;};
	return setMinimumConnections(modFunc, compFunc);
}



void ReadCompGraph::setJustBestConnection(bool doTies){
	auto compFunc = [](const comparison & comp, const comparison & best){
		return comp.distances_.getNumOfEvents(false) < best.distances_.getNumOfEvents(false);
	};
	auto equalFunc = [](const comparison & comp, const comparison & best){
		return comp.distances_.getNumOfEvents(false) == best.distances_.getNumOfEvents(false);
	};


	allDetermineBestDistanceWithComp(doTies, compFunc, equalFunc);

}

std::string & padAsNeeded(std::string & str, size_t compareLen) {
	if (compareLen > str.size()) {
		str.append(std::string(compareLen - str.size(), ' '));
	}
	return str;
}

std::string & padAsNeeded(std::string & str, const std::string & compare) {
	return padAsNeeded(str, compare.size());
}

std::string & padAsNeeded(std::string & str, const VecStr & compareStrs) {
	size_t maxCompareLen = 0;
	for (const auto & s : compareStrs) {
		if (s.size() > maxCompareLen) {
			maxCompareLen = s.size();
		}
	}
	return padAsNeeded(str, maxCompareLen);
}

std::string padAsNeededRet(std::string str, size_t compareLen) {
	if (compareLen > str.size()) {
		str.append(std::string(compareLen - str.size(), ' '));
	}
	return str;
}

std::string padAsNeededRet(std::string str, const std::string & compare) {
	return padAsNeededRet(str, compare.size());
}

std::string padAsNeededRet(std::string str, const VecStr & compareStrs) {
	size_t maxCompareLen = 0;
	for (const auto & s : compareStrs) {
		if (s.size() > maxCompareLen) {
			maxCompareLen = s.size();
		}
	}
	return padAsNeededRet(str, maxCompareLen);
}

std::string & padAsNeededHtml(std::string & str, size_t compareLen) {
	if (compareLen > str.size()) {
		str.append(std::string(compareLen - str.size(), ' '));
	}
	return str;
}

std::string & padAsNeededHtml(std::string & str, const std::string & compare) {
	return padAsNeededHtml(str, compare.size());
}

std::string & padAsNeededHtml(std::string & str, const VecStr & compareStrs) {
	size_t maxCompareLen = 0;
	for (const auto & s : compareStrs) {
		if (s.size() > maxCompareLen) {
			maxCompareLen = s.size();
		}
	}
	return padAsNeededHtml(str, maxCompareLen);
}

std::string padAsNeededHtmlRet(std::string str, size_t compareLen) {
	if (compareLen > str.size()) {
		str.append(repeatString("&nbsp;",compareLen - str.size()));
	}
	return str;
}

std::string padAsNeededHtmlRet(std::string str, const std::string & compare) {
	return padAsNeededHtmlRet(str, compare.size());
}

std::string padAsNeededHtmlRet(std::string str, const VecStr & compareStrs) {
	size_t maxCompareLen = 0;
	for (const auto & s : compareStrs) {
		if (s.size() > maxCompareLen) {
			maxCompareLen = s.size();
		}
	}
	return padAsNeededHtmlRet(str, maxCompareLen);
}


Json::Value createSnpNode(const std::string & snpName, uint32_t group,
		const mismatch & mis, const comparison & comp, uint32_t snpNodeSize = 20) {
	Json::Value ret;
	ret["color"] = "red";
	ret["type"] = "snp";
	ret["size"] = snpNodeSize;
	ret["name"] = bib::json::toJson(snpName);
	ret["group"] = bib::json::toJson(group);
	ret["ref"] = bib::json::toJson(comp.refName_);
	ret["refBase"] = bib::json::toJson(estd::to_string(mis.refBase));
	ret["refBaseQual"] = bib::json::toJson(estd::to_string(mis.refQual));
	ret["refPos"] = bib::json::toJson(mis.refBasePos);
	ret["seq"] = bib::json::toJson(comp.queryName_);
	ret["seqBase"] = bib::json::toJson(estd::to_string(mis.seqBase));
	ret["seqBaseQual"] = bib::json::toJson(estd::to_string(mis.seqQual));
	ret["seqPos"] = bib::json::toJson(mis.seqBasePos);
	return ret;
}

Json::Value createIndelNode(const std::string & indelName, uint32_t group,
		const gap & alnGap, const comparison & comp,
		uint32_t indelBaseNodeSize = 20) {
	Json::Value ret;
	ret["color"] = "yellow";
	ret["type"] = "indel";
	ret["size"] = indelBaseNodeSize * alnGap.size_;
	ret["name"] = bib::json::toJson(indelName);
	ret["group"] = bib::json::toJson(group);
	ret["ref"] = bib::json::toJson(comp.refName_);
	ret["seqPos"] = bib::json::toJson(alnGap.seqPos_);
	ret["seq"] = bib::json::toJson(comp.queryName_);
	ret["refPos"] = bib::json::toJson(alnGap.refPos_);
	ret["inRef"] = bib::json::toJson(alnGap.ref_);
	ret["gapSeq"] = bib::json::toJson(alnGap.gapedSequence_);
	//std::string gapSeqFiller = repeatString("&mdash;", alnGap.gapedSequence_.size());
	std::string gapSeqFiller = "";
	std::string refDisplay = "";
	std::string seqDispaly = "";
	if (alnGap.ref_) {
		seqDispaly = alnGap.gapedSequence_;
		refDisplay = gapSeqFiller;
	} else {
		seqDispaly = gapSeqFiller;
		refDisplay = alnGap.gapedSequence_;
	}
	ret["refDisplay"] = bib::json::toJson(refDisplay);
	ret["seqDisplay"] = bib::json::toJson(seqDispaly);

	return ret;
}



Json::Value createVariantNode(const seqInfo & seqBase, uint32_t group,
		double totalReadCount, const scale<double> & cntScale,
		const bib::color & color) {
	Json::Value ret;
	ret["name"] = bib::json::toJson(seqBase.name_);
	ret["frac"] = bib::json::toJson(
			roundDecPlaces(seqBase.cnt_ / totalReadCount, 3));
	ret["cnt"] = bib::json::toJson(seqBase.cnt_);
	ret["group"] = bib::json::toJson(group);
	ret["color"] = "#" + color.hexStr_;
	ret["size"] = cntScale.get(seqBase.cnt_);
	ret["type"] = "variant";
	return ret;
}

Json::Value createVariantNode(const seqInfo & seqBase, uint32_t group,
		double totalReadCount, const scale<double> & cntScale,
		const std::unordered_map<std::string, bib::color> & nameColors) {
	return createVariantNode(seqBase, group, totalReadCount, cntScale, nameColors.at(seqBase.name_));
}

Json::Value createLink(uint32_t source, uint32_t target, bool on,
		const bib::color & color, uint32_t value = 1) {
	Json::Value ret;
	ret["source"] = bib::json::toJson(source);
	ret["target"] = bib::json::toJson(target);
	ret["value"] = bib::json::toJson(value);
	ret["on"] = bib::json::toJson(on);
	ret["color"] = "#" + color.hexStr_;
	return ret;
}



Json::Value ReadCompGraph::toD3Json(bib::color backgroundColor,
		const std::unordered_map<std::string, bib::color> & nameColors) {
	Json::Value graphJson;
	graphJson["backgroundColor"] = "#" + backgroundColor.hexStr_;
	auto & totalDiffs = graphJson["totalDiffs"];
	auto & nodes = graphJson["nodes"];
	auto & links = graphJson["links"];
	double maxReadCnt = std::numeric_limits<double>::lowest();
	double total = 0;
	for (const auto & n : this->nodes_) {
		if (n->on_) {
			if (n->value_->cnt_ > maxReadCnt) {
				maxReadCnt = n->value_->cnt_;
			}
			total += n->value_->cnt_;
		}
	}
	scale<double> cntScale( { 0, maxReadCnt }, { 50.0, 1000.0 });
	uint32_t nCount = 0;
	std::unordered_map<std::string, uint64_t> nameToNewPos;
	uint64_t pos = 0;
	for (const auto & n : this->nodes_) {
		if (n->on_) {
			nameToNewPos[n->name_] = pos;
			++pos;
			nodes[nCount] = createVariantNode(*(n->value_), n->group_, total,
					cntScale, nameColors);
			++nCount;
		}
	}
	//std::unordered_map<std::string, uint32_t> totalDiffsCounter;
	uint32_t lCount = 0;
	for (const auto & e : this->edges_) {
		if (e->on_) {
			totalDiffs[e->dist_.refName_ + e->dist_.queryName_] = e->dist_.distances_.getNumOfEvents(false);
			if (e->dist_.distances_.getNumOfEvents(false) == 0) {
				auto lColor =
						getColsBetweenExcludeClosest(
								nameColors.at(
										this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_),
								nameColors.at(e->nodeToNode_.begin()->second.lock()->name_), 1);

				links[lCount] =
						createLink(
								nameToNewPos[this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_],
								nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_],
								e->on_, lColor.front());
				++lCount;
			} else {

				std::string lastName =
						this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_;
				auto lColors =
						getColsBetweenExcludeClosest(
								nameColors.at(
										this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_),
								nameColors.at(e->nodeToNode_.begin()->second.lock()->name_),
								e->dist_.distances_.getNumOfEvents(false) + 1);
				uint32_t evenNumber = 0;
				for (const auto & snp : e->dist_.distances_.mismatches_) {
					std::string snpName;
					if (e->dist_.refName_
							== e->nodeToNode_.begin()->second.lock()->name_) {
						snpName = vectorToString(
								toVecStr(snp.second.seqBasePos, snp.second.seqBase,
										snp.second.refBasePos, snp.second.refBase), "");
					} else {
						snpName = vectorToString(
								toVecStr(snp.second.refBasePos, snp.second.refBase,
										snp.second.seqBasePos, snp.second.seqBase), "");
					}
					std::string newName =
							this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_
									+ "_" + snpName + "_"
									+ e->nodeToNode_.begin()->second.lock()->name_;
					nameToNewPos[newName] = pos;
					++pos;
					nodes[nCount] = createSnpNode(newName,
							e->nodeToNode_.begin()->second.lock()->group_, snp.second,
							e->dist_);
					++nCount;
					links[lCount] = createLink(nameToNewPos[lastName],
							nameToNewPos[newName], true, lColors[evenNumber]);
					++lCount;
					lastName = newName;
					++evenNumber;
				}
				for (const auto & indel : e->dist_.distances_.alignmentGaps_) {
					std::string indelName;
					if (e->dist_.refName_
							== e->nodeToNode_.begin()->second.lock()->name_) {
						indelName = vectorToString(
								toVecStr(indel.second.seqPos_, indel.second.gapedSequence_,
										indel.second.refPos_), "");
					} else {
						indelName = vectorToString(
								toVecStr(indel.second.refPos_, indel.second.gapedSequence_,
										indel.second.seqPos_), "");
					}
					std::string newName =
							this->nodes_[this->nameToNodePos_[e->nodeToNode_.begin()->first]]->name_
									+ "_" + indelName + "_"
									+ e->nodeToNode_.begin()->second.lock()->name_;
					nameToNewPos[newName] = pos;
					++pos;
					nodes[nCount] = createIndelNode(newName,
							e->nodeToNode_.begin()->second.lock()->group_, indel.second,
							e->dist_);
					++nCount;
					links[lCount] = createLink(nameToNewPos[lastName],
							nameToNewPos[newName], true, lColors[evenNumber]);
					++lCount;
					lastName = newName;
					++evenNumber;
				}
				links[lCount] = createLink(nameToNewPos[lastName],
						nameToNewPos[e->nodeToNode_.begin()->second.lock()->name_], true,
						lColors[e->dist_.distances_.getNumOfEvents(false)]);
				++lCount;
			}
		}
	}
	return graphJson;
}



Json::Value ReadCompGraph::getSingleLineJsonOut(const ConnectedHaplotypeNetworkPars & netPars,
		const std::unordered_map<std::string, bib::color> & colorLookup){
	Json::Value outputJson;
		auto & nodes = outputJson["nodes"];
		auto & links = outputJson["links"];
		double minNodeCntObserved = std::numeric_limits<double>::max();
		double maxNodeCntObserved = std::numeric_limits<double>::lowest();
		for(const auto & n : nodes_){
			if(n->on_){
				if(n->value_->cnt_ > maxNodeCntObserved){
					maxNodeCntObserved = n->value_->cnt_;
				}
				if(n->value_->cnt_ < minNodeCntObserved){
					minNodeCntObserved = n->value_->cnt_;
				}
			}
		}

		scale<double> nodeSizeScale(
				std::make_pair(1, std::max(2.0, maxNodeCntObserved)),
				std::make_pair<double>(78.53982, 78.53982 * 20));

		for(const auto & n : nodes_){
			if(n->on_){
				Json::Value node;
				node["id"] = n->name_;
				if(!netPars.noLabel){
					if(nullptr != netPars.seqMeta && "" != netPars.labelField){
						auto labForName = netPars.seqMeta->groupData_.at(netPars.labelField)->getGroupForSample(n->name_);
						node["label"] = "NA" ==  labForName ? std::string("") : labForName ;
					}else{
						node["label"] = n->name_;
					}
				}else{
					node["label"] = "";
				}
				if(nullptr != netPars.seqMeta && "" != netPars.colorField){
					node["color"] = bib::mapAt(colorLookup, netPars.seqMeta->groupData_.at(netPars.colorField)->getGroupForSample(n->name_)).getHexStr();
				}else{
					node["color"] = "#00A0FA";
				}
				node["size"] = bib::json::toJson(nodeSizeScale.get(n->value_->cnt_));
				nodes.append(node);
			}
		}
		double minIdObserved = std::numeric_limits<double>::max();
		double maxIdObserved = std::numeric_limits<double>::lowest();
		for(const auto & e : edges_){
			if(netPars.setJustBest && !e->best_){
				continue;
			}
			if(e->on_){
				if(e->dist_.distances_.eventBasedIdentity_  < minIdObserved){
					minIdObserved = e->dist_.distances_.eventBasedIdentity_ ;
				}
				if(e->dist_.distances_.eventBasedIdentity_  > maxIdObserved){
					maxIdObserved = e->dist_.distances_.eventBasedIdentity_ ;
				}
			}
		}
		scale<double> linkSizeScale(std::make_pair(minIdObserved, maxIdObserved), std::make_pair<double>(1,20));
		for (const auto & e : edges_) {
			if (netPars.setJustBest && !e->best_) {
				continue;
			}
			if(e->on_){
				Json::Value link;
				link["source"] = e->nodeToNode_.begin()->first ;
				link["target"] = e->nodeToNode_.begin()->second.lock()->name_ ;
				link["value"] = linkSizeScale.get(e->dist_.distances_.eventBasedIdentity_);
				links.append(link);
			}
		}
		return outputJson;
}

void ReadCompGraph::writeAdjListPerId(const OutOptions & linksOutOpts,
		const ConnectedHaplotypeNetworkPars & netPars){
	OutputStream linksOut(linksOutOpts);
	for (const auto & e : edges_) {
		if (netPars.setJustBest && !e->best_) {
			continue;
		}
		if(e->on_){
			linksOut << e->nodeToNode_.begin()->first << " -- "
					<< e->nodeToNode_.begin()->second.lock()->name_
					<< " " << e->dist_.distances_.getNumOfEvents(true)
					<< " " << e->dist_.distances_.eventBasedIdentity_ << std::endl;
		}
	}
}

void ReadCompGraph::addEdgesBasedOnIdOrMinDif(const ConnectedHaplotypeNetworkPars & netPars,
		concurrent::AlignerPool & alnPool){

	std::vector<kmerInfo> kInfos;
	for(const auto & n : nodes_){
		kInfos.emplace_back(n->value_->seq_, netPars.matchPars.kmerLen_, false	);
	}
	PairwisePairFactory pFactory(nodes_.size());

	std::mutex graphMut;
	bib::ProgressBar pBar(pFactory.totalCompares_);
	auto addEdges = [this,&pFactory,&kInfos,&alnPool,&netPars,&graphMut,&pBar](){
		auto alignerObj = alnPool.popAligner();
		PairwisePairFactory::PairwisePairVec pairs;
		while(pFactory.setNextPairs(pairs, netPars.matchPars.batchAmount_)){
			std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<comparison>>> currentComps;
			for(const auto & pair : pairs.pairs_){
				if(nodes_[pair.col_]->value_->seq_ == nodes_[pair.row_]->value_->seq_){
					alignerObj->alignObjectA_ = baseReadObject(*nodes_[pair.row_]->value_);
					alignerObj->alignObjectB_ = baseReadObject(*nodes_[pair.col_]->value_);
				}else{
					if(kInfos[pair.row_].compareKmers(kInfos[pair.col_]).second < netPars.matchPars.kmerCutOff_){
						continue;
					}
					alignerObj->alignCacheGlobal(nodes_[pair.row_]->value_, nodes_[pair.col_]->value_);
				}
				alignerObj->profileAlignment(nodes_[pair.row_]->value_, nodes_[pair.col_]->value_, false, false, false);
				currentComps[nodes_[pair.row_]->value_->name_][nodes_[pair.col_]->value_->name_] = std::make_shared<comparison>(alignerObj->comp_);
			}
			{
				std::lock_guard<std::mutex> lock(graphMut);
				for(const auto & r1 : currentComps){
					for(const auto & r2 : r1.second){
						addEdge(r1.first, r2.first, *r2.second);
					}
				}
				if(netPars.verbose){
					pBar.outputProgAdd(std::cout, pairs.pairs_.size(), true);
				}
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < netPars.numThreads; ++t){
		threads.emplace_back(addEdges);
	}
	bib::concurrent::joinAllJoinableThreads(threads);

	if(netPars.setJustBest){
		//turn of same seq connections as so the best connection will not be just equal matches
		for(auto & e : edges_){
			if(e->nodeToNode_.begin()->second.lock()->value_->seq_ == e->nodeToNode_.rbegin()->second.lock()->value_->seq_){
				e->on_ = false;
			}
		}
		setJustBestConnection(netPars.doTies);
		//turn exact matches back on
		for(auto & e : edges_){
			if(e->nodeToNode_.begin()->second.lock()->value_->seq_ == e->nodeToNode_.rbegin()->second.lock()->value_->seq_){
				e->on_ = true;
				e->best_ = true;
			}
		}
	}

	if(std::numeric_limits<double>::lowest() != netPars.minId){
		for(auto & e : edges_){
			if(e->on_ && e->dist_.distances_.eventBasedIdentity_ < netPars.minId){
				e->on_ = false;
			}
		}
	}
	if(std::numeric_limits<uint32_t>::max() != netPars.minNumberOfEvents){
		for(auto & e : edges_){
			if(e->on_ && e->dist_.distances_.getNumOfEvents(true) > netPars.minNumberOfEvents){
				e->on_ = false;
			}
		}
	}
}


std::string ReadCompGraph::ConnectedHaplotypeNetworkPars::htmlPageForConnectHpaNet = R"FOOBAR(<!DOCTYPE html>
		<meta charset="utf-8">
		<svg id="chart" width="2500" height="2500"></svg>
		<script src="https://d3js.org/d3.v4.min.js"></script>
		<script
		  src="http://code.jquery.com/jquery-3.3.1.min.js"
		  integrity="sha256-FgpCb/KJQlLNfOu91ta32o/NMZxltwRo8QtmkMRdAu8="
		  crossorigin="anonymous"></script>
		<script>


		var svg = d3.select("svg"),
		    width = +svg.attr("width"),
		    height = +svg.attr("height");
		    
		function addSvgSaveButton(buttonId, topSvg, name) {
		    	name = name || "graph.svg"
		    	if(!name.endsWith('.svg')){
		    		name += ".svg";
		    	}
		    	d3.select(buttonId).append("a").attr("id", "imgDownload");
		    	d3.select(buttonId).on(
		    			"click",
		    			function() {
		    				var html = $(
		    						d3.select(topSvg).attr("version", 1.1).attr("xmlns",
		    								"http://www.w3.org/2000/svg").node()).clone()
		    						.wrap('<p/>').parent().html();
		    				// add the svg information to a and then click it to trigger the
		    				// download
		    				var imgsrc = 'data:image/svg+xml;base64,' + btoa(html);
		    				d3.select("#imgDownload").attr("download", name);
		    				d3.select("#imgDownload").attr("href", imgsrc);
		    				var a = $("#imgDownload")[0];
		    				a.click();
		    			});
		    }
		d3.select("body").append("button")
		    					.style("float", "top")
		    					.attr("class", "btn btn-success")
		    					.attr("id", "saveButton")
		    					.style("margin", "2px")
		    					.text("Save As Svg");
		addSvgSaveButton("#saveButton", "#chart", "connectedNetwork")
		//var color = d3.scaleOrdinal(d3.schemeCategory20);

		var simulation = d3.forceSimulation()
		    .force("link", d3.forceLink().id(function(d) { return d.id; }))
		    .force("charge", d3.forceManyBody().strength(-120))
			    .force("center", d3.forceCenter(width / 2, height / 2))
			    .force("y", d3.forceY(height / 2))
			    .force("x", d3.forceX(width / 2)) ;
		    //.force("charge", d3.forceManyBody())
		    //.force("center", d3.forceCenter(width / 2, height / 2));

		d3.json("tree.json", function(error, graph) {
		  if (error) throw error;

		  var initializerSimulation = d3.forceSimulation(graph.nodes)
		          .force("link", d3.forceLink(graph.links).id(function(d) { return d.id; }))
		          .force("charge", d3.forceManyBody().strength(-120))
		            .force("center", d3.forceCenter(width / 2, height / 2))
		            .force("y", d3.forceY(height / 2))
		            .force("x", d3.forceX(width / 2))
		            .stop(); ;
		  var n = Math.ceil(Math.log(initializerSimulation.alphaMin()) / Math.log(1 - initializerSimulation.alphaDecay()));
		  //n = n * 2
		  for (var i = 0; i < n; ++i) {
		    //postMessage({type: "tick", progress: i / n});
		    console.log(i,"", n);
		    initializerSimulation.tick();
		  }

		  var link = svg.append("g")
		      .attr("class", "links")
		    .selectAll("line")
		    .data(graph.links)
		    .enter().append("line")
		      .attr("stroke-width", function(d) { return Math.sqrt(d.value); })
		      .attr("stroke", "#999")
		      .attr("stroke-opacity", 0.6);


		  svg.selectAll(".node")
		    		.data(graph.nodes)
		    		.enter().append("g")
		    			.attr("class", function(d) {return "node";})
		    				.call(d3.drag()
		    			          .on("start", dragstarted)
		    			          .on("drag", dragged)
		    			          .on("end", dragended))
		    			.append("circle");
		  var node = svg.selectAll(".node");
		  node.select("circle")
		    		  .attr("r", function(d) { return Math.sqrt(d.size/Math.PI); })
		    	    .style("fill", function(d) { return d.color; })
		    	    .style("stroke", "#fff")
		    	    .style("stroke-width", "1.5px").append("title")
		              .text(function(d) { return d.id; });


		  /*
		  var node = svg.append("g")
		      .attr("class", "nodes")
		    .selectAll("circle")
		    .data(graph.nodes)
		    .enter().append("circle")
		      .attr("r", function(d) { return Math.sqrt(d.size/Math.PI); })
		      .attr("fill", function(d) { return d.color; })
		      .call(d3.drag()
		          .on("start", dragstarted)
		          .on("drag", dragged)
		          .on("end", dragended));
		          node.append("title")
		              .text(function(d) { return d.id; });
		          */

		  

		  simulation
		      .nodes(graph.nodes)
		      .on("tick", ticked);

		  simulation.force("link")
		      .links(graph.links);

		  node.append("text")
		      	  .attr("x", 12)
		      	  .attr("dy", ".35em")
		      	  .style("fill","#000")
		      	  .style("font-family", "\"HelveticaNeue-Light\", \"Helvetica Neue Light\", \"Helvetica Neue\", Helvetica, Arial, \"Lucida Grande\", sans-serif")
		      	  .style("font-size", "12px")
		      	  .style("font-weight","900")
		      	  .style("pointer-events", "none")
		      	  .style("stroke", "#000")
		      	  .style("stroke-width", "1px")
		      	  .text(function(d) {return d.label;});

		    function ticked() {
		  	    link.attr("x1", function(d) { return d.source.x; })
		  	        .attr("y1", function(d) { return d.source.y; })
		  	        .attr("x2", function(d) { return d.target.x; })
		  	        .attr("y2", function(d) { return d.target.y; });
		  	    node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")" ; });
		  	};
		  	//for drag objects start
		  	function dragstarted(d) {
		  	  if (!d3.event.active) simulation.alphaTarget(0.3).restart();
		  	  d.fx = d.x;
		  	  d.fy = d.y;
		  	};
		  	//for when dragging has been detected
		  	function dragged(d) {
		  	  d.fx = d3.event.x;
		  	  d.fy = d3.event.y;
		  	};
		  	//for when dragging has ended
		  	function dragended(d) {
		  	  if (!d3.event.active) simulation.alphaTarget(0);
		  	  d.fx = null;
		  	  d.fy = null;
		  	};
		    });
		  /*
		  function ticked() {
		    link
		        .attr("x1", function(d) { return d.source.x; })
		        .attr("y1", function(d) { return d.source.y; })
		        .attr("x2", function(d) { return d.target.x; })
		        .attr("y2", function(d) { return d.target.y; });

		    node
		        .attr("cx", function(d) { return d.x; })
		        .attr("cy", function(d) { return d.y; });
		  }
		});

		function dragstarted(d) {
		  if (!d3.event.active) simulation.alphaTarget(0.3).restart();
		  d.fx = d.x;
		  d.fy = d.y;
		}

		function dragged(d) {
		  d.fx = d3.event.x;
		  d.fy = d3.event.y;
		}

		function dragended(d) {
		  if (!d3.event.active) simulation.alphaTarget(0);
		  d.fx = null;
		  d.fy = null;
		}*/

		</script>


		)FOOBAR";

}  // namespace bibseq

