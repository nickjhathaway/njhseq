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

namespace bibseq {

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
	ret["refPos"] = bib::json::toJson(mis.refBasePos);
	ret["seq"] = bib::json::toJson(comp.queryName_);
	ret["seqBase"] = bib::json::toJson(estd::to_string(mis.seqBase));
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
	auto & nodes = graphJson["nodes"];
	auto & links = graphJson["links"];
	double maxReadCnt = std::numeric_limits<double>::min();
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
	uint32_t lCount = 0;
	for (const auto & e : this->edges_) {
		if (e->on_) {
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


}  // namespace bibseq

