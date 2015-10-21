//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
/*
 * UndirWeightedGraph.cpp
 *
 *  Created on: Dec 10, 2014
 *      Author: nickhathaway
 */

#include "UndirWeightedGraph.hpp"

namespace bibseq {
std::unordered_map<std::string,bib::color> getColorsForNames(const VecStr & popNames){
	auto popColors = bib::njhColors(popNames.size());
	bibseq::VecStr popColorsStrs(popColors.size());
	uint32_t count = 0;
	uint32_t halfCount = 0;
	for(const auto & cPos : iter::range(popColors.size())) {
		uint32_t pos = 0;
		if(cPos %2 == 0) {
			pos = popColors.size()/2 + halfCount;
			++halfCount;
		} else {
			pos = count;
			++count;
		}
		popColorsStrs[cPos] = "#" + popColors[pos].hexStr_;
	}
	std::unordered_map<std::string,bib::color> nameColors;
	for(auto pos : iter::range(popNames.size())){
		nameColors[popNames[pos]] = popColorsStrs[pos];
	}
	return nameColors;
}

void jsonTreeToDot(Json::Value treeData, std::ostream & outDot){
	bib::color bgColor = bib::color(treeData["backgroundColor"].asString());
	outDot << "graph G  { " << std::endl;
	outDot << "bgcolor =\"" << bgColor.getHexStr() <<"\"" << std::endl;
	outDot << "overlap = false; " << std::endl;
	outDot << "fixedsize = true; " << std::endl;
	bgColor.invert();
	outDot << "fontcolor = \"" << bgColor.getHexStr() <<"\"" << std::endl;
	outDot << "fontsize = 20" << std::endl;
	std::vector<double> sizes;
	for (const auto & node : treeData["nodes"]){
		if(node["color"].asString() != "red"){
			sizes.emplace_back(node["size"].asDouble());
		}
	}

	scale<double> cntScale({0, vectorMaximum(sizes)},{0.5, 3});
	//std::cout << cntScale.toJson() << std::endl;
	for (const auto & node : treeData["nodes"]){
		if(node["color"].asString() == "red"){
			// width = 0.15 , label =""
			outDot << node["name"].asString() << " [shape=circle,style=filled,fixesize "
			                                    "=true, color = \"#000000\", fillcolor ="
			          << "\"" << node["color"].asString() << "\""
			          << ", width = 0.15, label =\"\"]"
			          << std::endl;
		}else{
			outDot << node["name"].asString() << " [shape=circle,style=filled,fixesize "
			                                    "=true, color = \"#000000\", fillcolor ="
			          << "\"" << node["color"].asString() << "\""
			          << ", width = " << cntScale.get(node["size"].asDouble())<< "]"
			          << std::endl;
		}
	}
	for(const auto & link : treeData["links"]){
		outDot << treeData["nodes"][link["source"].asUInt()]["name"].asString()
				<< " -- " << treeData["nodes"][link["target"].asUInt()]["name"].asString()
          << " [penwidth=5, color=\""
          << link["color"].asString() << "\"]" << std::endl;
	}
	outDot << "}" << std::endl;
}

} /* namespace bibseq */
