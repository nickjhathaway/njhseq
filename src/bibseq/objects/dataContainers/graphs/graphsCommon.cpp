/*
 * graphsCommon.cpp
 *
 *  Created on: Dec 16, 2016
 *      Author: nick
 */


#include "graphsCommon.hpp"

namespace bibseq {

std::vector<bib::color> getColsBetweenExcludeClosest(bib::color first,
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


std::unordered_map<std::string,bib::color> getColorsForNames(const VecStr & popNames){
	auto popColors = bib::njhColors(popNames.size());
	bibseq::VecStr popColorsStrs(popColors.size(), "");
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
	outDot << "#overlap = false; " << std::endl;
	outDot << "fixedsize = true; " << std::endl;
	bgColor.invert();
	outDot << "fontcolor = \"" << bgColor.getHexStr() <<"\"" << std::endl;
	outDot << "fontsize = 20" << std::endl;
	outDot << "fontname = \"helvetica\"" << std::endl;
	std::vector<double> sizes;
	for (const auto & node : treeData["nodes"]){
		if(node["color"].asString() != "red"){
			sizes.emplace_back(node["size"].asDouble());
		}
	}

	scale<double> cntScale({0, vectorMaximum(sizes)},{0.5, 5});
	//std::cout << cntScale.toJson() << std::endl;
	uint32_t mCounts = 0;
	std::unordered_map<std::string, std::string> decoder;
	for (const auto & node : treeData["nodes"]){
		if(node["color"].asString() == "red" || node["color"].asString() == "yellow"){
			// width = 0.15 , label =""
			decoder[node["name"].asString() ] = "mis" + estd::to_string(mCounts);
			outDot << "mis" << mCounts << " [shape=circle,style=filled,fixedsize "
			                                    "=true, color = \"#000000\", fillcolor ="
			          << "\"" << node["color"].asString() << "\""
			          << ", width = 0.15, label =\"\"]"
			          << std::endl;
			++mCounts;
		}else{
			//dot doesn't allow for "." in the names, it splits the name into two different things
			decoder[node["name"].asString() ] = bib::replaceString(node["name"].asString(), ".","_");
			outDot << bib::replaceString(node["name"].asString(), ".","_") << " [shape=circle,style=filled,fixedsize "
			                                    "=true, color = \"#000000\", fillcolor ="
			          << "\"" << node["color"].asString() << "\""
			          << ", width = " << cntScale.get(node["size"].asDouble())<< "]"
			          << std::endl;
		}
	}
	for(const auto & link : treeData["links"]){
		outDot << decoder[treeData["nodes"][link["source"].asUInt()]["name"].asString()]
				<< " -- " << decoder[treeData["nodes"][link["target"].asUInt()]["name"].asString()]
          << " [penwidth=5, color=\""
          << link["color"].asString() << "\"]" << std::endl;
	}
	outDot << "}" << std::endl;
}

void genTreeHtml(std::ostream & out, const std::string & jsonFileName,
		const std::string & treeJsFilename) {
	out << "<!DOCTYPE html>" << std::endl;
	out << "<meta charset=\"utf-8\">" << std::endl;
	out << "<body style=\"background-color: #000\">" << std::endl;
	out << "<script src=\"http://d3js.org/d3.v3.min.js\">" << std::endl;
	out << "</script>" << std::endl;
	out
			<< "<script src=\"http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js\">"
			<< std::endl;
	out << "</script>" << std::endl;
	out << "<script src=\"" << treeJsFilename << "\">" << std::endl;
	out << "</script>" << std::endl;
	out << "<button id=\"save_svg\">" << std::endl;
	out << "Save as Svg</button>" << std::endl;
	out << "<div style = \"width: 2500px; margin:0 auto;text-align: center;\">"
			<< std::endl;
	out << "<svg id = \"main\">" << std::endl;
	out << "</svg>" << std::endl;
	out << "</div>" << std::endl;
	out << "<script>" << std::endl;
	out << "var jsonDat = \"" << jsonFileName << "\";" << std::endl;
	out << "var add = \"#main\";drawPsuedoMinTree(jsonDat, add);" << std::endl;
	out << "var bName = \"#save_svg\";addSvgSaveButton(bName, add);" << std::endl;
	out << "</script>" << std::endl;
	out << "</body>" << std::endl;
}

void genSimpleTreeJs(std::ostream & out) {
		out << "function drawPsuedoMinTree(jsonData, addTo){" << std::endl;
		out << "	var width = 2500," << std::endl;
		out << "    height = 2500;" << std::endl;
		out << "" << std::endl;
		out << "var force = d3.layout.force()" << std::endl;
		out << "    .charge(-120)" << std::endl;
		out << "    .linkDistance(30)" << std::endl;
		out << "    .size([width, height]);" << std::endl;
		out << "	d3.select(addTo)" << std::endl;
		out << "    	.attr(\"width\", width)" << std::endl;
		out << "    	.attr(\"height\", height);" << std::endl;
		out << "var svg = d3.select(addTo).append(\"svg\")" << std::endl;
		out << "    .attr(\"width\", width)" << std::endl;
		out << "    .attr(\"height\", height)" << std::endl;
		out << "    .attr(\"id\", \"chart\");" << std::endl;
		out << "" << std::endl;
		out << "d3.json(jsonData, function(error, graph) {" << std::endl;
		out << "  svg.append('rect')" << std::endl;
		out << "	    .attr('width',width)" << std::endl;
		out << "	    .attr('height', height)" << std::endl;
		out << "	    .attr('fill',graph.backgroundColor);" << std::endl;
		out << "  " << std::endl;
		out << "  force" << std::endl;
		out << "      .nodes(graph.nodes)" << std::endl;
		out << "      .links(graph.links)" << std::endl;
		out << "      .on(\"tick\",tick)" << std::endl;
		out << "      .start();" << std::endl;
		out << "" << std::endl;
		out << "  var link = svg.selectAll(\".link\")" << std::endl;
		out << "      .data(graph.links)" << std::endl;
		out << "    .enter().append(\"line\")" << std::endl;
		out << "      .attr(\"class\", \"link\")" << std::endl;
		out << "      .style(\"stroke\", function(d) { return d.color; })" << std::endl;
		out << "      .style(\"stroke-width\", function(d) { return 2; });" << std::endl;
		out << "var node = svg.selectAll(\".node\")" << std::endl;
		out << "    .data(force.nodes())" << std::endl;
		out << "  .enter().append(\"g\")" << std::endl;
		out << "    .attr(\"class\", \"node\")" << std::endl;
		out << "    .call(force.drag);" << std::endl;
		out << "" << std::endl;
		out << "// add the nodes" << std::endl;
		out << "	node.append(\"circle\")" << std::endl;
		out << "      .attr(\"r\", function(d) { return Math.sqrt(d.size/Math.PI); })" << std::endl;
		out << "      .style(\"fill\", function(d) { return d.color; })" << std::endl;
		out << "      .style(\"stroke\", function(d) { return d.color; })" << std::endl;
		out << "      .style(\"stroke-width\", \"1.5px\")" << std::endl;
		out << "      .call(force.drag);" << std::endl;
		out << "" << std::endl;
		out << "  node.append(\"title\")" << std::endl;
		out << "      .text(function(d) { return d.name; });" << std::endl;
		out << "" << std::endl;
		out << "function tick() {" << std::endl;
		out << "    link.attr(\"x1\", function(d) { return d.source.x; })" << std::endl;
		out << "        .attr(\"y1\", function(d) { return d.source.y; })" << std::endl;
		out << "        .attr(\"x2\", function(d) { return d.target.x; })" << std::endl;
		out << "        .attr(\"y2\", function(d) { return d.target.y; });" << std::endl;
		out << "" << std::endl;
		out << "    node.attr(\"transform\", function(d) { return \"translate(\" + d.x + \",\" + d.y + \")\"; });" << std::endl;
		out << "};" << std::endl;
		out << "" << std::endl;
		out << "});" << std::endl;
		out << "}" << std::endl;
		out << "" << std::endl;
		out << "function addSvgSaveButton(buttonId, topSvg) {" << std::endl;
		out << "	d3.select(buttonId).append(\"a\").attr(\"id\", \"imgDownload\");" << std::endl;
		out << "	d3.select(buttonId).on(" << std::endl;
		out << "			\"click\"," << std::endl;
		out << "			function() {" << std::endl;
		out << "				var html = $(" << std::endl;
		out << "						d3.select(topSvg).attr(\"version\", 1.1).attr(\"xmlns\"," << std::endl;
		out << "								\"http://www.w3.org/2000/svg\").node()).clone()" << std::endl;
		out << "						.wrap('<p/>').parent().html();" << std::endl;
		out << "				;" << std::endl;
		out << "				// add the svg information to a and then click it to trigger the" << std::endl;
		out << "				// download" << std::endl;
		out << "				var imgsrc = 'data:image/svg+xml;base64,' + btoa(html);" << std::endl;
		out << "				d3.select(\"#imgDownload\").attr(\"download\", \"graph.svg\");" << std::endl;
		out << "				d3.select(\"#imgDownload\").attr(\"href\", imgsrc);" << std::endl;
		out << "				var a = $(\"#imgDownload\")[0];" << std::endl;
		out << "				a.click();" << std::endl;
		out << "			});" << std::endl;
		out << "}" << std::endl;
}


}  // namespace bibseq


