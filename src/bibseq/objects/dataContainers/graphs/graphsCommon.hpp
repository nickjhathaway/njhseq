#pragma once
/*
 * graphsCommon.hpp
 *
 *  Created on: Dec 16, 2016
 *      Author: nick
 */

#include <bibcpp/jsonUtils.h>
#include <bibcpp/graphics.h>

#include "bibseq/common.h"
#include "bibseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "bibseq/seqToolsUtils/distCalc.hpp"

namespace bibseq {

std::vector<bib::color> getColsBetweenExcludeClosest(bib::color first,
		bib::color last, uint32_t num);

std::unordered_map<std::string, bib::color> getColorsForNames(
		const VecStr & popNames);
void jsonTreeToDot(Json::Value treeData, std::ostream & outDot);

void genTreeHtml(std::ostream & out, const std::string & jsonFileName,
		const std::string & treeJsFilename);

void genSimpleTreeJs(std::ostream & out);



}  // namespace bibseq



