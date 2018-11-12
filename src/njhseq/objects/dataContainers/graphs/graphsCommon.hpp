#pragma once
/*
 * graphsCommon.hpp
 *
 *  Created on: Dec 16, 2016
 *      Author: nick
 */

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

#include <njhcpp/jsonUtils.h>
#include <njhcpp/graphics.h>

#include "njhseq/common.h"
#include "njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp"
#include "njhseq/seqToolsUtils/distCalc.hpp"

namespace njhseq {

std::vector<njh::color> getColsBetweenExcludeClosest(njh::color first,
		njh::color last, uint32_t num);

std::unordered_map<std::string, njh::color> getColorsForNames(
		const VecStr & popNames);
void jsonTreeToDot(Json::Value treeData, std::ostream & outDot);

void genTreeHtml(std::ostream & out, const std::string & jsonFileName,
		const std::string & treeJsFilename);

void genSimpleTreeJs(std::ostream & out);



}  // namespace njhseq



