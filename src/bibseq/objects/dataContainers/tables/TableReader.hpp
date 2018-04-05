#pragma once
/*
 * TableReader.hpp
 *
 *  Created on: Jan 21, 2018
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
#include "bibseq/objects/dataContainers/tables/table.hpp"
#include "bibseq/IO/IOOptions.h"
#include "bibseq/IO/OutputStream.hpp"
#include "bibseq/IO/InputStream.hpp"

namespace bibseq {

class TableReader {
public:

	TableReader(const TableIOOpts & tabOpts);
	const TableIOOpts tabOpts_;
	table header_;

	std::unique_ptr<InputStream> in_;

	bool getNextRow(VecStr & row);

	VecStr extractCols(const VecStr & row, const VecStr & cols) const;

};

}  // namespace bibseq




