#pragma once
/*
 * TableReader.hpp
 *
 *  Created on: Jan 21, 2018
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
#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/IO/IOOptions.h"
#include "njhseq/IO/OutputStream.hpp"
#include "njhseq/IO/InputStream.hpp"

namespace njhseq {

class TableReader {
public:

	TableReader(const TableIOOpts & tabOpts);
	const TableIOOpts tabOpts_;
	table header_;

	std::unique_ptr<InputStream> in_;

	bool getNextRow(VecStr & row);

	VecStr extractCols(const VecStr & row, const VecStr & cols) const;

	void setHeaderlessHeader(uint32_t numOfCols) ;

	bool doNotCheckRowSizes = false;
};

}  // namespace njhseq




