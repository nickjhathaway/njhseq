#pragma once
/*
 * TableReader.hpp
 *
 *  Created on: Jan 21, 2018
 *      Author: nick
 */

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




