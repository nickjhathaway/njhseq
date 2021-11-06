#pragma once

/*
 * TableWriter.hpp
 *
 *  Created on: Nov 6, 2021
 *      Author: nick
 */


#include "njhseq/objects/dataContainers/tables/table.hpp"
#include "njhseq/IO/IOOptions.h"
#include "njhseq/IO/OutputStream.hpp"
#include "njhseq/IO/InputStream.hpp"

namespace njhseq {

class TableWriter {
public:

	TableWriter(const TableIOOpts & tabOpts, const VecStr & header);
	const TableIOOpts tabOpts_;
	table header_;

	std::unique_ptr<OutputStream> out_;

	void writeRow(VecStr & row);
	void writeRowLock(VecStr & row);

	std::mutex mut_;
	bool doNotCheckRowSizes = false;
};

}  // namespace njhseq


