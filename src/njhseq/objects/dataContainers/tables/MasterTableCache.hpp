#pragma once
/*
 * MasterTableCache.hpp
 *
 *  Created on: Feb 7, 2016
 *      Author: nick
 */

//
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

namespace njhseq {
class MasterTableCache {
public:
	MasterTableCache(const TableIOOpts & tabOpts, const std::string & dir,
			bool recursive, uint32_t depth, const std::regex & pattern);

	TableIOOpts tabOpts_;
	const std::string dir_;
	bool recursive_;
	uint32_t depth_;
	const std::regex pattern_;
	std::map<njh::files::bfs::path, bool> files_;
	table masterTab_;

	void collectFiles();
	bool needsUpdate() const;
	void updateTab();
	void writeTab();
};


}  // namespace njhseq

