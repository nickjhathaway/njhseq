#pragma once
/*
 * AllGroupDataPaths.hpp
 *
 *  Created on: Sep 13, 2016
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


#include "bibseq/objects/Meta/MultipleGroupMetaData.hpp"



namespace bibseq {

class AllGroupDataPaths {
public:
	class GroupDataPaths {
	public:
		class SubGroupDataPaths {
		public:
			SubGroupDataPaths(const bfs::path & mainDir);

			bfs::path mainDir_;
			bfs::path popFileFnp_;
			bfs::path sampFileFnp_;
			bfs::path hapIdTabFnp_;
			bfs::path subGroupNamesDataFnp_;

			VecStr readInPopUIDs() const;
			VecStr readInSampNames() const;

			Json::Value toJson() const;
		};

		GroupDataPaths(const bfs::path & mainDir,
				const std::set<std::string> & subGroups);

		bfs::path mainDir_;
		bfs::path groupInfoFnp_;
		std::unordered_map<std::string, SubGroupDataPaths> groupPaths_;

		Json::Value toJson() const;

	};
	AllGroupDataPaths(const bfs::path & mainDir,
			const std::unique_ptr<MultipleGroupMetaData> & groupMetaData);

	bfs::path mainDir_;
	std::unordered_map<std::string, GroupDataPaths> allGroupPaths_;


	Json::Value toJson() const;

};


}  // namespace bibseq
