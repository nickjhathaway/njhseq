#pragma once
/*
 * AllGroupDataPaths.hpp
 *
 *  Created on: Sep 13, 2016
 *      Author: nick
 */



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
