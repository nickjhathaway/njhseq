#pragma once
/*
 * CutOutputRefPositions.hpp
 *
 *  Created on: Dec 13, 2017
 *      Author: nick
 */

#include "bibseq/alignment/alignerUtils/comparison.hpp"

namespace bibseq {

struct CutOutputRefPositions {
	uint32_t refStart_ = 0;
	uint32_t refStartLength_ = 1;
	uint32_t refStop_ = std::numeric_limits<uint32_t>::max();
	uint32_t refStopLength_ = 1;
	std::string refStartStr_ = "";
	std::string refStopStr_ = "";

	std::string name_ = "";

	bool checkComp_ = false;
	comparison comp_;

	Json::Value toJson() const;

	static std::vector<CutOutputRefPositions> readInPositions(
			const bfs::path & positionsJsonFnp);

	void checkPositionsThrow(const std::string & funcName) const;
	void checkPositionsThrow(const std::string & funcName, const size_t refLength,
			const std::string & refName) const;

	static void setPositions(CutOutputRefPositions & refPositions,
			std::string funcName);

	std::string getId() const;

};


}  // namespace bibseq



