#pragma once
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
/*
 * GFFCore.hpp
 *
 *  Created on: Dec 28, 2015
 *      Author: nick
 */

#include "njhseq/common.h"
#include "njhseq/utils.h"


namespace njhseq {

class GFFCore {
public:
	GFFCore();
	GFFCore(const std::string & line);

	std::string seqid_;
	std::string source_;
	std::string type_;
	size_t start_; //counting from 1 as per the GFF file specifications (ugh)
	size_t end_; //counting from 1 as per the GFF file specifications (ugh)
	double score_;
	char strand_; //either + or - (. for NA)
	uint16_t phase_; //ether 0, 1, or 2 (uint16_t max for NA/.)
	std::unordered_map<std::string, std::string> attributes_;

	bool hasAttr(const std::string & attr) const;
	template<typename T>
	T getAttrAs(const std::string & attr) const {
		if(hasAttr(attr)){
			return njh::lexical_cast<T>(attributes_.at(attr));
		}
		return T{};
	}
	std::string getAttr(const std::string & attr) const;

	bool isReverseStrand() const;

	std::string getIDAttr() const;

	Json::Value toJson() const;

	void writeGffRecord(std::ostream & out) const;

};



}  // namespace njhseq


