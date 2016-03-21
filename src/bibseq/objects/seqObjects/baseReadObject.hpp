#pragma once
//
//  baseReadObject.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 8/31/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
//
// bibseq - A library for analyzing sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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
#include "bibseq/utils.h"
#include "bibseq/objects/seqObjects/seqInfo.hpp"
namespace bibseq {

class baseReadObject {

public:
	baseReadObject();
	baseReadObject(const seqInfo& seqBase);

	seqInfo seqBase_;
	using size_type = seqInfo::size_type;
	virtual Json::Value toJson() const;

	virtual ~baseReadObject(){}
};

template<>
inline baseReadObject::size_type len(const baseReadObject & read){
	return read.seqBase_.seq_.size();
}

}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "baseReadObject.cpp"
#endif
