#include "baseReadObject.hpp"
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
namespace njhseq {

baseReadObject::baseReadObject() :
		seqBase_(seqInfo()) {
}

baseReadObject::baseReadObject(const seqInfo& seqBase) :
		seqBase_(seqBase) {
}

Json::Value baseReadObject::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson("njhseq::baseReadObject");
	ret["seqBase_"] = njh::json::toJson(seqBase_);
	return ret;
}

baseReadObject::~baseReadObject() {
}

}  // namespace njhseq
