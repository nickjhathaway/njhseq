#include "baseReadObject.hpp"
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
namespace bibseq {

baseReadObject::baseReadObject() : seqBase_(seqInfo()) {}


baseReadObject::baseReadObject(const seqInfo& seqBase) : seqBase_(seqBase) {
  //std::cout << "baseReadObject constructor: " << std::endl;
  //std::cout << seqBase_.name_ << std::endl;
  //std::cout << seqBase_.cnt_ << std::endl;
  //std::cout << seqBase_.frac_ << std::endl;
}

Json::Value baseReadObject::toJson()const{
	Json::Value ret;
	ret["class"] = bib::json::toJson("bibseq::baseReadObject");
	ret["seqBase_"] = bib::json::toJson(seqBase_);
	return ret;
}
}  // namespace bib
