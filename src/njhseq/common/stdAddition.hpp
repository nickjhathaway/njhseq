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
//
//  stdAddition.hpp
//
//  Created by Nicholas Hathaway on 1/8/14.
//

#include "allSystemIncludes.h"
#include <njhcpp.h>

namespace njhseq {

template<typename T>
typename T::size_type len(const T & con){
	return con.size();
}

template<typename T>
const T & getRef(const T & val) {
	return val;
}

template<typename T>
const T & getRef(const std::shared_ptr<const T> & val) {
	return getRef(*val);
}

template<typename T>
const T & getRef(const std::unique_ptr<const T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(T & val) {
	return val;
}

template<typename T>
T & getRef(const std::shared_ptr<T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(const std::unique_ptr<T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(std::shared_ptr<T> & val) {
	return getRef(*val);
}

template<typename T>
T & getRef(std::unique_ptr<T> & val) {
	return getRef(*val);
}


}  // njh



