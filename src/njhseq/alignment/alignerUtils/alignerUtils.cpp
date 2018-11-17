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
 * alignerUtils.cpp
 *
 *  Created on: Nov 30, 2015
 *      Author: nick
 */

#include "alignerUtils.hpp"

namespace njhseq {



size_t getAlnPosForRealPos(const std::string & seq, size_t realSeqPos) {
	if (0 == seq.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": seq size should not be zero" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (realSeqPos >= seq.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": realSeqPos is greater than or equal to seq size, seq size: "
				<< seq.size() << ", realSeqPos: " << realSeqPos << std::endl;
		throw std::runtime_error { ss.str() };
	}
	size_t bases = 0;
	for (size_t i = 0; i < seq.size(); ++i) {
		if ('-' != seq[i]) {
			++bases;
		}
		if (bases == realSeqPos + 1) {
			return i;
		}
	}
	if (realSeqPos >= bases) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": realSeqPos is greater than or equal to the number of bases counted size, bases: "
				<< bases << ", realSeqPos: " << realSeqPos << std::endl;
		throw std::runtime_error { ss.str() };
	}
	return seq.size() - 1;
}

size_t getRealPosForAlnPos(const std::string & seq, size_t seqAlnPos) {
	if (0 == seq.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": seq size should not be zero" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if (seqAlnPos >= seq.size()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ": seqAlnPos is greater than or equal to seq size, seq size: "
				<< seq.size() << ", seqAlnPos: " << seqAlnPos << std::endl;
		throw std::runtime_error { ss.str() };
	}
	if(seq.find_first_not_of('-') >=seqAlnPos){
		return 0;
	}
	size_t realPos = 0;
	for (size_t i = 0; i <= seqAlnPos; ++i) {
		if ('-' != seq[i]) {
			++realPos;
		}
	}
	//to account for further gaps at this position
	if ('-' == seq[seqAlnPos] && 0 != realPos) {
		++realPos;
	}
	if (realPos == 0) {
		return std::numeric_limits<size_t>::max();
	}
	return realPos - 1;
}

size_t getRealPosForSeqForRefPos(const std::string & ref,
		const std::string & seq, size_t refPos) {
	return getRealPosForAlnPos(seq, getAlnPosForRealPos(ref, refPos));
}

size_t getRealPosForRefForSeqPos(const std::string & ref,
		const std::string & seq, size_t seqPos) {
	return getRealPosForAlnPos(ref, getAlnPosForRealPos(seq, seqPos));
}


}  // namespace njhseq

