#pragma once
//
//  urlUtils.hpp
//
//  Created by Nick Hathaway on 5/27/15.

//
//
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
#include "bibseq/common.h"
#include "bibseq/utils/bitSwaps.hpp"

namespace bibseq {



bool inline xdigit(int c) {
	//from cppcms 1.05
	return ('0' <= c && c <= '9') || ('a' <= c && c <= 'f')
			|| ('A' <= c && c <= 'F');
}

std::string urldecode(char const *begin, char const *end);

std::string urldecode(std::string const &s);


template<typename Iterator>
void urlencode_impl(char const *b, char const *e, Iterator out) {
	//from cppcms 1.05
	while (b != e) {
		char c = *b++;
		if (('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z')
				|| ('0' <= c && c <= '9')) {
			*out++ = c;
		} else {
			switch (c) {
			case '-':
			case '_':
			case '.':
			case '~':
				*out++ = c;
				break;
			default: {
				static char const hex[] = "0123456789abcdef";
				unsigned char uc = c;
				*out++ = '%';
				*out++ = hex[(uc >> 4) & 0xF];
				*out++ = hex[uc & 0xF];

			}
			};
		}
	};
}

void urlencode(char const *b, char const *e, std::ostream &out);
std::string urlencode(const std::string &s);

}  // namespace bibseq

